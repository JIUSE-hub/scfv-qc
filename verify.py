# -*- coding: utf-8 -*-
"""
verify.py — pAIM1 scFv QC 웹도구 자체 검사
==============================================================================
    python verify.py

검사별 통과/실패를 표로 출력하고, 하나라도 실패하면 종료코드 1 을 반환합니다.

  [A] 구문   core.py / xlsx_writer.py / verify.py 의 AST, index.html 의 JS 구문,
             % 포맷 자리표시자 대 인자 개수, 함수 안의 미정의 이름
  [B] 계약   index.html 이 참조하는 컬럼명·키가 core.py 가 실제로 주는 것인가
  [C] 누출   index.html 에 서열이 하드코딩되어 있지 않은가
  [D] 회귀   testdata/ 로 core.analyze 를 돌려 알려진 값이 재현되는가
  [E] 단위   합성 서열로 core.check_landmark 의 5 개 상태를 직접 확인하고,
             index.html 의 badge() 가 그 5 종을 모두 처리하는지 대조

표준 라이브러리만 사용합니다.
node 와 openpyxl 은 있으면 쓰고 없으면 해당 검사를 "건너뜀" 으로 표시합니다.
기대값은 별도 파일이 아니라 이 파일 안의 EXPECT 리터럴에 둡니다
(.gitignore 가 *.csv 를 막으므로).
"""

import ast
import builtins
import io
import os
import re
import shutil
import string
import subprocess
import sys
import tempfile
import traceback
import unicodedata

VERIFY_VERSION = "1.3"

ROOT = os.path.dirname(os.path.abspath(__file__))
PY_FILES = ("core.py", "xlsx_writer.py", "verify.py")
HTML_FILE = "index.html"
TESTDATA = "testdata"

# --- [B] 이슈 3 이 아직 열려 있는 동안 허용되는 DESIGN_DOC 누락 -----------------
# 이슈 3 을 해결하면 이 집합을 비워야 하고, 그러면 검사도 0 건을 요구하게 됩니다.
KNOWN_DESIGN_DOC_GAP = {"rna_bone_marrow", "rna_peripheral"}

# --- [C] 서열이 아니라 문자 클래스이므로 제외 --------------------------------
MARKSEQ_REGEX = r"/\b[ACGT]{8,}\b/g"

# --- [D] 기대값 --------------------------------------------------------------
EXPECT = {
    "param_hash": "3473927a",
    "clones": [
        {"id": "c01", "verdict": "NO_LINKER",  "insert": 276, "mix": 2.5},
        {"id": "c02", "verdict": "PASS",       "insert": 728, "d1": 366, "d2": 309, "mix": 2.1},
        {"id": "c03", "verdict": "FRAMESHIFT", "insert": 735, "d1": 346, "d2": 336, "mix": 3.9},
        {"id": "c04", "verdict": "PASS",       "insert": 746, "d1": 372, "d2": 321, "mix": 1.9},
    ],
    "cdr3_median": 10,
    "vl": [("IGKV4", 1), ("IGKV3", 1)],
    "jh": [("JH1|JH2", 1), ("JH4", 1)],
    "fasta_n": 2,
}

# JS 내장 프로퍼티 — 데이터 키가 아니므로 하위키 대조에서 제외합니다.
JS_PROPS = {
    "length", "forEach", "map", "filter", "find", "join", "split", "slice",
    "replace", "indexOf", "push", "pop", "concat", "sort", "reverse", "some",
    "every", "includes", "startsWith", "endsWith", "toFixed", "toString",
    "trim", "keys", "values", "entries", "has", "add", "delete", "get", "set",
    "size", "then", "catch", "flat", "reduce", "charAt", "substring", "padStart",
}


# =============================================================================
#  결과 수집과 표 출력
# =============================================================================
PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"
_MARK = {PASS: "통과", FAIL: "실패", SKIP: "건너뜀"}


class Report(object):
    def __init__(self):
        self.rows = []

    def add(self, section, name, status, detail=""):
        self.rows.append((section, name, status, detail))

    def ok(self, section, name, cond, detail_ok="", detail_ng=""):
        if cond:
            self.add(section, name, PASS, detail_ok)
        else:
            self.add(section, name, FAIL, detail_ng or detail_ok)

    def count(self, section, status):
        return sum(1 for r in self.rows if r[0] == section and r[2] == status)

    def failed(self):
        return any(r[2] == FAIL for r in self.rows)


def guard(report, section, label, fn, *args):
    """검사 하나가 예외로 죽어도 표를 끝까지 낸다. 예외는 그 자체로 실패."""
    try:
        fn(*args)
    except Exception as e:
        tb = traceback.extract_tb(sys.exc_info()[2])
        where = ""
        if tb:
            last = tb[-1]
            where = " @ %s:%d" % (os.path.basename(last.filename), last.lineno)
        report.add(section, label, FAIL,
                   "검사 중 예외 %s: %s%s" % (type(e).__name__, e, where))


def width(text):
    """한글을 2칸으로 세는 표시폭."""
    n = 0
    for ch in str(text):
        n += 2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1
    return n


def pad(text, w):
    text = str(text)
    return text + " " * max(0, w - width(text))


def clip(text, w):
    text = str(text).replace("\n", " ")
    if width(text) <= w:
        return text
    out = ""
    for ch in text:
        if width(out) + width(ch) > w - 1:
            return out + "…"
        out += ch
    return out


def print_table(report):
    head = ("구분", "검사", "결과", "상세")
    detail_w = 62
    body = [(r[0], r[1], _MARK[r[2]], clip(r[3], detail_w)) for r in report.rows]
    w0 = max([width(head[0])] + [width(r[0]) for r in body])
    w1 = max([width(head[1])] + [width(r[1]) for r in body])
    w2 = max([width(head[2])] + [width(r[2]) for r in body])
    w3 = max([width(head[3])] + [width(r[3]) for r in body])
    line = "  ".join(("-" * w0, "-" * w1, "-" * w2, "-" * w3))
    print(line)
    print("  ".join((pad(head[0], w0), pad(head[1], w1), pad(head[2], w2), head[3])))
    print(line)
    for r in body:
        print("  ".join((pad(r[0], w0), pad(r[1], w1), pad(r[2], w2), r[3])))
    print(line)


# =============================================================================
#  파일 읽기
# =============================================================================
def read_text(name):
    with io.open(os.path.join(ROOT, name), encoding="utf-8") as fh:
        return fh.read()


def parse_py(name):
    return ast.parse(read_text(name), filename=name)


# =============================================================================
#  index.html 조각 내기
# =============================================================================
_SCRIPT_RE = re.compile(r"<script\b([^>]*)>(.*?)</script>", re.S | re.I)


def script_blocks(html):
    """외부 src 가 아닌 <script> 블록의 본문 목록."""
    out = []
    for m in _SCRIPT_RE.finditer(html):
        if "src=" in m.group(1).lower():
            continue
        out.append(m.group(2))
    return out


def glue_source(js):
    """index.html 안의 GLUE 템플릿 리터럴(파이썬 접착부) 본문."""
    m = re.search(r"const\s+GLUE\s*=\s*`", js)
    if not m:
        return None
    start = m.end()
    end = js.find("`", start)
    if end < 0:
        return None
    return js[start:end]


_JS_STR_RE = re.compile(r"'((?:[^'\\\n]|\\.)*)'|\"((?:[^\"\\\n]|\\.)*)\"")


def js_strings(js):
    out = []
    for m in _JS_STR_RE.finditer(js):
        out.append(m.group(1) if m.group(1) is not None else m.group(2))
    return out


def js_array_literal(js, decl):
    """`const NAME = [...]` / `const NAME = new Set([...])` 의 문자열 항목."""
    m = re.search(re.escape(decl) + r"\s*=\s*(?:new\s+Set\s*\()?\s*\[", js)
    if not m:
        return None
    start = js.index("[", m.end() - 1)
    depth = 0
    for j in range(start, len(js)):
        if js[j] == "[":
            depth += 1
        elif js[j] == "]":
            depth -= 1
            if depth == 0:
                return js_strings(js[start:j + 1])
    return None


def js_function_body(js, name):
    """function NAME(...) { ... } 의 본문. 못 찾으면 None."""
    m = re.search(r"function\s+" + re.escape(name) + r"\s*\([^)]*\)\s*\{", js)
    if not m:
        return None
    start = m.end() - 1
    depth = 0
    for j in range(start, len(js)):
        if js[j] == "{":
            depth += 1
        elif js[j] == "}":
            depth -= 1
            if depth == 0:
                return js[start:j + 1]
    return None


def js_member_refs(js, obj):
    """`obj.key` 와 `obj.key.sub` 참조를 (key, sub) 쌍으로 모은다."""
    pat = re.compile(r"(?<![A-Za-z0-9_$.])" + re.escape(obj) +
                     r"\.([A-Za-z_$][A-Za-z0-9_$]*)(?:\.([A-Za-z_$][A-Za-z0-9_$]*))?")
    out = []
    for m in pat.finditer(js):
        out.append((m.group(1), m.group(2)))
    return out


# =============================================================================
#  [A] 구문
# =============================================================================
def check_ast(report):
    for name in PY_FILES:
        try:
            parse_py(name)
            report.add("A", name + " AST 파싱", PASS, "구문 오류 없음")
        except SyntaxError as e:
            report.add("A", name + " AST 파싱", FAIL,
                       "%s:%s %s" % (name, e.lineno, e.msg))


def node_error_line(stderr, tmpdir):
    """node --check 출력에서 사람이 읽을 한 줄만 뽑는다 (임시경로 제거)."""
    lines = [x.strip() for x in (stderr or "").splitlines() if x.strip()]
    pick = ""
    for x in lines:
        if re.match(r"^\w*Error\b", x):
            pick = x
            break
    if not pick:
        pick = lines[0] if lines else "node --check 실패"
    for x in lines:
        m = re.search(r"\.js:(\d+)", x)
        if m:
            pick = "블록 %s행 · %s" % (m.group(1), pick)
            break
    return pick.replace(tmpdir, "").replace(os.sep + "block", "block")


def check_js_syntax(report, html):
    node = shutil.which("node")
    blocks = script_blocks(html)
    if not blocks:
        report.add("A", "index.html JS 구문", FAIL, "<script> 블록을 찾지 못했습니다")
        return
    if not node:
        report.add("A", "index.html JS 구문", SKIP,
                   "node 없음 · <script> %d 블록 미검사" % len(blocks))
        return
    tmp = tempfile.mkdtemp(prefix="scfvqc_verify_")
    try:
        bad = []
        for i, src in enumerate(blocks):
            path = os.path.join(tmp, "block%d.js" % (i + 1))
            with io.open(path, "w", encoding="utf-8") as fh:
                fh.write(src)
            p = subprocess.run([node, "--check", path], capture_output=True,
                               text=True, encoding="utf-8", errors="replace")
            if p.returncode != 0:
                bad.append("블록%d: %s" % (i + 1, node_error_line(p.stderr, tmp)))
        report.ok("A", "index.html JS 구문", not bad,
                  "node --check 통과 · <script> %d 블록" % len(blocks),
                  " / ".join(bad))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


_PCT_RE = re.compile(
    r"%(%|[-+ #0]*(?:\*|\d+)?(?:\.(?:\*|\d+))?[hlL]?[diouxXeEfFgGcrsa])")

# 우변이 이 노드면 값 개수를 정적으로 알 수 없으므로 검사하지 않습니다.
UNKNOWN_RHS = (ast.Name, ast.Call, ast.Subscript, ast.Attribute,
               ast.GeneratorExp, ast.Dict, ast.DictComp, ast.Starred)


def pct_placeholders(text):
    """%% 를 제외한 자리표시자 개수. 개수를 알 수 없으면 None."""
    if "%(" in text:
        return None
    n = 0
    for m in _PCT_RE.finditer(text):
        if m.group(1) == "%":
            continue
        if "*" in m.group(1):
            return None
        n += 1
    return n


def rhs_arity(node):
    """% 우변이 공급하는 값의 개수. 알 수 없으면 None."""
    if isinstance(node, ast.Tuple):
        if any(isinstance(e, ast.Starred) for e in node.elts):
            return None
        return len(node.elts)
    if isinstance(node, UNKNOWN_RHS):
        return None
    return 1


def format_field_arity(text):
    """.format 의 위치 자리표시자 개수. 이름/속성 참조가 섞이면 None."""
    auto, indexed = 0, -1
    for _lit, field, _spec, _conv in string.Formatter().parse(text):
        if field is None:
            continue
        head = field.split(".")[0].split("[")[0]
        if head == "":
            auto += 1
        elif head.isdigit():
            indexed = max(indexed, int(head))
        else:
            return None
    if auto and indexed >= 0:
        return None
    return auto if auto else (indexed + 1)


def check_format_arity(report):
    bad, checked, skipped = [], 0, 0
    for name in PY_FILES:
        tree = parse_py(name)
        for node in ast.walk(tree):
            if (isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mod)
                    and isinstance(node.left, ast.Constant)
                    and isinstance(node.left.value, str)):
                need = pct_placeholders(node.left.value)
                have = rhs_arity(node.right)
                if need is None or have is None:
                    skipped += 1
                    continue
                checked += 1
                if need != have:
                    bad.append("%s:%d %% 자리 %d 대 인자 %d"
                               % (name, node.lineno, need, have))
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "format"
                    and isinstance(node.func.value, ast.Constant)
                    and isinstance(node.func.value.value, str)):
                if node.keywords or any(isinstance(a, ast.Starred) for a in node.args):
                    skipped += 1
                    continue
                need = format_field_arity(node.func.value.value)
                if need is None:
                    skipped += 1
                    continue
                checked += 1
                if need != len(node.args):
                    bad.append("%s:%d .format 자리 %d 대 인자 %d"
                               % (name, node.lineno, need, len(node.args)))
    report.ok("A", "포맷 자리표시자 대조", not bad,
              "대조 %d 건 · 개수 불명으로 제외 %d 건" % (checked, skipped),
              " / ".join(bad[:4]) + (" 외 %d" % (len(bad) - 4) if len(bad) > 4 else ""))


def _bind_targets(node, out):
    if isinstance(node, ast.Name):
        out.add(node.id)
    elif isinstance(node, (ast.Tuple, ast.List)):
        for e in node.elts:
            _bind_targets(e, out)
    elif isinstance(node, ast.Starred):
        _bind_targets(node.value, out)


def arg_names(args):
    out = set()
    for a in list(args.posonlyargs) + list(args.args) + list(args.kwonlyargs):
        out.add(a.arg)
    if args.vararg:
        out.add(args.vararg.arg)
    if args.kwarg:
        out.add(args.kwarg.arg)
    return out


def bound_names(body):
    """한 스코프가 묶는 이름 전부. 중첩 함수/람다 안으로는 들어가지 않는다."""
    out = set()
    stack = list(body)
    while stack:
        n = stack.pop()
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out.add(n.name)
            continue
        if isinstance(n, ast.Lambda):
            continue
        if isinstance(n, ast.Assign):
            for t in n.targets:
                _bind_targets(t, out)
        elif isinstance(n, (ast.AugAssign, ast.AnnAssign, ast.NamedExpr)):
            _bind_targets(n.target, out)
        elif isinstance(n, (ast.For, ast.AsyncFor)):
            _bind_targets(n.target, out)
        elif isinstance(n, (ast.With, ast.AsyncWith)):
            for it in n.items:
                if it.optional_vars is not None:
                    _bind_targets(it.optional_vars, out)
        elif isinstance(n, ast.ExceptHandler):
            if n.name:
                out.add(n.name)
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            for a in n.names:
                out.add((a.asname or a.name).split(".")[0])
        elif isinstance(n, (ast.Global, ast.Nonlocal)):
            out.update(n.names)
        elif isinstance(n, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)):
            for g in n.generators:
                _bind_targets(g.target, out)
        stack.extend(ast.iter_child_nodes(n))
    return out


def _scan_names(node, scope, undef, modname):
    """Load 로 쓰이는 이름이 스코프에 없으면 undef 에 담는다."""
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        for d in node.decorator_list:
            _scan_names(d, scope, undef, modname)
        for d in list(node.args.defaults) + [x for x in node.args.kw_defaults if x]:
            _scan_names(d, scope, undef, modname)
        # 중첩함수는 바깥 스코프를 그대로 물려받는다 (클로저)
        inner = set(scope) | arg_names(node.args) | bound_names(node.body)
        for st in node.body:
            _scan_names(st, inner, undef, modname)
        return
    if isinstance(node, ast.Lambda):
        for d in list(node.args.defaults) + [x for x in node.args.kw_defaults if x]:
            _scan_names(d, scope, undef, modname)
        _scan_names(node.body, set(scope) | arg_names(node.args), undef, modname)
        return
    if isinstance(node, ast.ClassDef):
        for d in node.decorator_list + list(node.bases):
            _scan_names(d, scope, undef, modname)
        inner = set(scope) | bound_names(node.body)
        for st in node.body:
            _scan_names(st, inner, undef, modname)
        return
    if isinstance(node, (ast.ListComp, ast.SetComp, ast.GeneratorExp, ast.DictComp)):
        inner = set(scope)
        for g in node.generators:
            _bind_targets(g.target, inner)
        for g in node.generators:
            _scan_names(g.iter, inner, undef, modname)
            for c in g.ifs:
                _scan_names(c, inner, undef, modname)
        if isinstance(node, ast.DictComp):
            _scan_names(node.key, inner, undef, modname)
            _scan_names(node.value, inner, undef, modname)
        else:
            _scan_names(node.elt, inner, undef, modname)
        return
    if isinstance(node, ast.ExceptHandler):
        inner = set(scope)
        if node.name:
            inner.add(node.name)
        if node.type is not None:
            _scan_names(node.type, inner, undef, modname)
        for st in node.body:
            _scan_names(st, inner, undef, modname)
        return
    if isinstance(node, ast.Name):
        if isinstance(node.ctx, ast.Load) and node.id not in scope:
            undef.append("%s:%d %s" % (modname, node.lineno, node.id))
        return
    for c in ast.iter_child_nodes(node):
        _scan_names(c, scope, undef, modname)


def check_undefined_names(report):
    base = set(dir(builtins)) | {"__name__", "__file__", "__doc__", "__spec__"}
    undef = []
    for name in PY_FILES:
        tree = parse_py(name)
        scope = base | bound_names(tree.body)
        for st in tree.body:
            _scan_names(st, scope, undef, name)
    report.ok("A", "함수 안 미정의 이름", not undef,
              "검사 %d 파일 · 미정의 0 건" % len(PY_FILES),
              " / ".join(undef[:5]) +
              (" 외 %d" % (len(undef) - 5) if len(undef) > 5 else ""))


# =============================================================================
#  [B] 계약
# =============================================================================
def _dict_literal_keys(d):
    out = {}
    for k, v in zip(d.keys, d.values):
        if isinstance(k, ast.Constant) and isinstance(k.value, str):
            out[k.value] = v
    return out


def _returned_dict(rv):
    if isinstance(rv, ast.Dict):
        return rv
    if isinstance(rv, ast.Call) and rv.args and isinstance(rv.args[0], ast.Dict):
        return rv.args[0]      # json.dumps({...}) 형태
    return None


def provided_keys(fn):
    """함수가 돌려주는 dict 의 키 집합과, 값이 dict 리터럴인 키의 하위키 맵."""
    keys = {}
    names = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and node.value is not None:
            d = _returned_dict(node.value)
            if d is not None:
                keys.update(_dict_literal_keys(d))
            elif isinstance(node.value, ast.Name):
                names.add(node.value.id)
    for node in ast.walk(fn):
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id in names and isinstance(node.value, ast.Dict):
                keys.update(_dict_literal_keys(node.value))
            if (isinstance(t, ast.Subscript) and isinstance(t.value, ast.Name)
                    and t.value.id in names and isinstance(t.slice, ast.Constant)
                    and isinstance(t.slice.value, str)):
                keys.setdefault(t.slice.value, node.value)
    sub = dict((k, set(_dict_literal_keys(v)))
               for k, v in keys.items() if isinstance(v, ast.Dict))
    return set(keys), sub


def find_function(tree, name):
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def check_summary_headers(report, js, core):
    headers = set(core.SUMMARY_HEADERS)
    total, missing = 0, []
    for decl in ("const COMPACT", "const NUMCOL", "const MONOCOL", "const BADGECOL"):
        items = js_array_literal(js, decl)
        label = decl.split()[-1]
        if items is None:
            report.add("B", label + " 대 SUMMARY_HEADERS", FAIL,
                       "index.html 에서 " + decl + " 선언을 찾지 못했습니다")
            continue
        total += len(items)
        gone = [x for x in items if x not in headers]
        missing += gone
        report.ok("B", label + " 대 SUMMARY_HEADERS", not gone,
                  "%d 항목 전부 존재" % len(items),
                  "SUMMARY_HEADERS 에 없음: " + ", ".join(gone))
    if total:
        report.add("B", "컬럼명 복제 합계", PASS if not missing else FAIL,
                   "대조 %d 항목 · 누락 %d" % (total, len(missing)))


def analyze_keyset(report, core, core_tree, reg_out, reg_note):
    """[B] 가 대조 기준으로 삼을 analyze / compose 의 키 집합.

    analyze 를 실제로 돌렸으면 그 결과의 키를 그대로 쓴다. analyze 에는 오류
    경로와 성공 경로의 return 이 따로 있어 AST 추정은 엉뚱한 쪽을 읽을 수 있다.
    testdata 가 없어 돌리지 못했을 때만 AST 로 폴백하고 그 사실을 표시한다.
    """
    try:
        cfg = core.build_config(None, [])
    except Exception:
        cfg = {}
    if reg_out is not None:
        out_sub = dict((k, set(v)) for k, v in reg_out.items() if isinstance(v, dict))
        comp = reg_out.get("composition") or {}
        return set(reg_out), out_sub, "실행 결과", set(comp), "실행 결과"

    fn = find_function(core_tree, "analyze")
    out_keys, out_sub = provided_keys(fn) if fn else (set(), {})
    out_sub = dict(out_sub)
    try:
        comp = core.compose([], cfg)
    except Exception:
        comp = {}
    out_sub["config"] = set(cfg)
    out_sub["composition"] = set(comp)
    if not out_keys:
        report.add("B", "core.analyze 반환 키 추출", FAIL,
                   "analyze 의 반환 dict 를 AST 로도 읽지 못했습니다")
    return (out_keys, out_sub, "AST 폴백 · " + (reg_note or "실행 결과 없음"),
            set(comp), "빈 입력 실행")


def check_js_keys(report, js, core, keyset, doc_keys, doc_sub):
    out_keys, out_sub, out_src, comp, comp_src = keyset
    scan = core.scan_inputs([], "")
    scan_sub = dict((k, set(v)) for k, v in scan.items() if isinstance(v, dict))

    targets = [
        ("OUT", "core.analyze", set(out_keys), out_sub, out_src),
        ("SCAN", "core.scan_inputs", set(scan), scan_sub, "실행 결과"),
        ("DOC", "js_docs (GLUE)", set(doc_keys), doc_sub, "AST"),
        ("C", "core.compose", set(comp), {}, comp_src),
    ]
    for obj, source, keys, sub, origin in targets:
        refs = js_member_refs(js, obj)
        if obj == "C" and not refs:
            report.add("B", "C.* 대 " + source, FAIL,
                       "composition 별칭 참조를 찾지 못했습니다")
            continue
        bad = []
        n = 0
        for key, subkey in refs:
            n += 1
            if key not in keys:
                bad.append("%s.%s" % (obj, key))
                continue
            if subkey and subkey not in JS_PROPS and key in sub and sub[key]:
                if subkey not in sub[key]:
                    bad.append("%s.%s.%s" % (obj, key, subkey))
        bad = sorted(set(bad))
        report.ok("B", obj + ".* 대 " + source, not bad,
                  "참조 %d 건 전부 존재 · 제공 키 %d (%s)" % (n, len(keys), origin),
                  "제공되지 않는 키: " + ", ".join(bad) + " (%s)" % origin)

    # fillBatchOptions 의 guess 인자 (SCAN.guess 하위키)
    body = js_function_body(js, "fillBatchOptions")
    if body is not None and isinstance(scan.get("guess"), dict):
        used = sorted(set(k for k, _ in js_member_refs(body, "guess")))
        gone = [k for k in used if k not in scan["guess"]]
        report.ok("B", "guess.* 대 scan_inputs['guess']", not gone,
                  "참조 %s 전부 존재" % ", ".join(used) if used else "참조 없음",
                  "없는 키: " + ", ".join(gone))


def check_readconfig_covers_design(report, js, core):
    keys = list(core.DESIGN_KEYS)
    parts, source = [], []
    for fname in ("buildDesignForm", "readConfig"):
        body = js_function_body(js, fname)
        if body is None:
            continue
        parts.append(body)
        source.append(fname)
    if not parts:
        parts, source = [js], ["script 전체 (함수 추출 실패)"]
    region = "\n".join(parts)
    literal = set(js_strings(region))
    dynamic = set()
    if re.search(r"DOC\.design", region):
        # DOC.design 을 순회해 폼과 cfg 를 만들면 DESIGN_DOC 키는 자동으로 덮인다
        dynamic = set(k for k, _lab, _memo in core.DESIGN_DOC)
    covered = literal | dynamic
    gone = [k for k in keys if k not in covered]
    report.ok("B", "readConfig 가 DESIGN_KEYS 를 덮는가", not gone,
              "%d/%d 키 · 근거 %s" % (len(keys) - len(gone), len(keys), "+".join(source)),
              "index.html 이 보내지 않는 키: " + ", ".join(gone))


# 이슈 2 에서 CONST / RULES 로 옮긴 옛 모듈 상수. 이름이 남아 있으면 이동이 덜 끝난 것.
MOVED_AWAY = ("_AL_MATCH", "_AL_MIS", "_AL_GAP", "_AL_FLANK",
              "_REPEAT_MAX_PERIOD", "_EXO_TRIM", "_OVERLONG_ZONE", "_FR4_MOTIF")


def sheets_for_check(core, reg_out):
    """시트 구조 검사용. 실행 결과가 있으면 그걸 쓰고, 없으면 빈 입력으로 만든다."""
    if reg_out is not None and reg_out.get("sheets"):
        return reg_out["sheets"]
    cfg = core.build_config(None, [])
    return core.build_sheets([], [], cfg, core.compose([], cfg), {})


def check_rules_doc(report, core):
    a = set(core.RULES)
    b = set(k for k, _d in core.RULES_DOC)
    report.ok("B", "RULES 대 RULES_DOC 키", a == b,
              "%d 키 일치" % len(a),
              "RULES 만: %s / RULES_DOC 만: %s"
              % (sorted(a - b) or "-", sorted(b - a) or "-"))


def check_moved_constants(report, core):
    left = [n for n in MOVED_AWAY if hasattr(core, n)]
    src = read_text("core.py")
    textual = [n for n in MOVED_AWAY if n in src]
    bad = sorted(set(left) | set(textual))
    report.ok("B", "옛 모듈 상수 잔존", not bad,
              "%d 개 전부 CONST/RULES 로 이동 완료" % len(MOVED_AWAY),
              "core.py 에 남아 있음: " + ", ".join(bad))


def check_rules_exposed(report, core, glue, sheets):
    """RULES 가 화면(js_docs)과 05_실행설정 시트에 실제로 노출되는지."""
    keys = [k for k, _d in core.RULES_DOC]

    doc_ok, doc_note = False, "GLUE 의 js_docs 를 읽지 못했습니다"
    if glue:
        fn = find_function(ast.parse(glue), "js_docs")
        if fn is not None:
            d = _returned_dict(fn.body[-1].value) if isinstance(fn.body[-1], ast.Return) else None
            got = set(_dict_literal_keys(d)) if d is not None else set()
            hit = sorted(k for k in got if "rule" in k.lower())
            doc_ok = bool(hit)
            doc_note = ("js_docs 키 " + ", ".join(hit)) if hit else \
                       "js_docs 반환에 RULES 관련 키가 없습니다"
    report.ok("B", "js_docs 가 RULES 를 내보내는가", doc_ok, doc_note, doc_note)

    rows = []
    for sh in sheets:
        if str(sh.get("title", "")).startswith("05"):
            rows = sh.get("rows") or []
            break
    listed = [r[1] for r in rows if len(r) > 1 and "알고리즘" in str(r[0])]
    gone = [k for k in keys if k not in listed]
    report.ok("B", "05_실행설정 에 RULES 행", rows and not gone,
              "'고정 상수 (알고리즘)' 행 %d 건 · RULES %d 키 전부 기록"
              % (len(listed), len(keys)),
              ("05_실행설정 시트를 찾지 못했습니다" if not rows
               else "시트에 없는 RULES 키: " + ", ".join(gone)))


def check_cfg_doc(report, core):
    a = set(k for k, _v in core.CFG_DEFAULTS)
    b = set(core.CFG_DOC)
    report.ok("B", "CFG_DEFAULTS 대 CFG_DOC 키", a == b,
              "%d 키 일치" % len(a),
              "DEFAULTS 만: %s / DOC 만: %s"
              % (sorted(a - b) or "-", sorted(b - a) or "-"))


def check_design_doc(report, core):
    a = [k for k, _v in core.DESIGN_DEFAULTS]
    b = set(k for k, _lab, _memo in core.DESIGN_DOC)
    gap = [k for k in a if k not in b]
    new = [k for k in gap if k not in KNOWN_DESIGN_DOC_GAP]
    extra = sorted(b - set(a))
    known = not new and not extra
    if new:
        detail = "DESIGN_DOC 에 없는 새 키 %d 건: %s" % (len(new), ", ".join(new))
    elif extra:
        detail = "DESIGN_DEFAULTS 에 없는데 DOC 에만 있음: " + ", ".join(extra)
    else:
        detail = "DEFAULTS %d · DOC %d · 알려진 누락 %d 건 (%s)" % (
            len(a), len(b), len(gap), ", ".join(gap) + " ← 이슈 3" if gap else "없음")
    report.add("B", "DESIGN_DEFAULTS 대 DESIGN_DOC 키", PASS if known else FAIL, detail)


# =============================================================================
#  [C] 누출
# =============================================================================
_SEQ_RE = re.compile(r"(?<![A-Za-z0-9_])[ACGT]{6,}(?![A-Za-z0-9_])")


def check_no_sequence(report, html):
    stripped = html.replace(MARKSEQ_REGEX, "")
    hits = []
    for m in _SEQ_RE.finditer(stripped):
        line = stripped.count("\n", 0, m.start()) + 1
        hits.append("%d행 %s" % (line, m.group(0)))
    report.ok("C", "index.html 6nt 이상 ACGT 리터럴", not hits,
              "0 건 (markSeq 정규식은 문자 클래스이므로 제외)",
              "발견: " + " / ".join(hits[:4]))


def check_no_const_leak(report, html, core):
    hits = []
    for k, v in core.CONST.items():
        if isinstance(v, str) and v and v in html:
            hits.append(k)
    report.ok("C", "CONST 서열이 index.html 에 등장", not hits,
              "문자열 상수 %d 개 전부 미등장"
              % sum(1 for v in core.CONST.values() if isinstance(v, str)),
              "index.html 에 등장: " + ", ".join(hits))


# =============================================================================
#  [D] 회귀
# =============================================================================
def testdata_inputs():
    d = os.path.join(ROOT, TESTDATA)
    if not os.path.isdir(d):
        return None, None
    ab1 = sorted(n for n in os.listdir(d) if n.lower().endswith(".ab1"))
    fa = [n for n in os.listdir(d) if n.lower().endswith((".fa", ".fasta", ".fas"))]
    if not ab1 or not fa:
        return None, None
    files = []
    for n in ab1:
        with open(os.path.join(d, n), "rb") as fh:
            files.append((n, fh.read()))
    pick = "scFv_primers.fa" if "scFv_primers.fa" in fa else sorted(fa)[0]
    with io.open(os.path.join(d, pick), encoding="utf-8", errors="ignore") as fh:
        return files, fh.read()


def norm_cell(v):
    """xlsx_writer 가 셀에 쓰는 값과 openpyxl 이 읽어온 값을 같은 자리로 옮긴다."""
    if v is None or v == "":
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, int):
        return v
    if isinstance(v, float):
        if v.is_integer() and abs(v) < 1e15:
            return int(v)
        return round(v, 10)
    return str(v)


D_NAMES = ["param_hash", "클론별 판정·길이", "클론별 혼합(%)", "CDR3-H3 중앙값",
           "조성 VL 순서", "조성 JH 순서", "fasta 통과 클론 수",
           "xlsx 바이트 재현성", "xlsx openpyxl 재검증"]


def run_analysis(core):
    """[B] 와 [D] 가 함께 쓰는 analyze 실행. (out, 사유, 상태) 를 반환한다.

    [B] 는 이 결과의 실제 키를 대조에 쓴다. analyze 에는 오류 경로와 성공 경로의
    return 이 따로 있어, AST 로 반환 키를 읽으면 엉뚱한 쪽을 볼 수 있기 때문이다.
    """
    files, ptext = testdata_inputs()
    if files is None:
        return None, "testdata/ 에 .ab1 또는 프라이머 FASTA 가 없습니다", SKIP
    meta = {"runtime": "verify.py " + VERIFY_VERSION,
            "timestamp": "0000-00-00 00:00:00", "primer_file": "scFv_primers.fa"}
    try:
        out = core.analyze(files, ptext, None, meta)
    except Exception as e:
        return None, "analyze 예외 %s: %s" % (type(e).__name__, e), FAIL
    if not out.get("ok"):
        return None, "analyze 중단: " + "; ".join(out.get("errors", [])), FAIL
    return out, "", PASS


def check_regression(report, out, note, status):
    if out is None:
        for n in D_NAMES:
            report.add("D", n, status, note)
        return

    import xlsx_writer

    report.ok("D", "param_hash", out["config"]["param_hash"] == EXPECT["param_hash"],
              EXPECT["param_hash"],
              "기대 %s 실제 %s" % (EXPECT["param_hash"], out["config"]["param_hash"]))

    qc = dict((q["id"], q) for q in out["qc"])
    bad = []
    for e in EXPECT["clones"]:
        q = qc.get(e["id"])
        if q is None:
            bad.append(e["id"] + " 없음")
            continue
        if q["verdict"] != e["verdict"]:
            bad.append("%s verdict %s(기대 %s)" % (e["id"], q["verdict"], e["verdict"]))
        if q["insert_bp"] != e["insert"]:
            bad.append("%s insert %s(기대 %d)" % (e["id"], q["insert_bp"], e["insert"]))
        for f in ("d1", "d2"):
            if f in e and q[f] != e[f]:
                bad.append("%s %s %s(기대 %d)" % (e["id"], f, q[f], e[f]))
    report.ok("D", "클론별 판정·길이", not bad,
              "%d 클론 · verdict / insert / d1 / d2 일치" % len(EXPECT["clones"]),
              " / ".join(bad))

    badm = []
    for e in EXPECT["clones"]:
        q = qc.get(e["id"]) or {}
        got = None if q.get("mix_pct") is None else round(q["mix_pct"], 1)
        if got != e["mix"]:
            badm.append("%s %s(기대 %s)" % (e["id"], got, e["mix"]))
    report.ok("D", "클론별 혼합(%)", not badm,
              " / ".join("%s %.1f" % (e["id"], e["mix"]) for e in EXPECT["clones"]),
              " / ".join(badm))

    comp = out["composition"]
    report.ok("D", "CDR3-H3 중앙값", comp["cdr3_median"] == EXPECT["cdr3_median"],
              str(EXPECT["cdr3_median"]),
              "기대 %s 실제 %s" % (EXPECT["cdr3_median"], comp["cdr3_median"]))
    for label, key in (("조성 VL 순서", "vl"), ("조성 JH 순서", "jh")):
        got = [tuple(x) for x in comp[key]]
        want = [tuple(x) for x in EXPECT[key]]
        report.ok("D", label, got == want, str(want), "기대 %s 실제 %s" % (want, got))

    n_fa = out["fasta"].count(">")
    report.ok("D", "fasta 통과 클론 수", n_fa == EXPECT["fasta_n"],
              str(EXPECT["fasta_n"]),
              "기대 %d 실제 %d" % (EXPECT["fasta_n"], n_fa))

    sheets = out["sheets"]
    b1 = xlsx_writer.write_xlsx(sheets)
    b2 = xlsx_writer.write_xlsx(sheets)
    report.ok("D", "xlsx 바이트 재현성", b1 == b2,
              "두 번 호출 결과 동일 · %d bytes" % len(b1),
              "두 번 호출 결과가 다릅니다 (%d vs %d bytes)" % (len(b1), len(b2)))

    check_xlsx_roundtrip(report, sheets, b1)


def check_xlsx_roundtrip(report, sheets, data):
    try:
        import openpyxl
    except ImportError:
        report.add("D", "xlsx openpyxl 재검증", SKIP, "openpyxl 없음")
        return
    wb = openpyxl.load_workbook(io.BytesIO(data), data_only=True)
    bad, ncell = [], 0
    if len(wb.worksheets) != len(sheets):
        report.add("D", "xlsx openpyxl 재검증", FAIL,
                   "시트 수 기대 %d 실제 %d" % (len(sheets), len(wb.worksheets)))
        return
    for si, sheet in enumerate(sheets):
        ws = wb.worksheets[si]
        headers = list(sheet.get("headers") or [])
        rows = list(sheet.get("rows") or [])
        note = sheet.get("note")
        ncol = max([len(headers)] + [len(r) for r in rows]) if (headers or rows) else 1
        top = 3 if note else 1
        want = []
        if note:
            want.append((1, [note] + [None] * (ncol - 1)))
        if headers:
            want.append((top, list(headers) + [None] * (ncol - len(headers))))
        for i, row in enumerate(rows):
            want.append((top + 1 + i, list(row) + [None] * (ncol - len(row))))
        for rn, vals in want:
            for j, v in enumerate(vals, start=1):
                ncell += 1
                got = norm_cell(ws.cell(row=rn, column=j).value)
                exp = norm_cell(v)
                if got != exp:
                    bad.append("%s!%s%d 기대 %r 실제 %r"
                               % (sheet["title"], openpyxl.utils.get_column_letter(j),
                                  rn, exp, got))
    report.ok("D", "xlsx openpyxl 재검증", not bad,
              "%d 시트 · %d 셀 전부 일치" % (len(sheets), ncell),
              " / ".join(bad[:3]) + (" 외 %d" % (len(bad) - 3) if len(bad) > 3 else ""))


# =============================================================================
#  [E] 단위 — 실측 testdata 가 밟지 않는 랜드마크 분기
# =============================================================================
# testdata 4 클론의 랜드마크 gap 은 {0, 1} 뿐이라 GAP# 분기가 한 번도 실행되지
# 않습니다. 합성 서열로 check_landmark 를 직접 불러 그 구간을 덮습니다.
UNIT_PRE = "T" * 20                 # 랜드마크와 무관한 5' 패딩
UNIT_SUF = "A" * 20                 # 랜드마크와 무관한 3' 패딩
UNIT_ALIEN = ("ATTC" * 12)[:45]     # 랜드마크 자리에 넣을 무관한 45 nt
UNIT_SUB2 = (6, 22)                 # 치환 위치 — 허용치(lm_max_sub) 이내
UNIT_SUB3 = (6, 22, 38)             # 치환 위치 — 허용치 초과
UNIT_DEL_AT = 12                    # 결실 시작 위치 (동종중합 구간 밖)
# 5' G-run 경계. 여기서 2 nt 를 지우면 갭 1 개와 치환 2 개가 박빙이라
# AL_GAP 이 -3 이 되면 S2G1/WARN 으로 내려앉는다. AL_GAP 회귀 감지용.
# 같은 위치라도 n=1 / n=3 은 격차가 커서 갈리지 않으므로 n=2 만 쓴다.
UNIT_DEL_TIGHT = 3
_TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}

#        태그   설명                      변형종류  인자                기대 status  기대 level
UNIT_CASES = [
    ("E1", "변형 없음",             "sub",   (),                    "OK",     "OK"),
    ("E2", "치환 2 개 (허용 이내)", "sub",   UNIT_SUB2,             "OK",     "OK"),
    ("E3", "치환 3 개 (허용 초과)", "sub",   UNIT_SUB3,             "S3G0",   "WARN"),
    ("E4", "1 nt 결실",             "del",   (UNIT_DEL_AT, 1),      "S0G1",   "WARN"),
    ("E5", "2 nt 결실",             "del",   (UNIT_DEL_AT, 2),      "GAP2",   "FAIL"),
    ("E6", "3 nt 결실",             "del",   (UNIT_DEL_AT, 3),      "GAP3",   "FAIL"),
    ("E7", "무관한 서열",           "alien", None,                  "ABSENT", "FAIL"),
    ("E8", "covered=False",         "na",    None,                  "NA",     "NA"),
    ("E9", "결실 위치 3 · 2 nt 결실", "del", (UNIT_DEL_TIGHT, 2),   "GAP2",   "FAIL"),
]


def mutate_landmark(lm, kind, arg):
    if kind == "sub":
        b = list(lm)
        for p in arg:
            b[p] = _TRANSITION[b[p]]
        return "".join(b)
    if kind == "del":
        pos, n = arg
        return lm[:pos] + lm[pos + n:]
    if kind == "alien":
        return UNIT_ALIEN
    return lm


def check_landmark_units(report, core):
    cfg = core.build_config(None, [])
    lm = core.CONST["QC2"]

    # 패딩과 대체서열이 랜드마크와 우연히 맞지 않는지 먼저 확인한다.
    su, _pu = core.ungapped_scan(lm, UNIT_PRE + UNIT_ALIEN + UNIT_SUF)
    report.ok("E", "패딩·대체서열 무관성", su > cfg["lm_max_sub"],
              "랜드마크 없는 서열의 최소 불일치 %d > 허용 치환 %d" % (su, cfg["lm_max_sub"]),
              "패딩이 랜드마크와 우연히 맞습니다 (불일치 %d ≤ 허용 %d)"
              % (su, cfg["lm_max_sub"]))

    for tag, label, kind, arg, want_st, want_lv in UNIT_CASES:
        seq = UNIT_PRE + mutate_landmark(lm, kind, arg) + UNIT_SUF
        r = core.check_landmark(lm, seq, [], kind != "na", cfg)
        good = (r["status"] == want_st and r["level"] == want_lv)
        detail = "status %s · level %s" % (want_st, want_lv)
        if tag == "E1":
            good = good and r["pos"] == len(UNIT_PRE)
            detail += " · 위치 %d" % r["pos"]
        report.ok("E", tag + " " + label, good, detail,
                  "기대 %s/%s · 실제 %s/%s (sub %d, gap %d, pos %d)"
                  % (want_st, want_lv, r["status"], r["level"],
                     r["sub"], r["gap"], r["pos"]))


# --- index.html 의 badge() 가 위 5 종을 모두 처리하는가 ------------------------
BADGE_EXPECT = [("OK", "b-pass"), ("S1G0", "b-warn"), ("GAP2", "b-fail"),
                ("ABSENT", "b-fail"), ("NA", "b-none")]

_BADGE_RULE_RE = re.compile(
    r"^\s*if\s*\((?P<cond>.+?)\)\s*return\s*'<span class=\"badge (?P<cls>b-[a-z]+)\"")
_BADGE_DEFAULT_RE = re.compile(
    r"^\s*return\s*'<span class=\"badge (?P<cls>b-[a-z]+)\"")


def _js_unquote(text):
    return text.replace('\\"', '"').replace("\\\\", "\\")


def eval_badge_cond(cond, s):
    """badge() 의 조건식을 제한 해석한다. 해석하지 못하면 None."""
    for part in cond.split("||"):
        part = part.strip()
        m = re.fullmatch(r's\s*===\s*"((?:[^"\\]|\\.)*)"', part)
        if m:
            if s == _js_unquote(m.group(1)):
                return True
            continue
        m = re.fullmatch(r"/(.+)/[a-z]*\.test\(s\)", part)
        if m:
            try:
                if re.search(m.group(1), s):
                    return True
            except re.error:
                return None
            continue
        return None
    return False


def badge_class(body, s):
    """badge() 본문을 위에서부터 따라가며 s 가 받을 클래스와 명시/기본 여부."""
    for line in body.splitlines():
        m = _BADGE_RULE_RE.match(line)
        if m:
            got = eval_badge_cond(m.group("cond"), s)
            if got is None:
                return None, "해석실패"
            if got:
                return m.group("cls"), "명시"
            continue
        m = _BADGE_DEFAULT_RE.match(line)
        if m:
            return m.group("cls"), "기본값"
    return None, "return 없음"


def check_badge_states(report, js):
    body = js_function_body(js, "badge")
    if body is None:
        report.add("E", "badge() 랜드마크 상태 5 종", FAIL,
                   "index.html 에서 badge() 를 찾지 못했습니다")
        return
    bad, shown = [], []
    for s, want in BADGE_EXPECT:
        cls, how = badge_class(body, s)
        shown.append("%s:%s%s" % (s, (cls or how).replace("b-", ""),
                                  "*" if how == "기본값" else ""))
        if cls != want:
            bad.append("%s 기대 %s 실제 %s(%s)" % (s, want, cls, how))
    report.ok("E", "badge() 랜드마크 상태 5 종", not bad,
              " ".join(shown) + "  (*기본값)", " / ".join(bad))


# =============================================================================
#  실행
# =============================================================================
def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    report = Report()
    html = read_text(HTML_FILE)
    blocks = script_blocks(html)
    js_all = "\n".join(blocks)
    glue = glue_source(js_all)
    js = js_all.replace(glue, "") if glue else js_all      # 파이썬 접착부 제외한 JS

    guard(report, "A", "AST 파싱", check_ast, report)
    guard(report, "A", "index.html JS 구문", check_js_syntax, report, html)
    guard(report, "A", "포맷 자리표시자 대조", check_format_arity, report)
    guard(report, "A", "함수 안 미정의 이름", check_undefined_names, report)

    sys.path.insert(0, ROOT)
    try:
        import core
    except Exception as e:
        core = None
        for sec in ("B", "C", "D"):
            report.add(sec, "core.py import", FAIL,
                       "core.py 를 불러올 수 없어 검사를 진행하지 못했습니다: %s: %s"
                       % (type(e).__name__, e))

    reg_out, reg_note, reg_status = None, "core.py 불러오기 실패", FAIL
    if core is not None:
        reg_out, reg_note, reg_status = run_analysis(core)

    if core is not None:
        core_tree = parse_py("core.py")
        if glue:
            gfn = find_function(ast.parse(glue), "js_docs")
            doc_keys, doc_sub = provided_keys(gfn) if gfn else (set(), {})
        else:
            doc_keys, doc_sub = set(), {}
        if not doc_keys:
            report.add("B", "js_docs 반환 키 추출", FAIL, "GLUE 의 js_docs 를 읽지 못했습니다")

        keyset = analyze_keyset(report, core, core_tree, reg_out, reg_note)

        guard(report, "B", "SUMMARY_HEADERS 대조", check_summary_headers, report, js, core)
        guard(report, "B", "JS 참조 키 대조", check_js_keys,
              report, js, core, keyset, doc_keys, doc_sub)
        guard(report, "B", "readConfig 가 DESIGN_KEYS 를 덮는가",
              check_readconfig_covers_design, report, js, core)
        guard(report, "B", "CFG_DEFAULTS 대 CFG_DOC 키", check_cfg_doc, report, core)
        guard(report, "B", "DESIGN_DEFAULTS 대 DESIGN_DOC 키", check_design_doc, report, core)
        guard(report, "B", "RULES 대 RULES_DOC 키", check_rules_doc, report, core)
        guard(report, "B", "옛 모듈 상수 잔존", check_moved_constants, report, core)
        guard(report, "B", "RULES 노출", check_rules_exposed,
              report, core, glue, sheets_for_check(core, reg_out))

        guard(report, "C", "index.html 6nt 이상 ACGT 리터럴",
              check_no_sequence, report, html)
        guard(report, "C", "CONST 서열이 index.html 에 등장",
              check_no_const_leak, report, html, core)

        guard(report, "D", "회귀", check_regression,
              report, reg_out, reg_note, reg_status)

        guard(report, "E", "check_landmark 단위", check_landmark_units, report, core)

    guard(report, "E", "badge() 랜드마크 상태 5 종", check_badge_states, report, js)

    ver = "core %s / notebook %s" % (core.CORE_VERSION, core.NB_VERSION) \
        if core is not None else "core.py 불러오기 실패"
    print("pAIM1 scFv QC — verify.py %s   (%s)" % (VERIFY_VERSION, ver))
    print("")
    print_table(report)
    print("")
    parts = []
    for s in ("A", "B", "C", "D", "E"):
        parts.append("[%s] 통과 %d 실패 %d 건너뜀 %d"
                     % (s, report.count(s, PASS), report.count(s, FAIL),
                        report.count(s, SKIP)))
    print("   ".join(parts))
    return 1 if report.failed() else 0


if __name__ == "__main__":
    sys.exit(main())
