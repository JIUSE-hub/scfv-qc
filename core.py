# -*- coding: utf-8 -*-
"""
core.py — pAIM1 scFv library Sanger read QC : 과학 로직 코어
==============================================================================
scFv_sequencing_QC.ipynb (노트북 v1.0) 의 Cell 1/3/4/5 판정 로직을
UI 와 파일 입출력에서 분리해 옮긴 모듈입니다.

설계 원칙
  1. 표준 라이브러리만 사용합니다.
     - .ab1 파싱      : Biopython SeqIO(abi) 를 struct 기반 파서로 대체
     - 아미노산 번역   : Biopython Seq.translate 를 코돈 표로 대체
                        (IUPAC 3375 코돈 전수 대조로 동치 확인)
     - 정렬           : 노트북과 동일한 자체 DP
  2. 전역 가변 상태를 두지 않습니다. 판정 수치는 항상 cfg 인자로 전달합니다.
     (노트북의 전역 CFG 를 인자로 바꾼 것 외에 판정 로직은 동일합니다.)
  3. 화면 출력을 하지 않습니다. 모든 결과를 JSON 직렬화 가능한 값으로 반환합니다.
  4. 코드에 고정된 값은 두 계층입니다.
       CONST  서열과 프레임 규칙
       RULES  알고리즘 규칙 (정렬 점수, 탐색 범위, 플래그 심각도)
     둘 다 UI 로 조절하지 않으며, 05_실행설정 시트와 화면 하단에 전부 노출됩니다.
     판정 임계값은 CFG_DEFAULTS 를 통해 UI 로 제어하며 param_hash 에 포함됩니다.

사용법
    import core
    scan = core.scan_inputs(files, primer_text)          # 폼 채우기용
    out  = core.analyze(files, primer_text, overrides)   # 분석
        files       : [(파일명, bytes), ...]
        primer_text : 프라이머 FASTA 전문 (str)
        overrides   : 노트북 Cell 2 UI 에 해당하는 설정 dict
"""

import hashlib
import json
import re
import struct

CORE_VERSION = "1.2"
NB_VERSION = "1.0"          # 기준 노트북 버전

# =============================================================================
#  1. 코드에 고정된 값 : CONST(서열·프레임) 와 RULES(알고리즘)
#     둘 다 UI 로 조절하지 않고 param_hash 에도 들어가지 않습니다.
#     값은 노트북 Cell 1 과 동일하며, 이번 계층 분리에서 이동만 했습니다.
# =============================================================================
CONST = {
    "NotI": "GCGGCCGC",
    "AscI": "GGCGCGCC",
    "PELB_ATG": "ATGAAATACCTATTGCCTACG",
    "QC1": "GCTCAACCCGCAATGGCGGCCGCA",
    "QC2": "GGTGGCGGAGGGTCTGGAGGTGGGGGCTCAGGGGGCGGTGGATCC",
    "QC3": "GGCGCGCCTGCGCACCATCACCATCACCATGGC",
    "QC4": "GATTACAAAGACCATGACGGAGACTACAAGGACCATGATATCGACTATAAGGATGATGATGACAAATAG",
    "STUFFER": "GGCTTTAATATCAAAGACACGTACATTCATTGG",
    "STUFFER_INSERT_BP": 386,
    "TAG_F1_For": "CGCAATGGCGGCCGC",
    "TAG_F2_Rev": "GGATCCACCGCCCCC",
    "TAG_F3_For": "TCAGGGGGCGGTGGATCC",
    "TAG_F3_Rev": "GTGCGCAGGCGCGCC",
    "MOTIF_F1_Rev": "ACAGTAATA",
    "MOTIF_F2_For": "TATTACTGT",
    "FR4_MOTIF": "WG.G",
    "FRAME_MOD": 2,
}

CONST_DOC = [
    ("NotI", "5' 클로닝 자리 인식서열"),
    ("AscI", "3' 클로닝 자리 인식서열"),
    ("PELB_ATG", "pelB 개시코돈 앵커 (ORF 시작점 탐색)"),
    ("QC1", "pelB 3' + NotI + frame-keeping A"),
    ("QC2", "(G4S)3 링커 전체"),
    ("QC3", "AscI + His6"),
    ("QC4", "3xFLAG + amber TAG"),
    ("STUFFER", "모클론 스터퍼 VH 검출용"),
    ("STUFFER_INSERT_BP", "스터퍼 보유 시 인서트 길이 (bp)"),
    ("TAG_F1_For", "프라이머 자동분류용 5' 태그"),
    ("TAG_F2_Rev", "프라이머 자동분류용 5' 태그"),
    ("TAG_F3_For", "프라이머 자동분류용 5' 태그"),
    ("TAG_F3_Rev", "프라이머 자동분류용 5' 태그"),
    ("MOTIF_F1_Rev", "프라이머 자동분류용 내부 모티프"),
    ("MOTIF_F2_For", "프라이머 자동분류용 내부 모티프"),
    ("FR4_MOTIF", "VH FR4 보존 모티프 (Trp-Gly-Xxx-Gly). CDR3-H3 끝 지점 결정"),
    ("FRAME_MOD", "프레임 유지 조건 : insert 길이 % 3 == 이 값"),
]

# --- 알고리즘 규칙 -----------------------------------------------------------
# 판정 임계값(CFG_DEFAULTS)과 달리 UI 로 조절하지 않지만, 판정 결과에는
# 직접 작용합니다. CONST 와 나란히 05_실행설정 시트와 화면 하단에 노출합니다.
RULES = {
    "AL_MATCH": 1,
    "AL_MIS": -1,
    "AL_GAP": -2,
    "AL_FLANK": 12,
    "REPEAT_MAX_PERIOD": 60,
    "EXO_TRIM": 3,
    "OVERLONG_ZONE": 5,
    "FLAG_SEV": {
        "CONCATEMER": 3, "PARENTAL": 3, "NO_NOTI": 3, "NO_ASCI": 3, "NO_LINKER": 3,
        "TOO_SHORT": 3, "TOO_LONG": 3, "FRAMESHIFT": 3, "INTERNAL_STOP": 3,
        "LINKER_DEL": 3, "QC_DEL": 3, "QC_ABSENT": 3,
        "ABERRANT_D1": 3, "ABERRANT_D2": 3, "TANDEM_REPEAT": 3, "MIXED": 3,
        "LONG_INSERT?": 2, "LOW_COVERAGE": 2,
        "QC_WARN": 1,
    },
}

RULES_DOC = [
    ("AL_MATCH", "랜드마크 정렬의 일치 점수. 치환/갭 개수 산출에 영향"),
    ("AL_MIS", "랜드마크 정렬의 불일치 점수"),
    ("AL_GAP", "랜드마크 정렬의 갭 점수. 값이 클수록 갭보다 치환으로 정렬됨"),
    ("AL_FLANK", "정렬 DP 창의 양쪽 여유 (nt). 검출 가능한 결실 폭의 상한"),
    ("REPEAT_MAX_PERIOD", "탠덤 반복 탐색의 최대 주기 (nt). 이보다 긴 반복은 검출하지 않음"),
    ("EXO_TRIM", "CDR3 경계 정합성 검사에서 제외하는 양끝 nt. "
                 "proofreading exonuclease 모자이크 보정"),
    ("OVERLONG_ZONE", "F1_For 불일치 편중을 점검하는 판별구간 앞쪽 범위 (nt)"),
    ("FLAG_SEV", "플래그별 심각도. 여러 플래그가 동시에 붙을 때 어느 것을 verdict 로 삼을지 결정"),
]


def sev_flags(level):
    """FLAG_SEV 에서 해당 심각도의 플래그 이름을 정의 순서대로 뽑는다."""
    return [k for k, v in RULES["FLAG_SEV"].items() if v == level]


def rule_value_text(key):
    """RULES 값을 표·화면에 넣을 문자열. dict 는 길어서 요약한다."""
    v = RULES[key]
    if isinstance(v, dict):
        levels = sorted(set(v.values()), reverse=True)
        return "%d 종 매핑 (심각도 %s)" % (len(v), "/".join(str(x) for x in levels))
    return str(v)


# =============================================================================
#  2. 판정 파라미터 기본값  (노트북 Cell 2 와 동일)
# =============================================================================
CFG_DEFAULTS = [
    ("trim_q", 20),
    ("trim_win", 20),
    ("insert_min", 685),
    ("insert_max", 840),
    ("d1_min", 340),
    ("d1_max", 440),
    ("d2_min", 295),
    ("d2_max", 345),
    ("lm_max_sub", 2),
    ("lm_gap_warn", 1),
    ("lm_gap_fail", 2),
    ("mix_ratio", 0.35),
    ("mix_pct", 8.0),
    ("repeat_min_len", 15),
    ("primer_max_mismatch", 2),
    ("primer_margin_min", 0),
]
THRESH_KEYS = [k for k, _v in CFG_DEFAULTS]
CFG_DEFAULT_MAP = dict(CFG_DEFAULTS)

CFG_DOC = {
    "trim_q": ("Q 트리밍 임계", "이 Q 미만 구간을 read 양끝에서 잘라냅니다"),
    "trim_win": ("Q 트리밍 창 크기", "이 길이의 이동창 평균 Q 로 판단합니다 (nt)"),
    "insert_min": ("인서트 길이 하한", "NotI 첫염기 ~ AscI 끝염기 (bp)"),
    "insert_max": ("인서트 길이 상한", "이 범위를 벗어나면 ABERRANT (bp)"),
    "d1_min": ("d1 하한", "NotI ~ 링커 시작 = VH 전체 (bp)"),
    "d1_max": ("d1 상한", "CDR3-H3 길이에 따라 변동 (bp)"),
    "d2_min": ("d2 하한", "링커 끝 ~ AscI = VL 전체 (bp)"),
    "d2_max": ("d2 상한", "(bp)"),
    "lm_max_sub": ("랜드마크 허용 치환", "이 개수까지는 시퀀싱 오독으로 간주"),
    "lm_gap_warn": ("랜드마크 갭 WARN", "이 개수 이상이면 WARN"),
    "lm_gap_fail": ("랜드마크 갭 FAIL", "이 개수 이상이면 FAIL (올리고 결실)"),
    "mix_ratio": ("혼합 판정 피크비", "2순위/1순위 피크 세기 비율"),
    "mix_pct": ("혼합 판정 위치비율", "위 비율을 넘는 위치가 이 % 이상이면 MIXED"),
    "repeat_min_len": ("탠덤 반복 최소 길이", "V 구간에서 이 길이 이상 반복 시 표시 (nt)"),
    "primer_max_mismatch": ("프라이머 허용 미스매치", "판별 구간에서 이 개수까지 허용 (개)"),
    "primer_margin_min": ("프라이머 모호성 마진", "1위-2위 점수차가 이 값 이하면 모호성 집합으로 보고"),
}

NOSEL = "(지정 안 함)"
POOL_OPTS = ["분리", "pool"]
CDNA_OPTS = ["oligo(dT)", "random hexamer", "혼합/불명"]

DESIGN_DEFAULTS = [
    ("f1_for_mode", "분리"),
    ("f1_rev_mode", "분리"),
    ("f2_for_mode", "분리"),
    ("f2_rev_mode", "pool"),
    ("f3_for_mode", "분리"),
    ("f3_rev_mode", "pool"),
    ("f3_product_pooled", True),
    ("batch_vh_family", NOSEL),
    ("batch_chain", NOSEL),
    ("cdna_frag1", "oligo(dT)"),
    ("cdna_frag2", "oligo(dT)"),
    ("cdna_frag3", "oligo(dT)"),
    ("rna_bone_marrow", True),
    ("rna_peripheral", True),
]
DESIGN_KEYS = [k for k, _v in DESIGN_DEFAULTS]
DESIGN_DEFAULT_MAP = dict(DESIGN_DEFAULTS)

# 05_실행설정에서는 아래 두 키를 개별 행으로 내지 않고 cfg["rna_source"] 로 합쳐
# 한 행에 씁니다. 합친 문자열이 두 값의 네 조합을 모두 구분하므로 정보 손실이 없습니다.
RNA_KEYS = ("rna_bone_marrow", "rna_peripheral")

DESIGN_DOC = [
    ("f1_for_mode", "Fragment 1 For", "VH FR1 프라이머를 family 별로 나눠 PCR 했는지"),
    ("f1_rev_mode", "Fragment 1 Rev", "CDR3 경계 프라이머"),
    ("f2_for_mode", "Fragment 2 For", "CDR3 경계 프라이머"),
    ("f2_rev_mode", "Fragment 2 Rev", "JH + 링커 프라이머"),
    ("f3_for_mode", "Fragment 3 For", "VL FR1 프라이머"),
    ("f3_rev_mode", "Fragment 3 Rev", "VL J 프라이머"),
    ("f3_product_pooled", "frag3 산물 pooling", "over-PCR 전에 섞었는지"),
    ("batch_vh_family", "배치 지정 VH family", "불일치 시 WRONG_FAMILY"),
    ("batch_chain", "배치 지정 경쇄", "불일치 시 WRONG_CHAIN"),
    ("cdna_frag1", "Fragment 1 cDNA", "판정에 미사용. 기록용"),
    ("cdna_frag2", "Fragment 2 cDNA", "판정에 미사용. 기록용"),
    ("cdna_frag3", "Fragment 3 cDNA", "판정에 미사용. 기록용"),
    ("rna_bone_marrow", "RNA 출처 · Bone marrow", "판정에 미사용. 기록용"),
    ("rna_peripheral", "RNA 출처 · Peripheral leukocytes", "판정에 미사용. 기록용"),
]


# =============================================================================
#  3. 저수준 유틸
# =============================================================================
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

_RC_TABLE = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def rc(s):
    """역상보."""
    return s.translate(_RC_TABLE)[::-1]


def hit(pat_base, obs_base):
    """축퇴염기를 고려한 1염기 일치 판정."""
    return obs_base in IUPAC.get(pat_base, pat_base)


_BASES = "TCAG"
_AA_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
for _i, _a in enumerate(_BASES):
    for _j, _b in enumerate(_BASES):
        for _k, _c in enumerate(_BASES):
            CODON[_a + _b + _c] = _AA_TABLE[_i * 16 + _j * 4 + _k]

# 축퇴 코돈이 한 가지 아미노산으로만 풀리면 그 값, 두 가지면 IUPAC 모호 코드
_AMBIG_AA = {frozenset("DN"): "B", frozenset("EQ"): "Z", frozenset("IL"): "J"}
_AA_CACHE = {}


def _aa_of(codon):
    a = _AA_CACHE.get(codon)
    if a is not None:
        return a
    a = CODON.get(codon)
    if a is None:
        opts = [IUPAC.get(ch) for ch in codon]
        if len(codon) == 3 and all(opts):
            s = set(CODON[x + y + z] for x in opts[0] for y in opts[1] for z in opts[2])
            a = s.pop() if len(s) == 1 else _AMBIG_AA.get(frozenset(s), "X")
        else:
            a = "X"
    _AA_CACHE[codon] = a
    return a


def translate(nt):
    """표준 유전암호 번역. Biopython Seq.translate 와 동치."""
    return "".join(_aa_of(nt[i:i + 3]) for i in range(0, len(nt) - len(nt) % 3, 3))


# =============================================================================
#  4. .ab1 파싱  (표준 라이브러리만 사용)
# =============================================================================
_TRACE_KEYS = ("DATA9", "DATA10", "DATA11", "DATA12")
_ENTRY = ">4slhhllll"          # name(4) number etype esize nelem dsize doffset dhandle
_ENTRY_SIZE = 28


def parse_ab1(data):
    """ABIF 바이너리에서 서열 / 품질값 / 트레이스 / 피크위치를 꺼낸다."""
    if len(data) < 34 or data[:4] != b"ABIF":
        raise ValueError("ABIF 형식이 아닙니다")
    _n, _num, _et, _es, nelem, _ds, doff, _dh = struct.unpack(_ENTRY, data[6:6 + _ENTRY_SIZE])
    tags = {}
    for i in range(nelem):
        o = doff + i * _ENTRY_SIZE
        if o + _ENTRY_SIZE > len(data):
            break
        name, num, etype, esize, ne, dsize, dofs, _h = struct.unpack(
            _ENTRY, data[o:o + _ENTRY_SIZE])
        raw = data[dofs:dofs + dsize] if dsize > 4 else data[o + 20:o + 20 + dsize]
        tags[(name.decode("ascii", "replace"), num)] = (etype, ne, raw)

    def get(tag, num):
        v = tags.get((tag, num))
        if v is None:
            return None
        etype, ne, raw = v
        if etype == 2:
            return raw
        if etype == 1:
            return list(raw)
        if etype == 3:
            return list(struct.unpack(">%dH" % ne, raw[:2 * ne]))
        if etype == 4:
            return list(struct.unpack(">%dh" % ne, raw[:2 * ne]))
        if etype == 5:
            return list(struct.unpack(">%dl" % ne, raw[:4 * ne]))
        return raw

    pbas = get("PBAS", 2)
    if pbas is None:
        raise ValueError("염기 호출(PBAS) 태그가 없습니다")
    seq = pbas.decode("ascii", "replace").upper()
    pcon = get("PCON", 2)
    qual = list(pcon) if pcon is not None else []
    ploc = get("PLOC", 2) or get("PLOC", 1) or []
    trace = {}
    for k in _TRACE_KEYS:
        v = get("DATA", int(k[4:]))
        if v is not None:
            trace[k] = v
    return {"seq": seq, "qual": qual, "trace": trace, "ploc": list(ploc)}


# =============================================================================
#  5. 파일명 파싱  (노트북 Cell 1 과 동일)
# =============================================================================
_FNAME_RE = re.compile(
    r"^(?P<date>\d{6})_(?P<batch>.+?)_(?P<clone>[cC]\d+)_(?P<primer>.+)$")


def _basename(path):
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def parse_read_name(filename):
    stem = _basename(filename)
    if "." in stem:
        stem = stem.rsplit(".", 1)[0]
    m = _FNAME_RE.match(stem)
    if m:
        info = m.groupdict()
        info["parsed"] = True
    else:
        info = {"date": "", "batch": "", "clone": stem, "primer": "", "parsed": False}
    probe = (info["primer"] or stem).lower()
    if "rev" in probe:
        info["direction"], info["dir_guessed"] = "R", False
    elif "for" in probe:
        info["direction"], info["dir_guessed"] = "F", False
    else:
        info["direction"], info["dir_guessed"] = "F", True
    info["stem"] = stem
    info["id"] = info["clone"] if info["parsed"] else stem
    return info


# =============================================================================
#  6. 프라이머 FASTA 파싱  (노트북 Cell 1 과 동일)
# =============================================================================
ANALYSIS_GROUPS = ("F1_For", "F1_Rev", "F2_For", "F2_Rev", "F3_For", "F3_Rev")


def classify_group(seq):
    """헤더에 group= 이 없을 때만 쓰는 서열 기반 분류."""
    s = seq.upper()
    t1 = CONST["TAG_F1_For"]
    t3r = CONST["TAG_F3_Rev"]
    if s.startswith(t1):
        return "OVER_For" if len(s) - len(t1) <= 6 else "F1_For"
    if s.startswith(t3r):
        return "OVER_Rev" if len(s) - len(t3r) <= 6 else "F3_Rev"
    if s.startswith(CONST["TAG_F2_Rev"]):
        return "F2_Rev"
    if s.startswith(CONST["TAG_F3_For"]):
        return "F3_For"
    if CONST["MOTIF_F1_Rev"] in s:
        return "F1_Rev"
    if CONST["MOTIF_F2_For"] in s:
        return "F2_For"
    return "UNKNOWN"


def _read_fasta_entries(text):
    entries = []
    name = None
    meta = {}
    buf = []
    last_key = None
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if name is not None:
                entries.append((name, meta, "".join(buf)))
            fields = [f.strip() for f in line[1:].split("|")]
            head = fields[0].split()
            name = head[0] if head else "unnamed"
            meta = {}
            last_key = None
            for f in fields[1:]:
                if "=" in f:
                    k, v = f.split("=", 1)
                    last_key = k.strip().lower()
                    meta[last_key] = v.strip()
                elif last_key is not None and f:
                    meta[last_key] = meta[last_key] + "," + f
            buf = []
        else:
            buf.append(re.sub(r"[^A-Za-z]", "", line).upper())
    if name is not None:
        entries.append((name, meta, "".join(buf)))
    return entries


def _common_prefix_len(seqs):
    if len(seqs) < 2:
        return 0
    shortest = min(len(s) for s in seqs)
    n = 0
    while n < shortest and len(set(s[n] for s in seqs)) == 1:
        n += 1
    return n


def parse_primer_fasta(text):
    """(primers, trims, warns) 반환. trims 는 {'group|chain': 공통태그길이}."""
    warns = []
    primers = []
    for name, meta, seq in _read_fasta_entries(text):
        if not seq:
            warns.append(name + " : 서열이 비어 있어 제외했습니다.")
            continue
        group = meta.get("group", "")
        if not group:
            group = classify_group(seq)
            warns.append(name + " : group= 없음 -> 서열 태그로 " + group + " 자동분류")
        chain = meta.get("chain", "")
        if not chain:
            low = name.lower()
            if "-k-" in low:
                chain = "kappa"
            elif "-l-" in low:
                chain = "lambda"
            elif group.startswith("F1") or group.startswith("F2"):
                chain = "heavy"
            else:
                chain = "unknown"
            warns.append(name + " : chain= 없음 -> " + chain + " 추정")
        family = meta.get("family", "")
        if not family:
            family = "unknown"
            warns.append(name + " : family= 없음 -> unknown 처리")
        fams = [x.strip() for x in re.split(r"[,|]", family) if x.strip()]
        fams = [x for x in fams if x.lower() not in ("none", "-", "unknown")]
        primers.append({
            "name": name, "seq": seq, "len": len(seq),
            "group": group, "chain": chain, "family": family, "families": fams,
            "target": meta.get("target", ""), "fragment": meta.get("fragment", ""),
            "dir": meta.get("dir", ""), "tm": meta.get("tm", ""),
        })

    buckets = {}
    for p in primers:
        if p["group"] in ANALYSIS_GROUPS:
            buckets.setdefault((p["group"], p["chain"]), []).append(p)
    trims = {}
    for key, plist in buckets.items():
        t = _common_prefix_len([p["seq"] for p in plist])
        trims[key[0] + "|" + key[1]] = t
        for p in plist:
            p["core_trim"] = t
            p["core"] = p["seq"][t:]
    for p in primers:
        if "core" not in p:
            p["core_trim"] = 0
            p["core"] = p["seq"]
    return primers, trims, warns


# =============================================================================
#  7. 설정 만들기  (노트북 Cell 2 의 get_config 와 동일한 규격)
# =============================================================================
def vh_family_options(primers):
    return sorted(set(f for p in primers if p["group"] == "F1_For" for f in p["families"]))


def chain_options(primers):
    return sorted(set(p["chain"] for p in primers if p["group"] == "F3_For"))


def guess_batch(batches, primers):
    """파일명 배치 문자열에서 VH family / 경쇄를 추정한다."""
    fams = vh_family_options(primers)
    chains = chain_options(primers)
    hv = set()
    for b in batches:
        up = b.upper()
        for f in fams:
            if f.upper() in up:
                hv.add(f)
    hc = set()
    for b in batches:
        low = b.lower()
        if ("vk" in low) or ("kappa" in low):
            hc.add("kappa")
        if ("vl" in low) or ("lambda" in low):
            hc.add("lambda")
    hc = hc & set(chains)
    return {"vh": list(hv)[0] if len(hv) == 1 else NOSEL,
            "chain": list(hc)[0] if len(hc) == 1 else NOSEL}


def _coerce(value, default):
    """기본값의 타입에 맞춰 변환한다. param_hash 재현성에 필요."""
    if isinstance(default, bool):
        if isinstance(value, str):
            return value.strip().lower() in ("true", "1", "y", "yes", "on")
        return bool(value)
    if isinstance(default, int):
        return int(round(float(value)))
    if isinstance(default, float):
        return float(value)
    return str(value)


def build_config(overrides=None, batches=None):
    """판정 파라미터 dict 를 만든다. 노트북 get_config() 반환 규격과 동일."""
    ov = dict(overrides or {})
    cfg = {}
    for k in THRESH_KEYS:
        d = CFG_DEFAULT_MAP[k]
        cfg[k] = _coerce(ov[k], d) if k in ov and ov[k] is not None else d
    for k in DESIGN_KEYS:
        d = DESIGN_DEFAULT_MAP[k]
        cfg[k] = _coerce(ov[k], d) if k in ov and ov[k] is not None else d

    src = []
    if cfg["rna_bone_marrow"]:
        src.append("Bone marrow")
    if cfg["rna_peripheral"]:
        src.append("Peripheral leukocytes")
    cfg["rna_source"] = " + ".join(src)
    cfg["batches"] = list(batches or [])
    cfg["nb_version"] = NB_VERSION
    cfg["core_version"] = CORE_VERSION
    cfg["nondefault"] = [k for k in THRESH_KEYS if cfg[k] != CFG_DEFAULT_MAP[k]]
    cfg["param_hash"] = hashlib.md5(
        json.dumps({k: cfg[k] for k in THRESH_KEYS}, sort_keys=True).encode("utf-8")
    ).hexdigest()[:8]

    err = []
    if cfg["insert_min"] >= cfg["insert_max"]:
        err.append("인서트 길이 하한 >= 상한")
    if cfg["d1_min"] >= cfg["d1_max"]:
        err.append("d1 하한 >= 상한")
    if cfg["d2_min"] >= cfg["d2_max"]:
        err.append("d2 하한 >= 상한")
    if cfg["lm_gap_warn"] > cfg["lm_gap_fail"]:
        err.append("랜드마크 갭 WARN > FAIL")
    if cfg["trim_win"] < 1:
        err.append("트리밍 창 크기가 1 미만")
    cfg["errors"] = err
    return cfg


# =============================================================================
#  8. 구조 QC  (노트북 Cell 3 과 동일. CFG 전역 -> cfg 인자)
# =============================================================================
def ungapped_scan(lm, seq):
    """갭 없이 lm 을 seq 전체에 훑어 최소 불일치 수와 위치를 반환."""
    m = len(lm)
    if len(seq) < m:
        return m, -1
    best, bpos = m + 1, -1
    for i in range(len(seq) - m + 1):
        c = 0
        for k in range(m):
            if not hit(lm[k], seq[i + k]):
                c += 1
                if c >= best:
                    break
        if c < best:
            best, bpos = c, i
            if best == 0:
                break
    return best, bpos


def dp_align(lm, win):
    """lm 전체를 win 안에 정렬 (win 양끝 갭 무료). (sub, gap, start) 반환."""
    m, n = len(lm), len(win)
    if n == 0:
        return m, 0, 0
    s_match, s_mis, s_gap = RULES["AL_MATCH"], RULES["AL_MIS"], RULES["AL_GAP"]
    D = [[0] * (n + 1) for _ in range(m + 1)]
    P = [[2] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        D[i][0] = s_gap * i
        P[i][0] = 1
    for i in range(1, m + 1):
        li = lm[i - 1]
        cur, prev, pc = D[i], D[i - 1], P[i]
        for j in range(1, n + 1):
            a = prev[j - 1] + (s_match if hit(li, win[j - 1]) else s_mis)
            b = prev[j] + s_gap
            c = cur[j - 1] + s_gap
            if a >= b and a >= c:
                cur[j], pc[j] = a, 0
            elif b >= c:
                cur[j], pc[j] = b, 1
            else:
                cur[j], pc[j] = c, 2
    jb = 0
    for j in range(n + 1):
        if D[m][j] > D[m][jb]:
            jb = j
    i, j, sub, gap = m, jb, 0, 0
    while i > 0:
        p = P[i][j] if j > 0 else 1
        if p == 0:
            if not hit(lm[i - 1], win[j - 1]):
                sub += 1
            i -= 1
            j -= 1
        elif p == 1:
            gap += 1
            i -= 1
        else:
            gap += 1
            j -= 1
    return sub, gap, j


def check_landmark(lm, seq, qual, covered, cfg):
    """랜드마크 하나를 판정. covered=False 면 read 범위 밖이라 판단 불가(NA)."""
    res = {"sub": 0, "gap": 0, "pos": -1, "minq": None, "status": "NA", "level": "NA"}
    if not covered:
        return res
    p = seq.find(lm)
    if p >= 0:
        res.update(pos=p, status="OK", level="OK")
        return res
    su, pu = ungapped_scan(lm, seq)
    if pu < 0:
        res.update(status="ABSENT", level="FAIL", sub=len(lm))
        return res
    if su <= cfg["lm_max_sub"]:
        res.update(sub=su, pos=pu, status="OK", level="OK")
    else:
        lo = max(0, pu - RULES["AL_FLANK"])
        hi = min(len(seq), pu + len(lm) + RULES["AL_FLANK"])
        sub, gap, st = dp_align(lm, seq[lo:hi])
        res.update(sub=sub, gap=gap, pos=lo + st)
        budget = cfg["lm_max_sub"] + cfg["lm_gap_fail"]
        if gap >= cfg["lm_gap_fail"]:
            res.update(status="GAP" + str(gap), level="FAIL")
        elif sub + gap > budget:
            res.update(status="ABSENT", level="FAIL")
        elif gap >= cfg["lm_gap_warn"] or sub > cfg["lm_max_sub"]:
            res.update(status="S" + str(sub) + "G" + str(gap), level="WARN")
        else:
            res.update(status="OK", level="OK")
    if res["pos"] >= 0 and qual:
        span = qual[res["pos"]:res["pos"] + len(lm)]
        if span:
            res["minq"] = min(span)
    return res


def find_tandem(seq, min_len):
    for period in range(min_len, min(RULES["REPEAT_MAX_PERIOD"], len(seq) // 2) + 1):
        for i in range(0, len(seq) - 2 * period + 1):
            if seq[i:i + period] == seq[i + period:i + 2 * period]:
                return {"pos": i, "period": period, "unit": seq[i:i + period]}
    return None


def trace_purity(trace, ploc, lo, hi, cfg):
    """2순위/1순위 피크비가 임계를 넘는 위치의 비율(%)과 검사 위치 수."""
    keys = sorted(trace.keys())
    if len(keys) != 4 or not ploc:
        return None, None
    chans = [trace[k] for k in keys]
    clen = min(len(c) for c in chans)
    ratios = []
    for i in range(lo, min(hi, len(ploc))):
        p = ploc[i]
        if p < 0 or p >= clen:
            continue
        v = sorted((c[p] for c in chans), reverse=True)
        if v[0] > 0:
            ratios.append(v[1] / float(v[0]))
    if not ratios:
        return None, None
    over = sum(1 for r in ratios if r > cfg["mix_ratio"])
    return 100.0 * over / len(ratios), len(ratios)


def trim_by_quality(seq, qual, q0, win):
    if not qual or len(qual) != len(seq) or len(qual) < win:
        return 0, len(seq)
    lo, hi = 0, len(seq)
    for i in range(len(qual) - win + 1):
        if sum(qual[i:i + win]) / float(win) >= q0:
            lo = i
            break
    for i in range(len(qual) - win, -1, -1):
        if sum(qual[i:i + win]) / float(win) >= q0:
            hi = i + win
            break
    return (lo, hi) if hi > lo else (0, len(seq))


_ANCHORS = ("PELB_ATG", "QC1", "NotI", "QC2", "AscI", "QC3")


def orient(seq):
    fwd = sum(1 for k in _ANCHORS if CONST[k] in seq)
    rev_seq = rc(seq)
    rev = sum(1 for k in _ANCHORS if CONST[k] in rev_seq)
    return ("R", rev_seq) if rev > fwd else ("F", seq)


def qc_one(read, cfg):
    """클론 1개의 구조 QC. read 는 make_read() 가 만든 dict."""
    r = {k: read.get(k) for k in
         ("id", "clone", "batch", "date", "primer", "filename", "raw_len")}
    raw, qraw = read["raw_seq"], read["raw_qual"]

    d, oriented = orient(raw)
    qor = list(reversed(qraw)) if (d == "R" and qraw) else list(qraw)
    lo, hi = trim_by_quality(oriented, qor, cfg["trim_q"], cfg["trim_win"])
    s = oriented[lo:hi]
    q = qor[lo:hi] if qor else []
    r.update(direction=d, trim_lo=lo, trim_hi=hi, used_bp=len(s), seq=s, qual=q)

    n_notI = s.count(CONST["NotI"])
    n_ascI = s.count(CONST["AscI"])
    n_link = s.count(CONST["QC2"])
    pos_n = s.find(CONST["NotI"])
    pos_a = s.find(CONST["AscI"])
    pos_l = s.find(CONST["QC2"])
    r.update(n_notI=n_notI, n_ascI=n_ascI, n_link=n_link,
             pos_notI=pos_n, pos_ascI=pos_a, pos_link=pos_l)

    flags = []
    notes = []

    if max(n_notI, n_ascI, n_link) > 1:
        flags.append("CONCATEMER")
        notes.append("NotI %d / AscI %d / 링커 %d 회 검출 (각 1 회여야 함)"
                     % (n_notI, n_ascI, n_link))

    r["stuffer"] = CONST["STUFFER"] in s
    if r["stuffer"]:
        flags.append("PARENTAL")
        notes.append("모클론 스터퍼 서열 검출 (미절단/단일절단 벡터)")

    q1off = CONST["QC1"].find(CONST["NotI"])
    cov = {"QC1": pos_n >= 0 and (pos_n - q1off) >= 0,
           "QC2": True, "QC3": pos_a >= 0, "QC4": pos_a >= 0}
    if pos_a >= 0:
        cov["QC3"] = pos_a + len(CONST["QC3"]) <= len(s)
        q4s = pos_a + len(CONST["QC3"])
        cov["QC4"] = q4s + len(CONST["QC4"]) <= len(s)
    qc = {}
    for k in ("QC1", "QC2", "QC3", "QC4"):
        qc[k] = check_landmark(CONST[k], s, q, cov[k], cfg)
    r["qc"] = qc

    for k in ("QC1", "QC2", "QC3", "QC4"):
        lv = qc[k]["level"]
        qtail = (", 최저 Q %d" % qc[k]["minq"]) if qc[k]["minq"] is not None else ""
        if lv == "FAIL":
            if k == "QC2":
                f = "NO_LINKER" if qc[k]["status"] == "ABSENT" else "LINKER_DEL"
            else:
                f = "QC_ABSENT" if qc[k]["status"] == "ABSENT" else "QC_DEL"
            if f not in flags:
                flags.append(f)
            notes.append("%s %s (치환 %d, 갭 %d%s)"
                         % (k, qc[k]["status"], qc[k]["sub"], qc[k]["gap"], qtail))
        elif lv == "WARN":
            if "QC_WARN" not in flags:
                flags.append("QC_WARN")
            notes.append("%s %s (치환 %d, 갭 %d%s) - 크로마토그램 확인 권장"
                         % (k, qc[k]["status"], qc[k]["sub"], qc[k]["gap"], qtail))
        elif lv == "NA":
            if "LOW_COVERAGE" not in flags:
                flags.append("LOW_COVERAGE")
            notes.append(k + " 는 read 범위 밖이라 판단 불가")

    insert = d1 = d2 = None
    if pos_n < 0:
        flags.append("NO_NOTI")
        notes.append("NotI 미검출")
    if pos_a < 0:
        if pos_l >= 0 and (len(s) - (pos_l + len(CONST["QC2"]))) >= cfg["d2_max"]:
            flags.append("NO_ASCI")
            notes.append("AscI 미검출 (VL 구간이 충분히 읽혔는데도 없음)")
        else:
            flags.append("LONG_INSERT?")
            notes.append("AscI 미검출 + 3' 커버리지 부족. 인서트가 길거나 read 가 짧음"
                         " -> 역방향 시퀀싱 권고")
    if pos_l < 0 and "NO_LINKER" not in flags:
        flags.append("NO_LINKER")
        notes.append("링커 미검출")

    if pos_n >= 0 and pos_a >= 0 and pos_a > pos_n:
        insert = pos_a + len(CONST["AscI"]) - pos_n
        if insert < cfg["insert_min"]:
            flags.append("TOO_SHORT")
            notes.append("인서트 %d bp (하한 %d 미만)" % (insert, cfg["insert_min"]))
        elif insert > cfg["insert_max"]:
            flags.append("TOO_LONG")
            notes.append("인서트 %d bp (상한 %d 초과)" % (insert, cfg["insert_max"]))
        if insert % 3 != CONST["FRAME_MOD"]:
            flags.append("FRAMESHIFT")
            notes.append("인서트 %d bp, %%3 = %d (유지 조건 %d)"
                         % (insert, insert % 3, CONST["FRAME_MOD"]))
    if pos_n >= 0 and pos_l >= 0 and pos_l > pos_n:
        d1 = pos_l - pos_n
        if not (cfg["d1_min"] <= d1 <= cfg["d1_max"]):
            flags.append("ABERRANT_D1")
            notes.append("d1 %d bp (범위 %d~%d 밖) - VH 구간 이상"
                         % (d1, cfg["d1_min"], cfg["d1_max"]))
    if pos_l >= 0 and pos_a >= 0 and pos_a > pos_l:
        d2 = pos_a - (pos_l + len(CONST["QC2"]))
        if not (cfg["d2_min"] <= d2 <= cfg["d2_max"]):
            flags.append("ABERRANT_D2")
            notes.append("d2 %d bp (범위 %d~%d 밖) - VL 구간 이상"
                         % (d2, cfg["d2_min"], cfg["d2_max"]))
    r.update(insert_bp=insert, d1=d1, d2=d2)

    prot, stop_ok = "", None
    ip = s.find(CONST["PELB_ATG"])
    if ip >= 0:
        orf = s[ip:ip + ((len(s) - ip) // 3) * 3]
        full = translate(orf)
        st = full.find("*")
        prot = full[:st] if st >= 0 else full
        if st >= 0 and pos_a >= 0:
            stop_ok = (ip + st * 3) >= pos_a
        elif st >= 0:
            stop_ok = None
        else:
            stop_ok = True
        if stop_ok is False:
            flags.append("INTERNAL_STOP")
            notes.append("AscI 이전에 종결코돈 (번역 %d aa 에서 중단)" % len(prot))
    r.update(prot=prot, aa_len=len(prot), stop_ok=stop_ok, orf_start=ip)

    rep = None
    for lab, a0, a1 in (("VH", pos_n, pos_l),
                        ("VL", (pos_l + len(CONST["QC2"])) if pos_l >= 0 else -1, pos_a)):
        if a0 >= 0 and a1 > a0:
            hitr = find_tandem(s[a0:a1], cfg["repeat_min_len"])
            if hitr:
                rep = dict(hitr, region=lab)
                break
    r["repeat"] = rep
    if rep:
        flags.append("TANDEM_REPEAT")
        notes.append("%s 구간에 %d nt 탠덤 반복 (%s) - PCR slipped-strand 산물 의심"
                     % (rep["region"], rep["period"], rep["unit"]))

    if d == "R":
        o_lo, o_hi = len(raw) - hi, len(raw) - lo
    else:
        o_lo, o_hi = lo, hi
    mix, npos = trace_purity(read.get("trace", {}), read.get("ploc", []),
                             o_lo, o_hi, cfg)
    r.update(mix_pct=mix, mix_n=npos)
    if mix is not None and mix > cfg["mix_pct"]:
        flags.append("MIXED")
        notes.append("2순위 피크 비율 %.2f 초과 위치가 %.1f%% (임계 %.1f%%) - 혼합 콜로니 의심"
                     % (cfg["mix_ratio"], mix, cfg["mix_pct"]))

    sev_map = RULES["FLAG_SEV"]
    sev = max([sev_map.get(f, 0) for f in flags], default=0)
    if sev >= 3:
        verdict = [f for f in flags if sev_map.get(f, 0) == 3][0]
    elif sev == 2:
        verdict = [f for f in flags if sev_map.get(f, 0) == 2][0]
    elif sev == 1:
        verdict = "WARN"
    else:
        verdict = "PASS"
    if verdict == "PASS" and any(qc[k]["level"] == "NA" for k in qc):
        verdict = "PASS*"
    r.update(flags=flags, notes=notes, severity=sev, verdict=verdict)
    return r


# =============================================================================
#  9. 프라이머 판별  (노트북 Cell 4 와 동일. CFG/PRIMERS 전역 -> 인자)
# =============================================================================
def locate_core(p, r):
    """프라이머 p 의 판별구간이 read 어디에 놓이는지 (start, core) 반환."""
    g, seq, t = p["group"], p["seq"], p["core_trim"]
    if g == "F1_For":
        if r["pos_notI"] < 0:
            return None
        j = seq.find(CONST["NotI"])
        if j < 0:
            return None
        return r["pos_notI"] - j + t, seq[t:]
    if g == "F3_For":
        if r["pos_link"] < 0:
            return None
        j = CONST["QC2"].find(seq[:t]) if t > 0 else -1
        base = r["pos_link"] + j if j >= 0 else r["pos_link"] + len(CONST["QC2"]) - t
        return base + t, seq[t:]
    if g == "F2_Rev":
        if r["pos_link"] < 0:
            return None
        rcs = rc(seq)
        j = rcs.find(CONST["QC2"])
        if j < 0:
            return None
        return r["pos_link"] - j, rcs[:len(seq) - t]
    if g == "F3_Rev":
        if r["pos_ascI"] < 0:
            return None
        rcs = rc(seq)
        j = rcs.find(CONST["AscI"])
        if j < 0:
            return None
        return r["pos_ascI"] - j, rcs[:len(seq) - t]
    return None


def _mm_positions(core, s, st):
    return [i for i in range(len(core)) if not hit(core[i], s[st + i])]


def score_group(r, group, primers, cfg, chain=None):
    """그룹 내 모든 프라이머를 채점. 최소 미스매치 후보 집합과 차순위 간격 반환."""
    s = r["seq"]
    cands = []
    for p in primers:
        if p["group"] != group:
            continue
        if chain is not None and p["chain"] != chain:
            continue
        loc = locate_core(p, r)
        if loc is None:
            continue
        st, core = loc
        if st < 0 or st + len(core) > len(s):
            continue
        mmp = _mm_positions(core, s, st)
        cands.append({"name": p["name"], "families": p["families"], "chain": p["chain"],
                      "mm": len(mmp), "mmpos": mmp, "len": len(core), "start": st})
    if not cands:
        return None
    cands.sort(key=lambda x: (x["mm"], -x["len"]))
    best = cands[0]["mm"]
    top = [c for c in cands if c["mm"] <= best + cfg["primer_margin_min"]]
    nxt = [c for c in cands if c["mm"] > best]
    fams = sorted(set(f for c in top for f in c["families"]))
    return {"mm": best,
            "delta": (nxt[0]["mm"] - best) if nxt else None,
            "names": [c["name"] for c in top],
            "families": fams,
            "chain": top[0]["chain"],
            "mmpos": cands[0]["mmpos"],
            "core_len": cands[0]["len"],
            "start": cands[0]["start"],
            "ok": best <= cfg["primer_max_mismatch"]}


def fmt_call(call, short=False):
    if call is None:
        return "-"
    if not call["ok"]:
        return "X(mm%d)" % call["mm"]
    fam = "|".join(call["families"]) if call["families"] else "?"
    if short:
        return fam
    names = call["names"]
    if len(names) == 1:
        tag = names[0]
    else:
        pres = set(n.rsplit("-", 1)[0] for n in names)
        if len(pres) == 1:
            tag = names[0].rsplit("-", 1)[0] + "-" \
                  + "|".join(n.rsplit("-", 1)[-1] for n in names)
        else:
            tag = "|".join(names)
    return fam + " (" + tag + ")"


def cdr3_boundary(r, vh_family, primers, cfg):
    """CDR3 경계 정합성 검사. 분류가 아니라 하나의 가설만 검정한다."""
    s = r["seq"]
    if r["pos_notI"] < 0 or r["pos_link"] < 0 or not vh_family:
        return {"status": "-"}
    win = s[r["pos_notI"]:r["pos_link"]]
    z = win.rfind(CONST["MOTIF_F2_For"])
    if z < 0:
        return {"status": "앵커없음"}
    z += r["pos_notI"]
    cand = [p for p in primers if p["group"] == "F2_For" and vh_family in p["families"]]
    if not cand:
        return {"status": "프라이머없음"}
    p = cand[0]
    off = p["seq"].find(CONST["MOTIF_F2_For"])
    if off < 0:
        return {"status": "모티프없음"}
    st = z - off
    trim = RULES["EXO_TRIM"]
    a, b = trim, len(p["seq"]) - trim
    if st + a < 0 or st + b > len(s):
        return {"status": "범위밖"}
    mm = sum(1 for i in range(a, b) if not hit(p["seq"][i], s[st + i]))
    return {"status": "OK" if mm <= cfg["primer_max_mismatch"] else "부정합",
            "mm": mm, "primer": p["name"], "cmp": b - a}


_LINKER_AA = translate(CONST["QC2"])          # = GGGGSGGGGSGGGGS


def cdr3_h3(prot):
    if not prot:
        return None
    li = prot.find(_LINKER_AA)
    seg = prot[:li] if li >= 0 else prot
    last = None
    for m in re.finditer(CONST["FR4_MOTIF"], seg):
        last = m
    if last is None:
        return None
    c = seg.rfind("C", 0, last.start())
    if c < 0:
        return None
    return seg[c + 1:last.start()]


def call_one(r, primers, cfg):
    """클론 1개의 프라이머 판별."""
    out = {"id": r["id"], "qc_verdict": r["verdict"]}
    flags = []
    notes = []

    vh = score_group(r, "F1_For", primers, cfg)
    jh = score_group(r, "F2_Rev", primers, cfg)
    out["vh"], out["jh"] = vh, jh
    if vh is None or not vh["ok"]:
        flags.append("NO_FRAG1")
        if vh is not None:
            notes.append("F1_For 최소 미스매치 %d (허용 %d). 최근접 %s, 불일치 위치(core) %s"
                         % (vh["mm"], cfg["primer_max_mismatch"],
                            ",".join(vh["names"]), [x + 1 for x in vh["mmpos"]]))
        else:
            notes.append("F1_For 판별 불가 (NotI 앵커 없음)")

    vl_k = score_group(r, "F3_For", primers, cfg, "kappa")
    vl_l = score_group(r, "F3_For", primers, cfg, "lambda")
    pool = [c for c in (vl_k, vl_l) if c is not None]
    vl = min(pool, key=lambda c: c["mm"]) if pool else None
    tie_chain = (vl_k is not None and vl_l is not None and vl_k["mm"] == vl_l["mm"])
    chain = vl["chain"] if (vl is not None and vl["ok"] and not tie_chain) else None
    out["vl"], out["chain"], out["chain_tie"] = vl, chain, tie_chain
    if vl is None or not vl["ok"]:
        flags.append("NO_VL")
        if vl is not None:
            notes.append("F3_For 최소 미스매치 %d (허용 %d). 최근접 %s"
                         % (vl["mm"], cfg["primer_max_mismatch"], ",".join(vl["names"])))
        else:
            notes.append("F3_For 판별 불가 (링커 앵커 없음)")
    elif tie_chain:
        notes.append("kappa/lambda 미스매치 동점 (%d) - 경쇄 판정 보류" % vl["mm"])

    vj = score_group(r, "F3_Rev", primers, cfg, chain) if chain else None
    out["vj"] = vj

    if vl is not None and not vl["ok"] and r["pos_link"] >= 0:
        s = r["seq"]
        st = r["pos_link"] + len(CONST["QC2"])
        for p in primers:
            if p["group"] != "F1_For":
                continue
            core = p["seq"][p["core_trim"]:]
            if st + len(core) <= len(s) and \
               len(_mm_positions(core, s, st)) <= cfg["primer_max_mismatch"]:
                notes.append("링커 뒤 서열이 VH FR1(%s) 과 일치 - VH-VH 조립 의심" % p["name"])
                break

    fam_matched = (vh["families"] if (vh is not None and vh["ok"]) else [])
    if cfg["batch_vh_family"] != NOSEL and fam_matched:
        if cfg["batch_vh_family"] not in fam_matched:
            flags.append("WRONG_FAMILY")
            notes.append("배치 지정 %s 인데 %s 로 판정 - 튜브 간 교차오염 의심"
                         % (cfg["batch_vh_family"], "|".join(fam_matched)))
    if cfg["batch_chain"] != NOSEL and chain:
        if chain != cfg["batch_chain"]:
            flags.append("WRONG_CHAIN")
            notes.append("배치 지정 %s 인데 %s 로 판정" % (cfg["batch_chain"], chain))

    if cfg["f1_for_mode"] == "분리" and cfg["f2_for_mode"] == "분리":
        vhf = fam_matched[0] if len(fam_matched) == 1 else None
        bd = cdr3_boundary(r, vhf, primers, cfg)
    else:
        bd = {"status": "pool"}
    out["boundary"] = bd
    if bd.get("status") == "부정합":
        flags.append("HETERO_JOIN?")
        notes.append("CDR3 경계가 %s 와 부정합 (미스매치 %d/%d) - frag1/frag2 이종 조립 의심 [약한 증거]"
                     % (bd["primer"], bd["mm"], bd["cmp"]))

    cd = cdr3_h3(r.get("prot", ""))
    out["cdr3"] = cd
    out["cdr3_len"] = len(cd) if cd else None
    out["flags"], out["notes"] = flags, notes
    return out


# =============================================================================
#  10. 배치 조성
# =============================================================================
def median(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def _fmt_med(m):
    if m is None:
        return None
    return int(m) if float(m).is_integer() else round(m, 1)


def compose(primer_results, cfg):
    """구조 QC 와 프라이머 판별을 모두 통과한 클론만 집계한다."""
    good = [p for p in primer_results
            if p["qc_verdict"] in ("PASS", "PASS*") and not p["flags"]]

    def tally(key):
        c = {}
        for g in good:
            call = g[key]
            if call is not None and call["ok"]:
                lab = "|".join(call["families"]) or "?"
                c[lab] = c.get(lab, 0) + 1
        # 노트북 Counter.most_common() 과 동일 : 건수 내림차순, 동점은 삽입(클론) 순
        return sorted(c.items(), key=lambda kv: -kv[1])

    chain_c = {}
    for g in good:
        if g["chain"]:
            chain_c[g["chain"]] = chain_c.get(g["chain"], 0) + 1
    lens = [g["cdr3_len"] for g in good if g["cdr3_len"] is not None]
    seqs = [g["cdr3"] for g in good if g["cdr3"]]
    seen = {}
    for x in seqs:
        seen[x] = seen.get(x, 0) + 1
    dup = [k for k, v in seen.items() if v > 1]

    vh_t = tally("vh")
    res = {
        "n_total": len(primer_results),
        "n_good": len(good),
        "vh": vh_t,
        "jh": tally("jh"),
        "chain": sorted(chain_c.items(), key=lambda kv: -kv[1]),
        "vl": tally("vl"),
        "vj": tally("vj"),
        "cdr3_lens": sorted(lens),
        "cdr3_median": _fmt_med(median(lens)),
        "cdr3_dup": dup,
        "batch_vh_match": None,
        "batch_chain_match": None,
    }
    if cfg["batch_vh_family"] != NOSEL:
        res["batch_vh_match"] = dict(vh_t).get(cfg["batch_vh_family"], 0)
    if cfg["batch_chain"] != NOSEL:
        res["batch_chain_match"] = chain_c.get(cfg["batch_chain"], 0)

    zone = [p["id"] for p in primer_results
            if p["vh"] is not None and p["vh"]["mmpos"]
            and all(x < RULES["OVERLONG_ZONE"] for x in p["vh"]["mmpos"])]
    res["overlong_suspect"] = zone if len(zone) >= 2 else []
    res["overlong_zone"] = RULES["OVERLONG_ZONE"]
    return res


# =============================================================================
#  11. 용어 설명  (xlsx 06_용어설명 시트 내용)
# =============================================================================
def glossary():
    L = CONST["QC2"]
    sev3 = ", ".join(sev_flags(3))
    sev2 = ", ".join(sev_flags(2))
    sev1 = ", ".join(sev_flags(1))
    return [
        ["파일 구성", "01_판정요약", "클론 1행. 이 노트북의 결론. 여기만 봐도 통과/실패를 알 수 있습니다."],
        ["파일 구성", "02_구조QC상세", "랜드마크 정렬 수치와 위치, 플래그, 비고 원문."],
        ["파일 구성", "03_프라이머판별", "그룹별 1행. 후보 프라이머·미스매치·Δ·불일치 위치."],
        ["파일 구성", "04_배치조성", "통과 클론의 family / JH / J 분포와 CDR3-H3 길이."],
        ["파일 구성", "05_실행설정", "판정에 쓰인 모든 수치와 코드에 고정된 상수 전체."],
        ["파일 구성", "06_용어설명", "이 시트."],
        ["파일 구성", "07_서열", "인서트 염기서열과 scFv 아미노산 서열."],

        ["대상 구조", "scFv 카세트",
         "pelB 리더 - NotI - VH - (G4S)3 링커 - VL - AscI - His6 - 3xFLAG - amber TAG - gIII. "
         "pelB 개시코돈부터 gIII 까지가 하나의 ORF 입니다."],
        ["대상 구조", "amber TAG",
         "scFv 와 pIII 사이의 종결코돈. supE 계열 균주(TG1, ER2738, DH5a)에서는 Gln 으로 읽혀 "
         "scFv-pIII 융합단백질이 만들어지고, 비억제 균주에서는 여기서 끊겨 가용성 scFv 가 됩니다."],
        ["대상 구조", "인서트", "NotI 인식서열 첫 염기부터 AscI 인식서열 마지막 염기까지."],
        ["대상 구조", "d1", "NotI 첫 염기 ~ 링커 시작. VH 전체(fragment 1 + fragment 2) 길이."],
        ["대상 구조", "d2", "링커 끝 ~ AscI 첫 염기. VL 전체(fragment 3) 길이."],

        ["기본 개념", "Phred 품질값 (Q)",
         "시퀀서가 각 염기 호출에 부여한 신뢰도. Q = -10 x log10(오류확률). "
         "Q20 = 오류 1/100 (정확도 99%), Q30 = 1/1000, Q40 = 1/10000. "
         "Q 값은 .ab1 파일에만 들어 있고 텍스트 .seq 에는 없습니다."],
        ["기본 개념", "Q 를 어디에 쓰는가",
         "(1) read 양끝 트리밍 - 설정한 창 크기의 이동평균 Q 가 임계 미만인 구간을 잘라냅니다. "
         "(2) 랜드마크가 '실제로 없는 것'인지 '품질이 낮아 못 읽은 것'인지 구분합니다. "
         "(3) 갭이 검출된 랜드마크의 최저 Q 를 함께 보고해 판단 근거를 남깁니다."],
        ["기본 개념", "트리밍구간",
         "원본 read 에서 잘라내고 남긴 구간(0-based, 끝 미포함). 이후 모든 위치는 이 구간 기준입니다."],
        ["기본 개념", "프레임 (%3)",
         "인서트 길이를 3 으로 나눈 나머지. 정상값은 2 입니다. pelB 개시코돈 기준으로 NotI 첫 염기는 "
         "코돈의 1번째 자리에서 시작하고 AscI 마지막 염기는 코돈의 2번째 자리에서 끝나므로, "
         "온전한 코돈 n 개 + 2 염기가 됩니다. 2 가 아니면 뒤에 이어지는 His6 / 3xFLAG / pIII 가 "
         "전부 다른 아미노산으로 번역되어 파지 표면 제시와 태그 검출이 모두 실패합니다."],
        ["기본 개념", "프레임 검사의 힘",
         "VH·VL 가변 구간에는 대조할 레퍼런스가 없어 랜드마크로는 이상을 잡을 수 없습니다. "
         "프레임 검사는 인서트 어디에서 생긴 삽입·결실이든 빠짐없이 잡아냅니다."],
        ["기본 개념", "내부종결",
         "amber TAG 보다 앞쪽에 종결코돈이 있는지. FAIL 이면 scFv 가 중간에서 끊깁니다."],

        ["QC 랜드마크", "무엇인가",
         "프라이머 설계상 모든 클론에서 100% 동일해야 하는 네 구간. 여기가 깨졌다면 조립이나 "
         "클로닝이 잘못된 것입니다. 반대로 VH·VL 가변 구간이 다른 것은 정상이며 그것이 다양성입니다."],
        ["QC 랜드마크", "QC1", "pelB 3' + NotI + 프레임 유지용 A. " + CONST["QC1"] + " (24 nt)"],
        ["QC 랜드마크", "QC2", "(G4S)3 링커 전체. " + L + " (45 nt)"],
        ["QC 랜드마크", "QC3", "AscI + His6. " + CONST["QC3"] + " (33 nt)"],
        ["QC 랜드마크", "QC4", "3xFLAG + amber TAG. " + CONST["QC4"] + " (69 nt)"],
        ["QC 랜드마크", "치환 vs 갭을 나누는 이유",
         "시퀀싱 오독은 길이를 바꾸지 않는 치환이고, 긴 올리고의 n-1 / n-2 산물은 길이를 바꾸는 "
         "갭입니다. 둘을 나누면 잡음과 진짜 결실이 갈립니다. QC2 는 63-mer Rev-2 프라이머가 "
         "만드는 구간이라 올리고 결실이 여기서 드러납니다."],
        ["QC 랜드마크", "상태 표기 OK", "일치. 또는 허용 치환 이내."],
        ["QC 랜드마크", "상태 표기 S#G#", "치환 # 개, 갭 # 개. 허용치를 넘었으나 랜드마크로는 인식됨 (WARN)."],
        ["QC 랜드마크", "상태 표기 GAP#", "갭 # 개로 FAIL 임계 이상. 올리고 결실 의심."],
        ["QC 랜드마크", "상태 표기 ABSENT",
         "치환 + 갭이 허용 예산(허용치환 + 갭FAIL)을 넘어 랜드마크로 인식되지 않음."],
        ["QC 랜드마크", "상태 표기 NA",
         "해당 구간이 read 범위 밖이라 판단 불가. FAIL 이 아닙니다."],

        ["트레이스 순도", "혼합(%)",
         "크로마토그램에서 각 염기 호출 위치의 2순위 피크 세기 / 1순위 피크 세기 비율을 계산해, "
         "설정한 비율을 넘는 위치가 전체의 몇 %인지. 단일 클론에서는 보통 수 % 이내입니다. "
         "임계를 넘으면 콜로니가 단일 클론이 아니거나 인서트가 두 종류 들어간 것을 의심합니다."],
        ["트레이스 순도", "검사위치수", "위 비율을 계산한 염기 위치의 개수."],

        ["탠덤 반복", "무엇인가",
         "VH·VL 구간에서 같은 서열이 곧바로 이어져 반복되는 구조. 항체 V 유전자에는 없는 형태로, "
         "PCR 중 slipped-strand mispairing 으로 생긴 인공 산물의 지표입니다. "
         "설정한 최소 길이 이상이 완전히 동일하게 2 회 이어질 때만 검출합니다."],

        ["프라이머 판별", "판별구간 (core)",
         "같은 그룹·사슬의 프라이머들이 공유하는 5' 불변 태그를 제거하고 남은 부분. "
         "실제 서열에서 계산하므로 프라이머를 추가하면 자동으로 다시 계산됩니다. "
         "예: Rev-2 계열 63-mer 는 앞 48 nt 가 링커 등으로 동일해 뒤 15 nt 만 JH 판별에 쓰입니다."],
        ["프라이머 판별", "미스매치",
         "판별구간과 read 서열이 다른 위치의 개수. 축퇴염기(R, Y, S, W, K, M, B, D, H, V, N)는 "
         "해당 염기 집합 안이면 일치로 셉니다. 프라이머 구간은 프라이머 서열이 그대로 복제된 "
         "자리이므로 오독 외에는 불일치가 나올 이유가 없습니다."],
        ["프라이머 판별", "Δ (차순위 간격)",
         "최상위 후보와 그 다음 후보의 미스매치 차이. 클수록 판정 근거가 두껍습니다. "
         "JH 는 Rev-2-123 과 Rev-2-4 가 판별구간 15 nt 중 단 1곳에서만 갈리므로 Δ 가 항상 1 입니다. "
         "해당 구간 Q 가 50 이상이면 오독 확률이 10^-5 수준이라 1 nt 차이도 신뢰할 만합니다. "
         "판정이 실패(X)한 경우의 Δ 는 실패한 후보들끼리의 순위차일 뿐이므로 읽지 않습니다."],
        ["프라이머 판별", "모호성 표기",
         "IGKV1 (For3-k-5|6) 처럼 표기합니다. 왼쪽 family 는 확정이고, 괄호 안은 축퇴 공간이 "
         "겹쳐 서로 구분할 수 없는 프라이머 후보 집합입니다. 하나로 골라 적는 것은 근거 없는 "
         "정보이므로 하지 않습니다. family 칸의 JH1|JH2 처럼 프라이머 하나가 여러 germline 을 "
         "표적하는 경우도 같은 기호를 쓰지만, 이때 후보 프라이머는 1 개입니다."],
        ["프라이머 판별", "모호성 클러스터",
         "For3-k-1/2/5/6 (전부 IGKV1), For3-k-8/9/12 (전부 IGKV2), For3-L-5/6 (IGLV2), "
         "For3-L-8/9 (IGLV3), For-1-3b/3c (VH3). 다섯 클러스터 모두 같은 family 안에서만 겹치므로 "
         "프라이머 이름이 모호해도 family 는 항상 유일하게 확정됩니다."],
        ["프라이머 판별", "그룹 F1_For", "VH FR1. NotI 를 앵커로 위치를 잡습니다. 판별 신뢰도 높음."],
        ["프라이머 판별", "그룹 F2_Rev", "JH + 링커. 링커 시작을 앵커로 잡는 역방향 프라이머. 신뢰도 높음."],
        ["프라이머 판별", "그룹 F3_For",
         "VL FR1. 링커 끝을 앵커로 잡습니다. kappa / lambda 를 각각 채점해 낮은 쪽을 택합니다."],
        ["프라이머 판별", "그룹 F3_Rev", "VL J. AscI 를 앵커로 잡는 역방향 프라이머. 신뢰도 높음."],
        ["프라이머 판별", "CDR3 경계 (F1_Rev / F2_For)",
         "분류에 쓰지 않고 정합성 검사로만 씁니다. 이유 두 가지 - (1) Rev-1-3 / For-2-3 의 축퇴 "
         "공간이 다른 family 프라이머를 완전히 포함해 유일 배정이 불가능합니다. "
         "(2) overlap extension 중 proofreading 중합효소의 3'->5' exonuclease 가 미스매치된 3' 말단을 "
         "제거하고 상대 fragment 를 주형으로 재연장하기 때문에, 이 구간은 프라이머 서열과 원 주형의 "
         "혼합이 됩니다. 그래서 양끝 3 nt 를 제외하고 채점합니다."],
        ["프라이머 판별", "CDR3경계 값 OK / 부정합 / pool / 앵커없음",
         "OK = FR1 이 부른 family 와 정합. 부정합 = frag1/frag2 이종 조립 의심(약한 증거). "
         "pool = fragment 1·2 를 pool 로 PCR 했다고 설정해 검사를 생략. "
         "앵커없음 = CDR3 경계의 보존 모티프를 찾지 못함."],
        ["프라이머 판별", "Over 프라이머",
         "For-Over / Rev-Over 는 fragment 프라이머의 5' 불변 태그와 동일한 부분집합이라 "
         "최종 산물에 자기만의 흔적을 남기지 않습니다. 따라서 판별 대상에서 제외됩니다."],

        ["CDR3-H3", "정의",
         "VH 의 세 번째 상보성 결정 부위. 보존된 Cys 다음 잔기부터 FR4 의 보존 Trp 앞 잔기까지 "
         "(IMGT 관례). 항체 다양성의 핵심 부위입니다."],
        ["CDR3-H3", "중앙값", "값이 짝수 개면 가운데 두 값의 평균입니다."],

        ["판정 코드", "PASS", "모든 항목 통과."],
        ["판정 코드", "PASS*", "통과. 단 일부 랜드마크가 read 커버리지 부족으로 판단 불가(NA)."],
        ["판정 코드", "WARN", "통과했으나 크로마토그램 육안 확인 권장."],
        ["판정 코드", "LINKER_DEL",
         "(G4S)3 링커에 결실. Rev-2 올리고의 n-1 / n-2 산물 의심. 올리고 정제 등급을 확인하세요."],
        ["판정 코드", "NO_LINKER", "링커 미검출. fragment 2 또는 3 이 빠진 조립 산물."],
        ["판정 코드", "FRAMESHIFT", "인서트 길이 %3 이 2 가 아님. 뒤쪽 태그와 pIII 가 전부 무의미해집니다."],
        ["판정 코드", "INTERNAL_STOP", "amber TAG 이전에 종결코돈."],
        ["판정 코드", "PARENTAL", "모클론 스터퍼 서열 검출. 미절단 또는 단일절단 벡터가 재결합한 배경 클론."],
        ["판정 코드", "CONCATEMER", "NotI / AscI / 링커가 2 회 이상 검출. 인서트가 여러 개 들어감."],
        ["판정 코드", "TOO_SHORT / TOO_LONG", "인서트 길이가 설정 범위 밖."],
        ["판정 코드", "ABERRANT_D1 / ABERRANT_D2", "VH 또는 VL 구간 길이가 설정 범위 밖."],
        ["판정 코드", "TANDEM_REPEAT", "V 구간에 탠덤 반복. PCR slipped-strand 산물 의심."],
        ["판정 코드", "MIXED", "트레이스 혼합. 콜로니를 다시 분리하세요."],
        ["판정 코드", "LONG_INSERT?",
         "AscI 미검출 + 3' 커버리지 부족. 인서트가 read 보다 길 가능성. 역방향 시퀀싱을 권합니다."],
        ["판정 코드", "LOW_COVERAGE", "일부 랜드마크가 read 범위 밖. FAIL 이 아닙니다."],
        ["판정 코드", "NO_FRAG1", "NotI 직후 서열이 어떤 F1_For 프라이머와도 허용 미스매치 안에서 맞지 않음."],
        ["판정 코드", "NO_VL", "링커 직후 서열이 어떤 VL 프라이머와도 맞지 않음."],
        ["판정 코드", "WRONG_FAMILY", "배치에 지정한 VH family 와 다름. 튜브 간 교차오염 신호."],
        ["판정 코드", "WRONG_CHAIN", "배치에 지정한 경쇄와 다름."],
        ["판정 코드", "HETERO_JOIN?",
         "CDR3 경계가 FR1 이 부른 family 와 부정합. fragment 1 / 2 이종 조립 의심 (약한 증거)."],

        ["고정 상수", "두 계층 CONST 와 RULES",
         "코드에 고정된 값은 두 종류다. CONST 는 서열과 프레임 규칙, RULES 는 알고리즘 규칙"
         "(정렬 점수·탐색 범위·플래그 심각도)이다. 둘 다 UI 로 바꿀 수 없고 param_hash 에도 "
         "들어가지 않는다. 05_실행설정 시트에 전체 목록이 기록된다."],
        ["고정 상수", "RULES 가 판정에 미치는 영향",
         "AL_* 는 랜드마크 정렬의 치환·갭 개수를 정하므로 QC1~QC4 의 OK/WARN/FAIL 을 좌우한다. "
         "AL_FLANK 는 DP 창 크기라 이보다 큰 결실은 검출되지 않는다. EXO_TRIM 은 CDR3 경계 "
         "정합성 검사 범위를, REPEAT_MAX_PERIOD 는 탠덤 반복 탐색 상한을 정한다."],
        ["고정 상수", "FLAG_SEV",
         "한 클론에 플래그가 여러 개 붙을 때 화면과 요약에 표시할 verdict 를 고르는 기준이다. "
         "심각도 3(FAIL) 이 있으면 그중 첫 번째, 없으면 2(확인 필요), 그다음 1(WARN), 아무것도 "
         "없으면 PASS 가 된다. 심각도 3: " + sev3 + ". 심각도 2: " + sev2 +
         ". 심각도 1: " + sev1 + "."],
        ["고정 상수", "FR4_MOTIF",
         "VH FR4 의 보존 모티프 " + CONST["FR4_MOTIF"] + " (Trp-Gly-Xxx-Gly). 정규식으로 쓰이며 "
         "CDR3-H3 의 끝 지점을 정한다. 이 모티프를 못 찾으면 CDR3-H3 가 추출되지 않는다."],

        ["재현성", "param_hash",
         "판정 임계값 16 개를 정렬해 만든 지문. 두 배치의 결과를 비교하려면 이 값이 같아야 합니다. "
         "실험 설계와 cDNA / RNA 출처는 해시에 포함되지 않고 05_실행설정 시트에 개별 기록됩니다."],
        ["재현성", "기본값과 동일",
         "05_실행설정 시트의 열. X 로 표시된 항목은 기본값에서 바뀐 것이며 판정 결과에 직접 영향을 줍니다."],
    ]


# =============================================================================
#  12. 결과 표 / 시트 조립
# =============================================================================
SUMMARY_HEADERS = [
    "clone", "batch", "date", "primer", "dir", "최종판정", "구조QC", "프라이머플래그",
    "사용길이(bp)", "인서트(bp)", "프레임(%3)", "d1(bp)", "d2(bp)",
    "QC1", "QC2", "QC3", "QC4", "내부종결", "스터퍼", "탠덤반복", "혼합(%)",
    "VH family", "JH", "경쇄", "VL-V family", "VL-J",
    "CDR3-H3", "CDR3-H3(aa)", "CDR3경계", "param_hash"]

_STOP_TXT = {True: "OK", False: "FAIL", None: "-"}


def _yn(v, yes="YES", no="-"):
    return yes if v else no


def final_verdict(r, p):
    if r["verdict"] not in ("PASS", "PASS*"):
        return "FAIL"
    if p["flags"]:
        return "FAIL"
    return "PASS" if r["verdict"] == "PASS" else "PASS*"


def _fam(call):
    return "|".join(call["families"]) if (call and call["ok"]) else "-"


def build_summary(qc_results, calls, cfg):
    rows = []
    by_id = {p["id"]: p for p in calls}
    for r in qc_results:
        p = by_id[r["id"]]
        ins = r["insert_bp"]
        rows.append([
            r["id"], r["batch"], r["date"], r["primer"], r["direction"],
            final_verdict(r, p), r["verdict"], ", ".join(p["flags"]) or "-",
            r["used_bp"], ins, (ins % 3) if ins is not None else None, r["d1"], r["d2"],
            r["qc"]["QC1"]["status"], r["qc"]["QC2"]["status"],
            r["qc"]["QC3"]["status"], r["qc"]["QC4"]["status"],
            _STOP_TXT[r["stop_ok"]], _yn(r["stuffer"]), _yn(bool(r["repeat"])),
            round(r["mix_pct"], 2) if r["mix_pct"] is not None else None,
            _fam(p["vh"]), _fam(p["jh"]), p["chain"] or "-",
            _fam(p["vl"]), _fam(p["vj"]),
            p["cdr3"] or "-", p["cdr3_len"], p["boundary"].get("status", "-"),
            cfg["param_hash"]])
    return rows


def build_sheets(qc_results, calls, cfg, comp, meta):
    """xlsx 시트 데이터를 순서대로 담은 리스트를 반환한다."""
    by_id = {p["id"]: p for p in calls}
    sheets = []

    # 01 판정요약
    sheets.append({
        "title": "01_판정요약", "headers": SUMMARY_HEADERS,
        "rows": build_summary(qc_results, calls, cfg), "wrap": [7], "maxw": 60,
        "note": "최종판정 = 구조QC 통과 AND 프라이머 판별 플래그 없음. "
                "각 컬럼의 정의는 06_용어설명 시트 참조."})

    # 02 구조QC상세
    h2 = ["clone", "구조QC", "플래그", "원본길이(bp)", "트리밍구간", "사용길이(bp)",
          "NotI위치", "NotI개수", "링커위치", "링커개수", "AscI위치", "AscI개수",
          "인서트(bp)", "프레임(%3)", "d1(bp)", "d2(bp)"]
    for k in ("QC1", "QC2", "QC3", "QC4"):
        h2 += [k + "상태", k + "치환", k + "갭", k + "최저Q"]
    h2 += ["번역길이(aa)", "내부종결", "스터퍼", "탠덤반복구간", "탠덤반복주기(nt)",
           "탠덤반복단위", "혼합(%)", "검사위치수", "비고"]
    r2 = []
    for r in qc_results:
        rep = r["repeat"] or {}
        row = [r["id"], r["verdict"], ", ".join(r["flags"]) or "-",
               r["raw_len"], "%d-%d" % (r["trim_lo"], r["trim_hi"]), r["used_bp"],
               r["pos_notI"] + 1 if r["pos_notI"] >= 0 else None, r["n_notI"],
               r["pos_link"] + 1 if r["pos_link"] >= 0 else None, r["n_link"],
               r["pos_ascI"] + 1 if r["pos_ascI"] >= 0 else None, r["n_ascI"],
               r["insert_bp"],
               (r["insert_bp"] % 3) if r["insert_bp"] is not None else None,
               r["d1"], r["d2"]]
        for k in ("QC1", "QC2", "QC3", "QC4"):
            row += [r["qc"][k]["status"], r["qc"][k]["sub"],
                    r["qc"][k]["gap"], r["qc"][k]["minq"]]
        row += [r["aa_len"], _STOP_TXT[r["stop_ok"]], _yn(r["stuffer"]),
                rep.get("region", "-"), rep.get("period"), rep.get("unit", "-"),
                round(r["mix_pct"], 2) if r["mix_pct"] is not None else None,
                r["mix_n"], " / ".join(r["notes"]) or "-"]
        r2.append(row)
    sheets.append({
        "title": "02_구조QC상세", "headers": h2, "rows": r2,
        "wrap": [len(h2) - 1], "maxw": 60,
        "note": "위치는 트리밍 후 서열 기준 1-based. 랜드마크 상태 표기"
                "(OK / S#G# / GAP# / ABSENT / NA)는 06_용어설명 참조."})

    # 03 프라이머판별
    h3 = ["clone", "그룹", "대상", "판정 family", "후보 프라이머", "후보수",
          "미스매치", "Δ(차순위간격)", "판별구간(nt)", "불일치위치(1-based)",
          "read상 시작위치", "허용 미스매치", "통과"]
    grp = [("vh", "F1_For", "VH family"), ("jh", "F2_Rev", "JH"),
           ("vl", "F3_For", "VL V-gene"), ("vj", "F3_Rev", "VL J-gene")]
    r3 = []
    for p in calls:
        for key, gname, target in grp:
            c = p[key]
            if c is None:
                r3.append([p["id"], gname, target, "-", "-", 0, None, None, None,
                           "-", None, cfg["primer_max_mismatch"], "앵커없음"])
                continue
            r3.append([p["id"], gname, target,
                       "|".join(c["families"]) or "?", "|".join(c["names"]),
                       len(c["names"]), c["mm"], c["delta"], c["core_len"],
                       ", ".join(str(x + 1) for x in c["mmpos"]) or "-",
                       c["start"] + 1, cfg["primer_max_mismatch"],
                       "OK" if c["ok"] else "FAIL"])
        b = p["boundary"]
        r3.append([p["id"], "F2_For", "CDR3 경계 정합성", "-", b.get("primer", "-"), 0,
                   b.get("mm"), None, b.get("cmp"), "-", None,
                   cfg["primer_max_mismatch"], b.get("status", "-")])
    sheets.append({
        "title": "03_프라이머판별", "headers": h3, "rows": r3, "wrap": [], "maxw": 60,
        "note": "그룹별 1행. Δ 는 최상위와 차순위 후보의 미스매치 간격이며 "
                "클수록 판정 근거가 두껍습니다."})

    # 04 배치조성
    h4 = ["구분", "항목", "값", "비고"]
    r4 = [["개요", "배치", ", ".join(cfg["batches"]) or "-", ""],
          ["개요", "분석 클론 수", comp["n_total"], ""],
          ["개요", "구조QC + 프라이머판별 모두 통과", comp["n_good"],
           "아래 분포는 이 클론들만 집계"]]
    if comp["batch_vh_match"] is not None:
        r4.append(["대조", "지정 VH family 일치",
                   "%d / %d" % (comp["batch_vh_match"], comp["n_good"]),
                   "지정 %s. 불일치는 튜브 간 교차오염 신호" % cfg["batch_vh_family"]])
    if comp["batch_chain_match"] is not None:
        r4.append(["대조", "지정 경쇄 일치",
                   "%d / %d" % (comp["batch_chain_match"], comp["n_good"]),
                   "지정 %s" % cfg["batch_chain"]])
    for lab, key, memo in (("VH family", "vh", "배치별 고정. 다양성 지표 아님"),
                           ("JH", "jh", "Rev-2 pool -> 배치 내 자유 변수"),
                           ("경쇄", "chain", "배치별 고정"),
                           ("VL V-gene", "vl", "frag3 pool -> 배치 내 자유 변수"),
                           ("VL J-gene", "vj", "Rev3 pool -> 배치 내 자유 변수")):
        items = comp[key]
        if items:
            for k, v in items:
                r4.append(["분포", lab, "%s : %d" % (k, v), memo])
            r4.append(["분포", lab + " (종 수)", len(items), ""])
        else:
            r4.append(["분포", lab, "-", memo])
    if comp["cdr3_lens"]:
        r4.append(["CDR3-H3", "길이 목록 (aa)",
                   ", ".join(str(x) for x in comp["cdr3_lens"]), ""])
        r4.append(["CDR3-H3", "중앙값 (aa)", comp["cdr3_median"], ""])
    r4.append(["CDR3-H3", "중복 서열 수", len(comp["cdr3_dup"]),
               ", ".join(comp["cdr3_dup"]) if comp["cdr3_dup"]
               else "동일 CDR3-H3 를 가진 클론 없음"])
    if comp["overlong_suspect"]:
        r4.append(["점검", "F1_For 불일치 앞쪽 편중",
                   ", ".join(comp["overlong_suspect"]),
                   "판별구간 앞 %d nt 안에만 불일치가 몰린 클론. "
                   "For-Over-long 프라이머 3' 말단 설계 영향 가능성"
                   % comp["overlong_zone"]])
    sheets.append({
        "title": "04_배치조성", "headers": h4, "rows": r4, "wrap": [3], "maxw": 60,
        "note": "구조QC 와 프라이머 판별을 모두 통과한 클론만 집계합니다."})

    # 05 실행설정
    h5 = ["구분", "항목", "값", "기본값", "기본값과 동일", "설명"]
    r5 = [["실행", "core 버전", cfg["core_version"], "", "", "기준 노트북 v" + cfg["nb_version"]],
          ["실행", "실행 환경", meta.get("runtime", ""), "", "", ""],
          ["실행", "실행 시각", meta.get("timestamp", ""), "", "", ""],
          ["실행", "param_hash", cfg["param_hash"], "", "",
           "판정 임계값 16개의 지문. 배치 간 비교 시 이 값이 같아야 동일 기준"],
          ["실행", "프라이머 FASTA", meta.get("primer_file", ""), "", "",
           "%d 종" % meta.get("primer_n", 0)],
          ["실행", "입력 파일", ", ".join(meta.get("files", [])), "", "", ""]]
    for k, lab, memo in DESIGN_DOC:
        if k in RNA_KEYS:
            continue          # 아래 "RNA 출처" 한 행이 두 값을 합쳐 보여준다
        r5.append(["실험 설계", lab, cfg[k], "", "", memo])
    r5.append(["실험 설계", "RNA 출처", cfg["rna_source"] or "(미지정)", "", "",
               "판정에 미사용. 기록용. 선택한 RNA 출처를 합친 값"])
    for k in THRESH_KEYS:
        r5.append(["판정 임계값", CFG_DOC[k][0], cfg[k], CFG_DEFAULT_MAP[k],
                   "O" if cfg[k] == CFG_DEFAULT_MAP[k] else "X  <-- 변경됨",
                   CFG_DOC[k][1]])
    for k, desc in CONST_DOC:
        r5.append(["고정 상수 (서열)", k, str(CONST[k]), "", "(코드 고정)", desc])
    for k, desc in RULES_DOC:
        r5.append(["고정 상수 (알고리즘)", k, rule_value_text(k), "", "(코드 고정)", desc])
    sheets.append({
        "title": "05_실행설정", "headers": h5, "rows": r5, "wrap": [5], "maxw": 60,
        "note": "이 배치의 판정에 실제로 쓰인 값 전부입니다. "
                "'기본값과 동일' 이 X 인 항목은 결과에 직접 영향을 줍니다."})

    # 06 용어설명
    sheets.append({
        "title": "06_용어설명", "headers": ["구분", "용어 / 컬럼", "설명"],
        "rows": glossary(), "wrap": [2], "maxw": 110,
        "note": "이 도구를 처음 보는 사람도 이 시트만으로 모든 컬럼과 판정을 "
                "이해할 수 있도록 정리했습니다."})

    # 07 서열
    h7 = ["clone", "최종판정", "인서트 길이(bp)", "인서트 염기서열 (NotI~AscI)",
          "scFv 길이(aa)", "scFv 아미노산 서열", "CDR3-H3"]
    r7 = []
    for r in qc_results:
        p = by_id[r["id"]]
        r7.append([r["id"], final_verdict(r, p), r["insert_bp"], insert_seq(r),
                   r["aa_len"], r["prot"], p["cdr3"] or "-"])
    sheets.append({
        "title": "07_서열", "headers": h7, "rows": r7, "wrap": [3, 5], "maxw": 60,
        "note": "인서트는 NotI 인식서열 첫 염기부터 AscI 인식서열 마지막 염기까지입니다."})
    return sheets


def insert_seq(r):
    if r["pos_notI"] >= 0 and r["pos_ascI"] > r["pos_notI"]:
        return r["seq"][r["pos_notI"]:r["pos_ascI"] + len(CONST["AscI"])]
    return ""


def build_fasta(qc_results, calls, line_width=60):
    """최종판정을 통과한 클론의 인서트 염기서열 FASTA."""
    by_id = {p["id"]: p for p in calls}
    out = []
    for r in qc_results:
        p = by_id[r["id"]]
        fin = final_verdict(r, p)
        s = insert_seq(r)
        if fin in ("PASS", "PASS*") and s:
            out.append(">" + r["id"] + " " + (r["batch"] or "-") + " " + fin
                       + " len=" + str(len(s)))
            for i in range(0, len(s), line_width):
                out.append(s[i:i + line_width])
    return "\n".join(out) + ("\n" if out else "")


# =============================================================================
#  13. 최상위 진입점
# =============================================================================
def make_read(filename, data):
    """(파일명, bytes) 에서 분석용 read dict 를 만든다."""
    info = parse_read_name(filename)
    ab = parse_ab1(data)
    info["filename"] = _basename(filename)
    info["raw_seq"] = ab["seq"]
    info["raw_qual"] = ab["qual"]
    info["trace"] = ab["trace"]
    info["ploc"] = ab["ploc"]
    info["raw_len"] = len(ab["seq"])
    info["mean_q"] = (sum(ab["qual"]) / len(ab["qual"])) if ab["qual"] else 0.0
    info["has_qual"] = (len(ab["qual"]) == len(ab["seq"])) and len(ab["qual"]) > 0
    info["has_trace"] = (len(ab["trace"]) == 4) and (len(ab["ploc"]) >= len(ab["seq"]))
    return info


def _load_reads(files):
    reads, errors = [], []
    for name, data in files:
        try:
            reads.append(make_read(name, data))
        except Exception as e:
            errors.append(_basename(name) + " : " + str(e))
    reads.sort(key=lambda r: r["filename"])
    return reads, errors


def _read_brief(r):
    return {"id": r["id"], "filename": r["filename"], "clone": r["clone"],
            "batch": r["batch"], "date": r["date"], "primer": r["primer"],
            "direction": r["direction"], "dir_guessed": r["dir_guessed"],
            "parsed": r["parsed"], "raw_len": r["raw_len"],
            "mean_q": round(r["mean_q"], 1),
            "has_qual": r["has_qual"], "has_trace": r["has_trace"]}


def scan_inputs(files, primer_text=""):
    """폼을 채우기 위한 입력 요약. 판정은 하지 않는다."""
    reads, errors = _load_reads(files)
    primers, trims, pwarns = parse_primer_fasta(primer_text) if primer_text else ([], {}, [])
    batches = sorted(set(r["batch"] for r in reads if r["batch"]))
    groups = {}
    for p in primers:
        key = p["group"] + " / " + p["chain"]
        g = groups.setdefault(key, {"group": p["group"], "chain": p["chain"],
                                    "n": 0, "trim": p["core_trim"],
                                    "core_min": 10 ** 6, "core_max": 0})
        g["n"] += 1
        g["core_min"] = min(g["core_min"], len(p["core"]))
        g["core_max"] = max(g["core_max"], len(p["core"]))
    return {
        "reads": [_read_brief(r) for r in reads],
        "read_errors": errors,
        "batches": batches,
        "primer_n": len(primers),
        "primer_groups": sorted(groups.values(), key=lambda g: (g["group"], g["chain"])),
        "primer_warnings": pwarns,
        "primer_unknown": [p["name"] for p in primers if p["group"] == "UNKNOWN"],
        "options": {"vh_family": [NOSEL] + vh_family_options(primers),
                    "chain": [NOSEL] + chain_options(primers),
                    "pool": list(POOL_OPTS), "cdna": list(CDNA_OPTS)},
        "guess": guess_batch(batches, primers),
        "const": [{"key": k, "value": str(CONST[k]), "desc": d} for k, d in CONST_DOC],
        "cfg_defaults": [{"key": k, "label": CFG_DOC[k][0], "desc": CFG_DOC[k][1],
                          "default": CFG_DEFAULT_MAP[k],
                          "type": "float" if isinstance(CFG_DEFAULT_MAP[k], float) else "int"}
                         for k in THRESH_KEYS],
        "design_doc": [{"key": k, "label": lab, "desc": memo} for k, lab, memo in DESIGN_DOC],
    }


def analyze(files, primer_text="", overrides=None, meta=None):
    """전체 분석. JSON 직렬화 가능한 dict 를 반환한다."""
    meta = dict(meta or {})
    reads, read_errors = _load_reads(files)
    if not reads:
        return {"ok": False, "errors": ["분석할 .ab1 이 없습니다."] + read_errors,
                "read_errors": read_errors}
    primers, trims, pwarns = parse_primer_fasta(primer_text) if primer_text else ([], {}, [])
    batches = sorted(set(r["batch"] for r in reads if r["batch"]))
    cfg = build_config(overrides, batches)
    if cfg["errors"]:
        return {"ok": False, "errors": cfg["errors"], "config": cfg}

    qc_results = [qc_one(r, cfg) for r in reads]
    if primers:
        calls = [call_one(r, primers, cfg) for r in qc_results]
    else:
        calls = [{"id": r["id"], "qc_verdict": r["verdict"], "vh": None, "jh": None,
                  "vl": None, "vj": None, "chain": None, "chain_tie": False,
                  "boundary": {"status": "-"}, "cdr3": None, "cdr3_len": None,
                  "flags": [], "notes": ["프라이머 FASTA 가 없어 판별을 건너뜀"]}
                 for r in qc_results]
    comp = compose(calls, cfg)

    meta.setdefault("primer_file", meta.get("primer_file", ""))
    meta["primer_n"] = len(primers)
    meta["files"] = [r["filename"] for r in reads]
    sheets = build_sheets(qc_results, calls, cfg, comp, meta)
    summary_rows = sheets[0]["rows"]

    pub_qc = []
    for r in qc_results:
        pub_qc.append({k: r[k] for k in
                       ("id", "clone", "batch", "date", "primer", "filename",
                        "direction", "raw_len", "used_bp", "trim_lo", "trim_hi",
                        "n_notI", "n_ascI", "n_link", "pos_notI", "pos_ascI", "pos_link",
                        "insert_bp", "d1", "d2", "aa_len", "stop_ok", "stuffer",
                        "mix_pct", "mix_n", "flags", "notes", "severity", "verdict")}
                      )
        pub_qc[-1]["qc"] = r["qc"]
        pub_qc[-1]["repeat"] = r["repeat"]
    pub_calls = []
    for p in calls:
        d = {k: p[k] for k in ("id", "qc_verdict", "chain", "chain_tie",
                               "boundary", "cdr3", "cdr3_len", "flags", "notes")}
        for k in ("vh", "jh", "vl", "vj"):
            d[k] = p[k]
            d[k + "_txt"] = fmt_call(p[k], short=(k in ("jh", "vj")))
        pub_calls.append(d)

    return {
        "ok": True,
        "errors": [],
        "read_errors": read_errors,
        "primer_warnings": pwarns,
        "primer_trims": trims,
        "version": {"core": CORE_VERSION, "notebook": NB_VERSION},
        "config": cfg,
        "reads": [_read_brief(r) for r in reads],
        "qc": pub_qc,
        "calls": pub_calls,
        "composition": comp,
        "summary": {"headers": SUMMARY_HEADERS, "rows": summary_rows},
        "sheets": sheets,
        "fasta": build_fasta(qc_results, calls),
    }
