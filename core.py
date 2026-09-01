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

CORE_VERSION = "3.2"
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
    ("STUFFER_INSERT_BP",
     "스터퍼 보유 시 인서트 길이 (bp). 서열 미검출 시 이 길이로 PARENTAL? 판정"),
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
    "LM_LOWQ_MARK": 20,
    "COVERAGE_ALPHA": 0.05,
    "FLAG_SEV": {
        "CONCATEMER": 3, "PARENTAL": 3, "NO_NOTI": 3, "NO_ASCI": 3, "NO_LINKER": 3,
        "TOO_SHORT": 3, "TOO_LONG": 3, "FRAMESHIFT": 3, "INTERNAL_STOP": 3,
        "LINKER_DEL": 3, "QC_DEL": 3, "QC_ABSENT": 3,
        "ABERRANT_D1": 3, "ABERRANT_D2": 3, "TANDEM_REPEAT": 3, "MIXED": 3,
        "LONG_INSERT?": 2, "LOW_COVERAGE": 2, "PARENTAL?": 2,
        "AMBIG_FAMILY?": 2, "AMBIG_CALL?": 2,
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
    ("LM_LOWQ_MARK", "랜드마크 이상을 저품질로 표시하는 최저 Q 기준. "
                     "표시 전용이며 판정에는 쓰이지 않는다"),
    ("COVERAGE_ALPHA", "미관측 프라이머를 dropout 이라 말할 수 없다고 보는 확률 기준. "
                       "균등 사용 가정에서 (1-1/S)^n 이 이 값 이상이면 표본 부족. "
                       "표시 전용이며 판정에는 쓰이지 않는다"),
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

# 분석 모드. 라이브러리를 만드는 중에는 VH·VL 을 아는 배치를 검증하고(assigned),
# 완성 후에는 지정 없이 라이브러리 전체를 판정합니다(library).
MODE_ASSIGNED = "assigned"
MODE_LIBRARY = "library"
MODE_NEGCTRL = "negctrl"
# 앞 두 값의 문자열과 순서는 바꾸지 않습니다. 바꾸면 기존 design_hash 가 깨집니다.
MODE_OPTS = [MODE_ASSIGNED, MODE_LIBRARY, MODE_NEGCTRL]

DESIGN_DEFAULTS = [
    ("analysis_mode", MODE_ASSIGNED),
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

# 설계 항목 중 판정 분기를 실제로 만드는 것만 모읍니다. design_hash 의 대상입니다.
#   analysis_mode
#       assigned 일 때만 배치 지정 VH family / 경쇄와 대조한다. library 는 건너뛴다.
#   batch_vh_family / batch_chain
#       call_one 에서 WRONG_FAMILY / WRONG_CHAIN 플래그를 만들고, compose 의
#       batch_vh_match / batch_chain_match 집계와 04_배치조성 "대조" 행까지 좌우한다.
#   f1_for_mode / f2_for_mode
#       둘 다 "분리" 일 때만 CDR3 경계 정합성 검사를 돌린다 (HETERO_JOIN?).
#
# 제외한 설계 항목은 두 부류이며, 둘 다 05_실행설정에 개별 기록됩니다.
#   cfg 에서 아예 읽히지 않음 :
#       f1_rev_mode, f2_rev_mode, f3_for_mode, f3_rev_mode,
#       f3_product_pooled, cdna_frag1~3
#   읽히지만 판정 분기가 아님 :
#       rna_bone_marrow, rna_peripheral — build_config 에서 rna_source
#       표시 문자열을 조립할 때만 쓰인다.
JUDGMENT_DESIGN_KEYS = ["analysis_mode", "batch_vh_family", "batch_chain",
                        "f1_for_mode", "f2_for_mode"]

DESIGN_DOC = [
    ("analysis_mode", "분석 모드",
     "assigned = VH family 와 경쇄를 지정하고 대조한다. "
     "library = 지정 없이 프라이머로만 판정한다. "
     "negctrl = 벡터만 ligation 한 음성 대조군. 기대값이 정반대라 판정 어휘가 다르다"),
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
#  5. 파일명 -> 클론 ID
# =============================================================================
def _basename(path):
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def read_id_from_name(filename):
    """확장자만 뗀 파일명을 그대로 클론 ID 로 씁니다.

    라이브러리 완성 후에는 수백 개를 시퀀싱하므로 파일명을 규칙에 맞춰
    바꾸는 것이 불가능합니다. 그래서 파일명에서 정보를 추출하지 않습니다.
      배치명 / 날짜 : 폼에서 받아 meta["batch_label"] / meta["batch_date"] 로 주입
      판독 방향     : qc_one 의 orient() 가 서열의 앵커 개수로 판정
    """
    stem = _basename(filename)
    if "." in stem:
        stem = stem.rsplit(".", 1)[0]
    return {"id": stem, "clone": stem, "stem": stem,
            "date": "", "batch": "", "primer": ""}


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
    # 판정 임계값과 성격이 달라 param_hash 에 합치지 않고 따로 둡니다.
    # 임계값은 배치 간에 같아야 하고, 설계는 배치마다 달라야 정상입니다.
    cfg["design_hash"] = hashlib.md5(
        json.dumps({k: cfg[k] for k in JUDGMENT_DESIGN_KEYS},
                   sort_keys=True).encode("utf-8")
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

    # 링커 위치 보정 : 완전 일치 탐색(s.find)이 놓쳤어도 랜드마크 검사가 허용
    # 치환 안에서 찾았다면 그 위치를 씁니다. 상태가 OK 일 때만 씁니다 —
    # GAP / ABSENT / S#G# 는 랜드마크가 깨진 것이라 위치를 신뢰할 수 없습니다.
    #
    # NotI 와 AscI 에는 같은 보정을 할 수 없습니다. QC1 은 pos_notI 가, QC3·QC4 는
    # pos_ascI 가 있어야 covered 로 검사되므로(위 cov), 그 자리가 비면 랜드마크도
    # NA 가 되어 되살릴 근거가 없습니다. 즉 이 어긋남은 QC2 에서만 생깁니다.
    if pos_l < 0 and qc["QC2"]["status"] == "OK" and qc["QC2"]["pos"] >= 0:
        pos_l = qc["QC2"]["pos"]
        r["pos_link"] = pos_l
        notes.append("링커를 완전 일치로는 못 찾았으나 랜드마크 검사가 치환 %d 개로 "
                     "위치 %d 에서 찾음 - 그 위치로 d1/d2 를 계산합니다"
                     % (qc["QC2"]["sub"], pos_l + 1))

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
        # 스터퍼 서열은 read 앞쪽이라 품질 저하로 트리밍되면 놓칠 수 있습니다.
        # 인서트 길이는 NotI / AscI 만 잡히면 계산되므로 더 넓게 덮습니다.
        # 다만 서열 근거가 없으므로 확정(PARENTAL)이 아니라 확인 필요 등급으로 둡니다.
        # 스터퍼는 벡터에 고정된 서열이라 길이가 정확하므로 허용 오차를 두지 않습니다.
        if not r["stuffer"] and insert == CONST["STUFFER_INSERT_BP"]:
            flags.append("PARENTAL?")
            notes.append("인서트 길이가 스터퍼 보유 클론과 같음 (%d bp). 스터퍼 서열은 "
                         "검출되지 않았으므로 확정이 아님 - read 앞쪽 품질 저하로 놓쳤을 "
                         "수 있고, 우연히 같은 길이인 다른 산물일 수도 있음. "
                         "크로마토그램 확인 권장" % CONST["STUFFER_INSERT_BP"])
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


def ambiguity_context(primers, cfg):
    """analyze 시작 시 한 번 계산해 call_one 에 넘기는 조회용 묶음.

    클론마다 다시 계산하면 프라이머 수의 제곱에 비례해 낭비입니다.
    """
    pairs = primer_ambiguity(primers, cfg) if primers else []
    return {"pairs": pairs,
            "by_name": dict((frozenset((p["a"], p["b"])), p) for p in pairs)}


def _amb_brief(p):
    return {"a": p["a"], "b": p["b"], "k": p["incompatible"],
            "tie": p["tie"], "families": "|".join(p["families"])}


def family_ambiguity(amb, group, fam_a, fam_b):
    """두 family 를 함께 아우르는 모호성 쌍. 없으면 빈 리스트."""
    out = []
    for p in (amb or {}).get("pairs", []):
        if p["group"] != group or p["same_family"]:
            continue
        fams = set(p["families"])
        if fam_a in fams and fam_b in fams:
            out.append(p)
    return out


def score_group(r, group, primers, cfg, chain=None, amb=None):
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

    # 판정 자체는 위에서 이미 끝났습니다. 아래는 그 판정이 알려진 프라이머
    # 모호성 때문인지 알려주는 근거이며 판정값을 바꾸지 않습니다.
    idx = (amb or {}).get("by_name", {})
    names = [c["name"] for c in top]
    seen, in_top = set(), []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            p = idx.get(frozenset((names[i], names[j])))
            if p is not None and not p["same_family"]:
                key = (p["a"], p["b"])
                if key not in seen:
                    seen.add(key)
                    in_top.append(_amb_brief(p))
    vs_next, seen2 = [], set()
    if nxt:
        lo = nxt[0]["mm"]
        for a in names:
            for c in nxt:
                if c["mm"] != lo:
                    continue
                p = idx.get(frozenset((a, c["name"])))
                if p is not None and (p["a"], p["b"]) not in seen2:
                    seen2.add((p["a"], p["b"]))
                    vs_next.append(_amb_brief(p))

    return {"mm": best,
            "delta": (nxt[0]["mm"] - best) if nxt else None,
            "names": names,
            "families": fams,
            "chain": top[0]["chain"],
            "mmpos": cands[0]["mmpos"],
            "core_len": cands[0]["len"],
            "start": cands[0]["start"],
            "ok": best <= cfg["primer_max_mismatch"],
            "ambiguity": in_top,
            "runner_up_pairs": vs_next}


# call_one 이 실제로 채점하는 그룹. 나머지 두 그룹의 모호성은 판정을 가르지 않습니다.
SCORED_GROUPS = ("F1_For", "F2_Rev", "F3_For", "F3_Rev")


def primer_ambiguity(primers, cfg):
    """판별구간이 서로 겹쳐 동점이 날 수 있는 프라이머 쌍을 프라이머 FASTA 에서 계산한다.

    같은 group·chain 안의 모든 쌍에 대해, 두 판별구간의 축퇴 집합이 서로 겹치지
    않는 위치("비호환")를 셉니다. 허용 미스매치를 T 라 할 때

      비호환 <= T        certain : 한쪽 프라이머와 같은 read 가 이미 두 프라이머
                                   모두와 T 안에서 맞으므로 동점이 확실히 가능
      T < 비호환 <= 2T   split   : read 가 비호환 위치를 두 프라이머에 나눠
                                   부담하면 동점이 가능 (각각 최대 T 까지)

    2T 를 넘으면 어떤 read 도 두 프라이머 모두와 T 안에서 맞을 수 없습니다.
    길이가 다른 판별구간은 짧은 쪽까지만 비교하고 truncated 로 표시합니다.
    """
    tol = cfg["primer_max_mismatch"]
    buckets = {}
    for p in primers:
        if p["group"] in ANALYSIS_GROUPS:
            buckets.setdefault((p["group"], p["chain"]), []).append(p)
    out = []
    for (group, chain), plist in buckets.items():
        for i in range(len(plist)):
            for j in range(i + 1, len(plist)):
                a, b = plist[i], plist[j]
                ca, cb = a["core"], b["core"]
                n = min(len(ca), len(cb))
                bad = sum(1 for k in range(n)
                          if not (set(IUPAC.get(ca[k], ca[k]))
                                  & set(IUPAC.get(cb[k], cb[k]))))
                if bad > 2 * tol:
                    continue
                fams = sorted(set(a["families"]) | set(b["families"]))
                out.append({
                    "group": group, "chain": chain,
                    "a": a["name"], "b": b["name"],
                    "families_a": list(a["families"]), "families_b": list(b["families"]),
                    "families": fams,
                    # 같은 family 면 이름만 모호하고 family 는 확정됩니다.
                    # 다르면 family 자체가 모호해져 판정의 뜻이 달라집니다.
                    "same_family": len(fams) == 1,
                    "cmp_len": n, "len_a": len(ca), "len_b": len(cb),
                    "truncated": len(ca) != len(cb),
                    "incompatible": bad,
                    "tie": "certain" if bad <= tol else "split",
                    "scored": group in SCORED_GROUPS,
                })
    out.sort(key=lambda d: (d["tie"] != "certain", d["same_family"],
                            d["group"], d["chain"], d["a"], d["b"]))
    return out


def ambiguity_summary(pairs):
    """primer_ambiguity 결과를 시트 한 줄로 줄인 요약."""
    cert = [p for p in pairs if p["tie"] == "certain"]
    split = [p for p in pairs if p["tie"] == "split"]
    cross = [p for p in pairs if not p["same_family"]]
    return {"total": len(pairs), "certain": len(cert), "split": len(split),
            "same_family": len(pairs) - len(cross), "cross_family": len(cross),
            "cross_scored": [p for p in cross if p["scored"]]}


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


def call_one(r, primers, cfg, amb=None):
    """클론 1개의 프라이머 판별. amb 는 analyze 가 한 번 계산해 넘기는 모호성 정보."""
    out = {"id": r["id"], "qc_verdict": r["verdict"]}
    flags = []
    notes = []

    vh = score_group(r, "F1_For", primers, cfg, None, amb)
    jh = score_group(r, "F2_Rev", primers, cfg, None, amb)
    out["vh"], out["jh"] = vh, jh
    if vh is None or not vh["ok"]:
        flags.append("NO_FRAG1")
        if vh is not None:
            notes.append("F1_For 최소 미스매치 %d (허용 %d). 최근접 %s, 불일치 위치(core) %s"
                         % (vh["mm"], cfg["primer_max_mismatch"],
                            ",".join(vh["names"]), [x + 1 for x in vh["mmpos"]]))
        else:
            notes.append("F1_For 판별 불가 (NotI 앵커 없음)")

    vl_k = score_group(r, "F3_For", primers, cfg, "kappa", amb)
    vl_l = score_group(r, "F3_For", primers, cfg, "lambda", amb)
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

    vj = score_group(r, "F3_Rev", primers, cfg, chain, amb) if chain else None
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

    # library 모드는 배치에 무엇이 들어 있는지 모르는 상태로 판정하므로 대조하지 않습니다.
    assigned = cfg["analysis_mode"] == MODE_ASSIGNED
    fam_matched = (vh["families"] if (vh is not None and vh["ok"]) else [])
    if assigned and cfg["batch_vh_family"] != NOSEL and fam_matched:
        if cfg["batch_vh_family"] not in fam_matched:
            # 배치 지정 family 와 판정 family 가 알려진 모호성 쌍으로 이어져 있으면
            # 교차오염으로 단정할 수 없습니다. 플래그를 갈라 구분합니다.
            hits = []
            for f in fam_matched:
                hits += family_ambiguity(amb, "F1_For", cfg["batch_vh_family"], f)
            if hits:
                flags.append("AMBIG_FAMILY?")
                notes.append("배치 지정 %s 인데 %s 로 판정. 두 family 는 %s 로 "
                             "구분되지 않으므로 교차오염으로 단정할 수 없음 - "
                             "03_프라이머판별의 모호성 열 확인"
                             % (cfg["batch_vh_family"], "|".join(fam_matched),
                                ", ".join("%s/%s(%s, k%d)"
                                          % (h["a"], h["b"], h["tie"], h["incompatible"])
                                          for h in hits)))
            else:
                flags.append("WRONG_FAMILY")
                notes.append("배치 지정 %s 인데 %s 로 판정 - 튜브 간 교차오염 의심. "
                             "두 family 를 잇는 모호성 쌍은 없음"
                             % (cfg["batch_vh_family"], "|".join(fam_matched)))
    # 판정 후보가 여럿이고 family 도 여러 종이면 family 가 확정되지 않은 것입니다.
    # 배치 지정 여부와 무관하게 붙으므로, 배치 지정이 판정 집합 안에 있어
    # WRONG_FAMILY / AMBIG_FAMILY? 분기에 닿지 않던 경우까지 덮습니다 (이슈 14).
    amb_calls = []
    for _k, _g, _lab in CALLED_GROUPS:
        c = out[_k]
        if c is not None and c["ok"] and len(c["names"]) > 1 and len(c["families"]) > 1:
            amb_calls.append("%s %s (후보 %d : %s)"
                             % (_lab, "|".join(c["families"]), len(c["names"]),
                                ", ".join(c["names"])))
    if amb_calls:
        flags.append("AMBIG_CALL?")
        notes.append("family 가 확정되지 않은 판정 - " + " · ".join(amb_calls) +
                     ". 어느 프라이머였는지 모르므로 family 를 하나로 좁힐 수 없습니다. "
                     "03_프라이머판별의 모호성 열을 보세요")

    if assigned and cfg["batch_chain"] != NOSEL and chain:
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
    assigned = cfg["analysis_mode"] == MODE_ASSIGNED
    if assigned and cfg["batch_vh_family"] != NOSEL:
        res["batch_vh_match"] = dict(vh_t).get(cfg["batch_vh_family"], 0)
    if assigned and cfg["batch_chain"] != NOSEL:
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
def fmt_ambiguity_pair(p):
    return "%s/%s (%s, k%d)" % (p["a"], p["b"], "|".join(p["families"]), p["incompatible"])


def _ambiguity_doc(primers, cfg):
    """모호성 클러스터 설명을 프라이머 FASTA 에서 계산해 만든다. 이름 하드코딩 없음."""
    if not primers or not cfg:
        return ("프라이머 FASTA 를 읽어야 계산됩니다. 이 실행에는 프라이머가 없어 "
                "모호성 쌍을 계산하지 못했습니다.")
    pairs = primer_ambiguity(primers, cfg)
    s = ambiguity_summary(pairs)
    combos = sorted(set("|".join(p["families"]) for p in pairs
                        if not p["same_family"] and p["scored"]
                        and p["tie"] == "certain"))
    return ("프라이머 FASTA 에서 실행 시마다 계산합니다. 같은 group·chain 안에서 판별구간의 "
            "축퇴 집합이 겹쳐 한 read 가 두 프라이머 모두와 허용 미스매치 %d 안에서 맞을 수 "
            "있는 쌍입니다. 프라이머 %d 종에서 certain %d 쌍(같은 family %d · 다른 family %d), "
            "split %d 쌍(같은 family %d · 다른 family %d)이 나옵니다. 등급의 뜻은 '모호성 등급' "
            "항목을, 전체 목록은 05_실행설정의 '모호성' 행을 보세요. 채점 그룹(%s)에서 "
            "certain 등급으로 혼동될 수 있는 family 조합은 %s 입니다."
            % (cfg["primer_max_mismatch"], len(primers),
               s["certain"], sum(1 for p in pairs if p["tie"] == "certain" and p["same_family"]),
               sum(1 for p in pairs if p["tie"] == "certain" and not p["same_family"]),
               s["split"], sum(1 for p in pairs if p["tie"] == "split" and p["same_family"]),
               sum(1 for p in pairs if p["tie"] == "split" and not p["same_family"]),
               ", ".join(SCORED_GROUPS),
               " · ".join(combos) if combos else "없음"))


# 판별 성공분 집계에 쓰는 그룹 대응. call_one 이 채점하는 네 그룹입니다.
CALLED_GROUPS = (("vh", "F1_For", "VH family"),
                 ("jh", "F2_Rev", "JH"),
                 ("vl", "F3_For", "VL V-gene"),
                 ("vj", "F3_Rev", "VL J-gene"))


def _tally(labels):
    c = {}
    for lab in labels:
        c[lab] = c.get(lab, 0) + 1
    # 노트북 Counter.most_common() 과 동일 : 건수 내림차순, 동점은 삽입 순
    return sorted(c.items(), key=lambda kv: -kv[1])


def compose_called(calls, cfg):
    """프라이머 판별에 성공한 클론만 그룹별로 집계한다.

    compose 는 '구조QC 와 판별을 모두 통과한' 클론만 세므로, 구조가 깨졌어도
    판별 구간은 멀쩡한 클론의 정보를 버립니다. 두 질문이 다릅니다.
      compose         : 쓸 수 있는 클론이 얼마나 다양한가
      compose_called  : 프라이머 pool 이 고르게 작동했나
    한 클론이 F1_For 는 성공하고 F3_For 는 실패할 수 있어 그룹마다 모수가
    다릅니다. 그래서 그룹별 n 을 반드시 함께 돌려줍니다.
    """
    groups = {}
    for key, group, label in CALLED_GROUPS:
        ok = [p[key] for p in calls if p[key] is not None and p[key]["ok"]]
        # 동점은 후보가 2 개 이상인 판정입니다. 프라이머 하나가 여러 germline 을
        # 표적해 라벨에 | 가 들어가는 경우(예: JH1|JH2)는 동점이 아닙니다.
        certain, possible, ties = set(), set(), {}
        for c in ok:
            fams = set(c["families"])
            possible |= fams
            if len(c["names"]) > 1:
                lab = "|".join(c["families"]) or "?"
                ties[lab] = ties.get(lab, 0) + 1
            else:
                certain |= fams
        groups[key] = {"group": group, "label": label, "n": len(ok),
                       "tally": _tally(["|".join(c["families"]) or "?" for c in ok]),
                       # 클론 수 분포(tally)는 동점 항목도 그대로 셉니다. 종 수만
                       # 확정/가능으로 나눕니다.
                       "species_certain": len(certain),
                       "species_possible": len(possible),
                       "species_only_possible": sorted(possible - certain),
                       "ambiguous_items": sorted(ties.items(), key=lambda kv: -kv[1])}
        # 화면과 시트가 같은 문구를 쓰도록 표기까지 core 에서 만듭니다.
        groups[key]["species_text"] = fmt_species(groups[key])
    ch = [p["chain"] for p in calls if p["chain"]]
    lens = [p["cdr3_len"] for p in calls if p["cdr3_len"] is not None]
    return {"n_total": len(calls), "groups": groups,
            "chain": {"label": "경쇄", "n": len(ch), "tally": _tally(ch)},
            "cdr3": {"n": len(lens), "lens": sorted(lens),
                     "min": min(lens) if lens else None,
                     "median": _fmt_med(median(lens)),
                     "max": max(lens) if lens else None}}


def batch_assigned(cfg):
    """cfg 에서 이 배치의 지정을 뽑는다. assigned 모드가 아니면 None."""
    if cfg.get("analysis_mode") != MODE_ASSIGNED:
        return None
    fams = [] if cfg["batch_vh_family"] == NOSEL else [cfg["batch_vh_family"]]
    chains = [] if cfg["batch_chain"] == NOSEL else [cfg["batch_chain"]]
    return {"vh_families": fams, "chains": chains}


def _expected_primer(p, group, chain, assigned):
    """이 프라이머가 이 배치에서 나올 수 있는 대상인가.

    그룹마다 배치 지정이 미치는 범위가 다릅니다.
      F1_For   배치 VH family 가 직접 결정한다. families 가 지정 목록과 겹치는
               프라이머만 기대 대상이다.
      F3_For   배치 경쇄가 결정한다. 지정된 chain 의 버킷만 기대 대상이다.
      F3_Rev   F3_For 와 같다.
      F2_Rev   JH 는 Rev-2 pool 을 통째로 쓰므로 배치 지정과 무관하다. 전부 기대 대상.
      F1_Rev / F2_For  채점 대상이 아니라 커버리지 계산에서 아예 제외한다.
    지정이 비어 있는 축(예: batch_chain 이 "(지정 안 함)")은 제한하지 않습니다.
    """
    if assigned is None:
        return True
    if group == "F1_For":
        fams = assigned.get("vh_families") or []
        if not fams:
            return True
        # family 를 못 읽은 프라이머는 배제할 근거가 없으므로 기대 대상으로 둡니다.
        return (not p["families"]) or bool(set(p["families"]) & set(fams))
    if group in ("F3_For", "F3_Rev"):
        chains = assigned.get("chains") or []
        return (not chains) or (chain in chains)
    return True


def primer_coverage(calls, primers, cfg, assigned=None):
    """그룹·사슬별로 한 번이라도 판정에 나온 프라이머와 미관측 프라이머.

    동점으로 여러 후보가 나온 경우 그 후보를 전부 관측으로 셉니다. 어느
    프라이머였는지 모르므로 어느 쪽도 dropout 이라 말할 수 없기 때문입니다.
    첫 번째만 세면 실제로 관측된 프라이머가 dropout 으로 잘못 보고됩니다.

    표본 부족 판정은 균등 사용을 가정한 하한 기준입니다. 종 수 S 인 pool 에서
    n 개를 뽑아 특정 한 종이 한 번도 안 나올 확률 (1 - 1/S)^n 이
    RULES["COVERAGE_ALPHA"] 이상이면 미관측을 dropout 이라 말할 수 없습니다.
    실제 germline 사용 빈도는 균등하지 않으므로 실제 확률은 이보다 큽니다.
    """
    seen, n_by = {}, {}
    for key, group, _lab in CALLED_GROUPS:
        for p in calls:
            c = p[key]
            if c is None or not c["ok"]:
                continue
            b = (group, c["chain"])
            n_by[b] = n_by.get(b, 0) + 1
            seen.setdefault(b, set()).update(c["names"])
    buckets = {}
    for p in primers:
        if p["group"] in SCORED_GROUPS:
            buckets.setdefault((p["group"], p["chain"]), []).append(p)
    out = []
    def brief(plist):
        return [{"name": x["name"], "families": list(x["families"])} for x in plist]

    for (group, chain), plist in sorted(buckets.items()):
        names = seen.get((group, chain), set())
        n, S = n_by.get((group, chain), 0), len(plist)
        exp_names = set(p["name"] for p in plist
                        if _expected_primer(p, group, chain, assigned))
        exp = [p for p in plist if p["name"] in exp_names]
        miss = [p for p in plist if p["name"] not in names]
        miss_exp = [p for p in miss if p["name"] in exp_names]
        miss_oth = [p for p in miss if p["name"] not in exp_names]
        E = len(exp)
        # 검정력 모수는 전체가 아니라 기대 대상 수입니다. 지정하지 않은 프라이머를
        # 모수에 넣으면 종 수가 부풀려져 검정력이 과소평가됩니다.
        power = ((1.0 - 1.0 / E) ** n) if E else None
        out.append({"group": group, "chain": chain, "total": S, "n": n,
                    "observed": S - len(miss),
                    "unobserved": brief(miss),
                    # 아래가 배치 지정을 반영한 값입니다.
                    "total_n": S, "expected_n": E,
                    "observed_expected": E - len(miss_exp),
                    "missing_expected": brief(miss_exp),
                    "missing_other": brief(miss_oth),
                    "power_p": power,
                    # p_miss 는 예전 이름이며 power_p 와 같은 값입니다.
                    "p_miss": power,
                    "underpowered": bool(E) and power >= RULES["COVERAGE_ALPHA"]})
    return out


# --- 용어설명 예시를 프라이머 FASTA 에서 고른다 -------------------------------
# 개념 설명이되 예시가 있어야 읽히는 항목들입니다. 이름을 코드에 박으면 프라이머
# 세트가 바뀔 때 설명이 사실과 달라지므로, 조건을 정해 실행 시 고릅니다.
# 조건에 맞는 것이 없으면 예시 없이 개념만 서술합니다.
def _pick_delta_pair(primers, cfg):
    """Δ 예시 : F2_Rev 에서 비호환 위치가 가장 적은 쌍 (이름순 tie-break)."""
    cand = [p for p in primer_ambiguity(primers, cfg) if p["group"] == "F2_Rev"]
    if not cand:
        return None
    return sorted(cand, key=lambda p: (p["incompatible"], p["a"], p["b"]))[0]


def _pick_same_family_pair(primers, cfg):
    """모호성 표기 예시 : 채점 그룹의 certain 등급이면서 같은 family 인 쌍."""
    cand = [p for p in primer_ambiguity(primers, cfg)
            if p["tie"] == "certain" and p["same_family"] and p["scored"]]
    if not cand:
        return None
    return sorted(cand, key=lambda p: (p["incompatible"], p["group"], p["a"], p["b"]))[0]


def _pick_superset_pair(primers):
    """CDR3 경계 예시 : F1_Rev / F2_For 에서 축퇴 공간이 다른 family 프라이머를
    완전히 포함하는 쌍. (a 의 모든 위치가 b 를 집합으로 포함)"""
    best = None
    for group in ("F1_Rev", "F2_For"):
        pl = [p for p in primers if p["group"] == group]
        for a in pl:
            for b in pl:
                if a["name"] == b["name"]:
                    continue
                if set(a["families"]) == set(b["families"]):
                    continue
                ca, cb = a["core"], b["core"]
                n = min(len(ca), len(cb))
                if n and all(set(IUPAC.get(cb[i], cb[i])) <= set(IUPAC.get(ca[i], ca[i]))
                             for i in range(n)):
                    key = (group, a["name"], b["name"])
                    if best is None or key < best[0]:
                        best = (key, a, b, n)
    return best


def _pick_multi_target(primers):
    """프라이머 하나가 여러 germline 을 표적하는 예시 (후보 1 개인데 라벨에 | 가 붙는 경우)."""
    cand = [p for p in primers
            if p["group"] in SCORED_GROUPS and len(p["families"]) > 1]
    return sorted(cand, key=lambda p: p["name"])[0] if cand else None


def _doc_pipe_meanings(primers, cfg):
    """family 칸의 | 가 갖는 두 가지 뜻. 예시는 프라이머 세트에서 고릅니다."""
    multi = _pick_multi_target(primers) if primers else None
    tie = None
    if primers and cfg:
        cand = [p for p in primer_ambiguity(primers, cfg)
                if not p["same_family"] and p["scored"]]
        tie = sorted(cand, key=lambda p: (p["incompatible"], p["group"],
                                          p["a"], p["b"]))[0] if cand else None
    a = ("예: %s 는 %s 한 프라이머가 여러 germline 을 함께 잡도록 설계된 것입니다. "
         % ("|".join(multi["families"]), multi["name"])) if multi else ""
    b = ("예: %s 는 %s 와 %s 가 같은 미스매치로 걸린 것입니다. "
         % ("|".join(tie["families"]), tie["a"], tie["b"])) if tie else ""
    return ("family 칸에 여러 값이 | 로 적히는 경우는 두 가지이고 뜻이 다릅니다. "
            "(a) 프라이머 하나가 여러 germline 을 표적하는 경우 — " + a +
            "(b) 프라이머 여러 개가 서로 구분되지 않아 동점인 경우 — " + b +
            "03_프라이머판별의 '후보수' 열로 구분합니다. 1 이면 (a), 2 이상이면 (b) 입니다. "
            "(b) 라면 같은 시트의 '모호성' 열에 어느 쌍이 어느 등급으로 관여했는지 나옵니다.")


def _doc_delta(primers, cfg):
    head = "최상위 후보와 그 다음 후보의 미스매치 차이. 클수록 판정 근거가 두껍습니다. "
    p = _pick_delta_pair(primers, cfg) if primers and cfg else None
    if p:
        head += ("이번 프라이머 세트에서는 %s 과 %s 가 판별구간 %d nt 중 %d 곳에서만 "
                 "갈리므로 그 그룹의 Δ 가 작습니다. " % (p["a"], p["b"], p["cmp_len"],
                                                    p["incompatible"]))
    return (head + "해당 구간 Q 가 50 이상이면 오독 확률이 10^-5 수준이라 1 nt 차이도 "
            "신뢰할 만합니다. 판정이 실패(X)한 경우의 Δ 는 실패한 후보들끼리의 순위차일 "
            "뿐이므로 읽지 않습니다.")


def _doc_ambig_notation(primers, cfg):
    p = _pick_same_family_pair(primers, cfg) if primers and cfg else None
    if p:
        shown = fmt_call({"ok": True, "families": p["families"], "names": [p["a"], p["b"]]})
        head = "%s 처럼 표기합니다. " % shown
    else:
        head = "family (프라이머 후보들) 형태로 표기합니다. "
    return (head + "왼쪽 family 는 확정이고, 괄호 안은 축퇴 공간이 겹쳐 서로 구분할 수 "
            "없는 프라이머 후보 집합입니다. 하나로 골라 적는 것은 근거 없는 정보이므로 "
            "하지 않습니다. family 칸에 여러 germline 이 | 로 적히는 경우도 같은 기호를 "
            "쓰지만, 프라이머 하나가 여러 germline 을 표적하는 것이라면 후보 프라이머는 "
            "1 개입니다.")


def _doc_cdr3_boundary(primers):
    hit = _pick_superset_pair(primers) if primers else None
    if hit:
        (_g, a_name, b_name), a, b, n = hit
        why = ("%s(%s) 의 축퇴 공간이 %s(%s) 를 판별구간 %d nt 전체에서 완전히 포함해 "
               "유일 배정이 불가능합니다" % (a_name, "|".join(a["families"]) or "?",
                                        b_name, "|".join(b["families"]) or "?", n))
    else:
        why = ("일부 프라이머의 축퇴 공간이 다른 family 프라이머를 완전히 포함해 "
               "유일 배정이 불가능합니다")
    return ("분류에 쓰지 않고 정합성 검사로만 씁니다. 이유 두 가지 - (1) " + why +
            ". (2) overlap extension 중 proofreading 중합효소의 3'->5' exonuclease 가 "
            "미스매치된 3' 말단을 제거하고 상대 fragment 를 주형으로 재연장하기 때문에, "
            "이 구간은 프라이머 서열과 원 주형의 혼합이 됩니다. 그래서 양끝 %d nt 를 "
            "제외하고 채점합니다." % RULES["EXO_TRIM"])


def fmt_species(g):
    """확정 종 수 · 가능 종 수를 사람이 읽을 한 줄로."""
    cert, poss = g["species_certain"], g["species_possible"]
    ties = g["ambiguous_items"]
    if not ties:
        return "%d 종 (동점 없음)" % cert
    tie_txt = " · ".join("%s %d 클론" % (lab, n) for lab, n in ties)
    if cert == poss:
        note = "새 family 없음"
    else:
        note = "%s 는 동점으로만 관측" % ", ".join(g["species_only_possible"])
    return "확정 %d 종 · 가능 %d 종 (동점 %s — %s)" % (cert, poss, tie_txt, note)


def _pct(v, n):
    """모수를 반드시 붙인 비율 표기. 모수 없는 비율은 오해를 부릅니다."""
    return "%d (%.1f%%, n=%d)" % (v, 100.0 * v / n, n) if n else "%d (n=0)" % v


def glossary(primers=None, cfg=None):
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
        ["QC 랜드마크", "이상 랜드마크 최저Q",
         "상태가 OK 나 NA 가 아닌 랜드마크들의 최저 Phred 값이다. 이 값이 낮으면 "
         "그 갭이나 치환이 실제 결실이 아니라 판독 오류일 수 있다. QC4 는 read 3' 끝에 "
         "있어 품질이 먼저 무너지고, QC2 는 read 중간이라 상대적으로 신뢰할 만하다. "
         "판정에는 반영되지 않으므로 FAIL 이라도 이 값이 낮으면 크로마토그램을 확인한다. "
         "이상이 없는 클론에도 낮은 Q 는 흔하므로, 낮은 Q 자체가 아니라 '이상이 난 자리의 "
         "Q 가 낮다'는 것이 신호다. 한 클론에 이상이 둘 이상이면 이 열은 그중 최저만 "
         "보여주므로 개별 값은 02_구조QC상세에서 확인한다. 저품질 표시 기준은 "
         "RULES 의 LM_LOWQ_MARK(%d)이며 표시 전용이다." % RULES["LM_LOWQ_MARK"]],

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
        ["프라이머 판별", "Δ (차순위 간격)", _doc_delta(primers, cfg)],
        ["프라이머 판별", "모호성 표기", _doc_ambig_notation(primers, cfg)],
        ["프라이머 판별", "모호성 클러스터", _ambiguity_doc(primers, cfg)],
        ["프라이머 판별", "family 칸의 | 가 갖는 두 가지 뜻",
         _doc_pipe_meanings(primers, cfg)],        ["프라이머 판별", "AMBIG_CALL?",
         "판정 후보 프라이머가 2 개 이상이고 그 후보들이 서로 다른 family 를 가리키는 "
         "경우다. 결과의 family 칸에 VH4|VH6 처럼 여러 값이 적히며, 어느 쪽인지 이 판별만 "
         "으로는 좁힐 수 없다. 배치 지정 여부와 무관하게 붙는다 — 배치 지정 family 가 "
         "후보 집합 안에 있으면 WRONG_FAMILY 도 AMBIG_FAMILY? 도 붙지 않아 family 가 "
         "확정되지 않았다는 사실이 드러나지 않던 사각지대를 메운다. 프라이머 하나가 여러 "
         "germline 을 표적해 | 가 붙는 경우(예: JH1|JH2)는 후보가 1 개이므로 해당하지 "
         "않는다. 어느 쌍이 겹쳤는지는 03_프라이머판별의 모호성 열에 있다."],
        ["프라이머 판별", "AMBIG_FAMILY?",
         "배치 지정 VH family 와 판정 family 가 다르지만, 두 family 가 알려진 모호성 "
         "쌍으로 이어져 있어 교차오염으로 단정할 수 없는 경우입니다. 모호성 쌍이 없으면 "
         "WRONG_FAMILY 가 붙습니다. 두 플래그를 가르는 것은 프라이머 설계상 구분 "
         "가능한가이지 증거의 세기가 아니므로, AMBIG_FAMILY? 라도 교차오염이 아니라는 "
         "뜻은 아닙니다. 03_프라이머판별의 모호성 열과 서열을 함께 보세요."],
        ["프라이머 판별", "모호성 등급 certain / split",
         "비호환 위치 수 k 와 허용 미스매치 T 로 나눕니다. certain (k <= T) 은 한쪽 "
         "프라이머와 완전히 같은 read 가 다른 쪽과도 허용 범위 안에 들어, 동점이 "
         "일상적으로 발생합니다. split (T < k <= 2T) 은 read 가 양쪽 프라이머 모두에서 "
         "벗어날 때만 동점이 되므로 가능하지만 조건부입니다. k > 2T 면 어떤 read 도 두 "
         "프라이머 모두와 허용 범위 안에서 맞을 수 없어 동점이 불가능합니다."],
        ["프라이머 판별", "모호성 · family 가 다른 쌍",
         "같은 family 안의 동점은 프라이머 이름만 모호하고 판정되는 family 는 확정됩니다. "
         "family 가 다른 동점은 판정 자체가 갈립니다 — 결과에 VH4|VH6 처럼 여러 family 가 "
         "함께 적히고, 어느 쪽인지는 이 판별만으로 결정할 수 없습니다. 이런 쌍이 있는지, "
         "어느 조합인지는 프라이머 FASTA 마다 다르므로 실행 시 계산해 05_실행설정에 "
         "기록합니다."],
        ["프라이머 판별", "그룹 F1_For", "VH FR1. NotI 를 앵커로 위치를 잡습니다. 판별 신뢰도 높음."],
        ["프라이머 판별", "그룹 F2_Rev", "JH + 링커. 링커 시작을 앵커로 잡는 역방향 프라이머. 신뢰도 높음."],
        ["프라이머 판별", "그룹 F3_For",
         "VL FR1. 링커 끝을 앵커로 잡습니다. kappa / lambda 를 각각 채점해 낮은 쪽을 택합니다."],
        ["프라이머 판별", "그룹 F3_Rev", "VL J. AscI 를 앵커로 잡는 역방향 프라이머. 신뢰도 높음."],
        ["프라이머 판별", "CDR3 경계 (F1_Rev / F2_For)", _doc_cdr3_boundary(primers)],
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
        ["판정 코드", "PARENTAL",
         "모클론 스터퍼 서열 검출. 미절단 또는 단일절단 벡터가 재결합한 배경 클론. "
         "길이만 일치하고 서열이 없으면 PARENTAL? 로 따로 표시한다."],
        ["판정 코드", "PARENTAL?",
         "인서트 길이가 스터퍼 보유 클론과 정확히 같지만(%d bp) 스터퍼 서열은 "
         "검출되지 않은 경우. 확정이 아니라 확인 필요 등급이다. read 앞쪽 품질이 "
         "나빠 스터퍼 구간이 트리밍됐을 수도 있고, 우연히 같은 길이인 다른 산물일 "
         "수도 있다. 크로마토그램으로 확인한다." % CONST["STUFFER_INSERT_BP"]],
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

        ["대조군 판정", "음성 대조군 모드",
         "벡터만 ligation 한 대조군을 읽는 모드다. 기대값이 일반 배치와 정반대여서 "
         "판정 어휘가 다르다. 일반 배치에서 PASS 는 성공이지만 대조군에서 완전한 "
         "scFv 가 나오면 오염이다. 구조 QC(qc_one)는 동일하게 돌리고 그 결과를 "
         "대조군 관점으로 다시 읽는다. 아래 순서로 먼저 맞는 것을 판정으로 쓴다."],
        ["대조군 판정", "CONCATEMER",
         "NotI / AscI / 링커가 2 회 이상 검출. 벡터가 여러 번 이어붙은 구조 이상이라 "
         "다른 어떤 해석보다 먼저 잡는다."],
        ["대조군 판정", "MIXED",
         "트레이스 2 순위 피크 비율이 임계를 넘음. 콜로니 하나를 집지 못한 것이므로 "
         "인서트 해석 자체를 신뢰할 수 없다."],
        ["대조군 판정", "EMPTY_VECTOR",
         "NotI 또는 AscI 미검출. 클로닝 자리가 파괴된 채 재결합한 벡터다. "
         "음성 대조군에서 기대되는 배경 유형이다."],
        ["대조군 판정", "PARENTAL",
         "스터퍼 서열 검출. 미절단 또는 단일절단 벡터가 그대로 살아남은 것이다. "
         "벡터 준비물을 새로 만들면 주된 유형이 될 수 있다."],
        ["대조군 판정", "PARENTAL?",
         "인서트 길이가 스터퍼 보유 클론과 같으나(%d bp) 스터퍼 서열은 미검출. "
         "read 앞쪽 품질 저하로 놓쳤을 수 있어 확정하지 않는다."
         % CONST["STUFFER_INSERT_BP"]],
        ["대조군 판정", "CONTAMINATED",
         "링커 검출 + 프레임 정상 + 인서트가 설정 범위 내. 즉 완전한 scFv 가 "
         "들어갔다. 음성 대조군에는 인서트를 넣지 않았으므로 ligation 단계에서 "
         "인서트가 섞여 들어갔다는 뜻이다. 대조군 판정 중 가장 심각하며, 같은 "
         "반응에서 나온 본 실험 결과 전체의 신뢰도에 영향을 준다."],
        ["대조군 판정", "CONTAMINATED?",
         "링커 검출 + 인서트가 설정 범위 내인데 프레임만 밀린 경우. 물음표는 오염 "
         "여부가 아니라 프레임 이상을 가리킨다. 대조군에서 이것은 '오염이 아니다' 가 "
         "아니라 '오염된 인서트에 프레임 문제도 있다' 는 뜻이므로 CONTAMINATED 와 "
         "같은 등급으로 다룬다. 프레임 정상만 오염으로 보면 오염 신호를 놓친다."],
        ["대조군 판정", "PARTIAL_INSERT",
         "링커는 검출되는데 인서트가 길이 하한 미만인 경우. VL 쪽이 잘린 부분 조립 "
         "산물이다. scFv 유래 물질이 들어간 것은 맞지만 완전한 scFv 는 아니므로 "
         "CONTAMINATED 계열과 구분한다. 길이는 06_서열에서 직접 확인한다."],
        ["대조군 판정", "CARRYOVER",
         "인서트는 있으나 링커 미검출. scFv 가 아닌 무언가가 들어간 것이다. "
         "CONTAMINATED 가 ligation 단계에서 인서트가 섞인 것이라면, 이쪽은 벡터 "
         "준비물이나 공용 시약에 이미 다른 DNA 가 있었을 가능성을 가리킨다. "
         "들어간 것이 무엇인지는 06_서열의 인서트 염기서열로 직접 확인해야 한다."],
        ["대조군 판정", "CHECK",
         "위 어디에도 맞지 않음. 억지로 분류하지 않고 남긴다. 이 값이 나오면 "
         "분류 규칙이 실제 벡터 거동을 못 따라가고 있다는 신호이므로 서열을 "
         "직접 보고 유형을 추가할지 판단한다."],
        ["대조군 판정", "인서트 md5",
         "인서트 염기서열의 지문. 음성 대조군에는 인서트를 넣지 않았으므로, 같은 "
         "지문이 본 실험 결과에도 있으면 같은 분자가 양쪽에 있다는 뜻이다. 모드별로 "
         "파일이 분리되어 자동 대조는 하지 않으니 사람이 대조한다."],

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

        ["배치 조성", "기대 대상과 배치 무관",
         "커버리지는 배치에 지정된 프라이머만 기대 대상으로 센다. 지정하지 않은 family 나 "
         "경쇄의 프라이머는 나올 수 없으므로 미관측이어도 dropout 이 아니다. 검정력 계산의 "
         "모수도 기대 대상 수를 쓴다 — 전체 수를 쓰면 검정력이 과소평가된다. JH 는 pool 을 "
         "통째로 쓰므로 배치 지정과 무관하게 전부 기대 대상이다. Library 와 Negative control "
         "모드에는 지정이 없으므로 전체가 기대 대상이 된다."],
        ["배치 조성", "확정 종 수와 가능 종 수",
         "동점 판정(VH4|VH6 등)은 어느 family 인지 확정되지 않았으므로 종 수를 두 값으로 "
         "센다. 확정 종 수는 단독 판정만으로 세고, 가능 종 수는 동점 항목의 구성원까지 "
         "합집합으로 센다. 두 값이 같으면 동점이 새 family 를 만들지 않는다는 뜻이고, "
         "다르면 그 차이만큼 다양성이 불확실하다. 클론 수 분포에서는 동점 항목도 "
         "그대로 세므로 두 수치의 합계가 다를 수 있다. 동점은 후보 프라이머가 2 개 "
         "이상인 판정을 말하며, 프라이머 하나가 여러 germline 을 표적해 라벨에 | 가 "
         "들어가는 경우(예: JH1|JH2)는 동점이 아니라 단독 판정으로 센다."],

        ["재현성", "param_hash",
         "판정 임계값 16 개를 정렬해 만든 지문. 두 배치의 결과를 비교하려면 이 값이 같아야 합니다. "
         "실험 설계와 cDNA / RNA 출처는 해시에 포함되지 않고 05_실행설정 시트에 개별 기록됩니다. "
         "설계 파라미터는 design_hash 로 따로 관리합니다."],
        ["재현성", "design_hash",
         "판정에 직접 관여하는 실험 설계 %d 개의 지문. batch_vh_family 와 batch_chain 은 "
         "WRONG_FAMILY / WRONG_CHAIN 플래그를, f1_for_mode 와 f2_for_mode 는 CDR3 경계 "
         "정합성 검사의 수행 여부를 결정한다. 임계값은 배치 간에 같아야 하고 설계는 "
         "배치마다 달라야 정상이므로 param_hash 와 분리해 둔다. 두 배치를 비교할 때 "
         "param_hash 는 같아야 하고, design_hash 가 다르면 무엇이 다른지 05_실행설정에서 "
         "확인해야 한다." % len(JUDGMENT_DESIGN_KEYS)],
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
    "CDR3-H3", "CDR3-H3(aa)", "CDR3경계", "param_hash",
    "이상 랜드마크 최저Q", "최저Q 랜드마크"]

_STOP_TXT = {True: "OK", False: "FAIL", None: "-"}


def _yn(v, yes="YES", no="-"):
    return yes if v else no


def final_verdict(r, p):
    if r["verdict"] not in ("PASS", "PASS*"):
        return "FAIL"
    if p["flags"]:
        return "FAIL"
    return "PASS" if r["verdict"] == "PASS" else "PASS*"


def fmt_call_ambiguity(c):
    """한 판정에 관여한 모호성 쌍을 03_프라이머판별 한 칸으로."""
    if not c:
        return "-"
    out = []
    for x in c.get("ambiguity") or []:
        out.append("동점 %s/%s (%s, %s, k%d)"
                   % (x["a"], x["b"], x["families"], x["tie"], x["k"]))
    for x in c.get("runner_up_pairs") or []:
        out.append("차순위 %s/%s (%s, %s, k%d)"
                   % (x["a"], x["b"], x["families"], x["tie"], x["k"]))
    return " · ".join(out) or "-"


def _fam(call):
    return "|".join(call["families"]) if (call and call["ok"]) else "-"


_LM_KEYS = ("QC1", "QC2", "QC3", "QC4")


def landmark_quality(r, cfg=None):
    """이미 계산된 랜드마크 결과를 읽어 품질 근거로 요약한다.

    qc_one / check_landmark 는 건드리지 않습니다. 판정에도 쓰이지 않고
    사람이 갭의 진위를 판단할 근거를 모으기만 합니다.

    worst_minq 는 '이상이 난 자리' 의 최저 Q 입니다. 이상이 없으면 None 이며,
    이상 여부와 무관한 전체 최저 Q(all_minq)와 혼동하면 안 됩니다. 이상이 없는
    클론에도 낮은 Q 는 흔하므로 낮은 Q 자체는 신호가 아닙니다.
    한 클론에 이상이 둘 이상이면 worst_minq 는 그중 최저 하나만 가리키므로
    개별 값은 02_구조QC상세에서 봐야 합니다.
    cfg 는 지금 쓰지 않지만 임계값이 cfg 로 옮겨갈 때를 위해 받아 둡니다.
    """
    anomalies = []
    for k in _LM_KEYS:
        q = r["qc"][k]
        if q["status"] in ("OK", "NA"):
            continue
        anomalies.append({"lm": k, "status": q["status"], "sub": q["sub"],
                          "gap": q["gap"], "minq": q["minq"]})
    scored = [a for a in anomalies if a["minq"] is not None]
    worst = min(scored, key=lambda a: a["minq"]) if scored else None
    every = [r["qc"][k]["minq"] for k in _LM_KEYS if r["qc"][k]["minq"] is not None]
    return {"anomalies": anomalies,
            "worst_minq": worst["minq"] if worst else None,
            "worst_lm": worst["lm"] if worst else None,
            "all_minq": min(every) if every else None}


def build_summary(qc_results, calls, cfg, meta=None):
    """batch / date 는 파일명이 아니라 폼에서 받은 meta 값을 씁니다."""
    meta = meta or {}
    batch = meta.get("batch_label", "")
    date = meta.get("batch_date", "")
    rows = []
    by_id = {p["id"]: p for p in calls}
    for r in qc_results:
        p = by_id[r["id"]]
        ins = r["insert_bp"]
        lq = landmark_quality(r, cfg)
        rows.append([
            r["id"], batch, date, r["primer"], r["direction"],
            final_verdict(r, p), r["verdict"], ", ".join(p["flags"]) or "-",
            r["used_bp"], ins, (ins % 3) if ins is not None else None, r["d1"], r["d2"],
            r["qc"]["QC1"]["status"], r["qc"]["QC2"]["status"],
            r["qc"]["QC3"]["status"], r["qc"]["QC4"]["status"],
            _STOP_TXT[r["stop_ok"]], _yn(r["stuffer"]), _yn(bool(r["repeat"])),
            round(r["mix_pct"], 2) if r["mix_pct"] is not None else None,
            _fam(p["vh"]), _fam(p["jh"]), p["chain"] or "-",
            _fam(p["vl"]), _fam(p["vj"]),
            p["cdr3"] or "-", p["cdr3_len"], p["boundary"].get("status", "-"),
            cfg["param_hash"], lq["worst_minq"], lq["worst_lm"] or ""])
    return rows


def build_sheets(qc_results, calls, cfg, comp, meta, primers=None, coverage=None):
    """xlsx 시트 데이터를 순서대로 담은 리스트를 반환한다."""
    if cfg.get("analysis_mode") == MODE_NEGCTRL:
        return build_negctrl_sheets(qc_results, cfg, meta, primers)
    by_id = {p["id"]: p for p in calls}
    sheets = []

    # 01 판정요약
    sheets.append({
        "title": "01_판정요약", "headers": SUMMARY_HEADERS,
        "rows": build_summary(qc_results, calls, cfg, meta), "wrap": [7], "maxw": 60,
        "note": "최종판정 = 구조QC 통과 AND 프라이머 판별 플래그 없음. "
                "각 컬럼의 정의는 06_용어설명 시트 참조."})

    sheets.append(_sheet_struct_qc(qc_results, "02_구조QC상세"))
    # 03 프라이머판별
    h3 = ["clone", "그룹", "대상", "판정 family", "후보 프라이머", "후보수",
          "미스매치", "Δ(차순위간격)", "판별구간(nt)", "불일치위치(1-based)",
          "read상 시작위치", "허용 미스매치", "통과", "모호성"]
    grp = [("vh", "F1_For", "VH family"), ("jh", "F2_Rev", "JH"),
           ("vl", "F3_For", "VL V-gene"), ("vj", "F3_Rev", "VL J-gene")]
    r3 = []
    for p in calls:
        for key, gname, target in grp:
            c = p[key]
            if c is None:
                r3.append([p["id"], gname, target, "-", "-", 0, None, None, None,
                           "-", None, cfg["primer_max_mismatch"], "앵커없음", "-"])
                continue
            r3.append([p["id"], gname, target,
                       "|".join(c["families"]) or "?", "|".join(c["names"]),
                       len(c["names"]), c["mm"], c["delta"], c["core_len"],
                       ", ".join(str(x + 1) for x in c["mmpos"]) or "-",
                       c["start"] + 1, cfg["primer_max_mismatch"],
                       "OK" if c["ok"] else "FAIL", fmt_call_ambiguity(c)])
        b = p["boundary"]
        r3.append([p["id"], "F2_For", "CDR3 경계 정합성", "-", b.get("primer", "-"), 0,
                   b.get("mm"), None, b.get("cmp"), "-", None,
                   cfg["primer_max_mismatch"], b.get("status", "-"), "-"])
    sheets.append({
        "title": "03_프라이머판별", "headers": h3, "rows": r3, "wrap": [13], "maxw": 60,
        "note": "그룹별 1행. Δ 는 최상위와 차순위 후보의 미스매치 간격이며 "
                "클수록 판정 근거가 두껍습니다."})

    # 04 배치조성
    h4 = ["구분", "항목", "값", "비고"]
    r4 = [["개요", "배치", ", ".join(cfg["batches"]) or "-", ""],
          ["개요", "분석 클론 수", comp["n_total"], ""],
          ["개요", "구조QC + 프라이머판별 모두 통과", comp["n_good"],
           "'분포(통과)' 는 이 클론들만 집계"]]
    called = compose_called(calls, cfg)
    for key, group, label in CALLED_GROUPS:
        g = called["groups"][key]
        r4.append(["개요", "%s 판별 성공 (%s)" % (label, group),
                   _pct(g["n"], called["n_total"]),
                   "구조 결함과 무관하게 이 그룹의 판별에 성공한 클론 수. "
                   "'분포(판별)' 의 모수"])
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
                r4.append(["분포(통과)", lab, "%s : %s" % (k, _pct(v, comp["n_good"])), memo])
            r4.append(["분포(통과)", lab + " (종 수)", len(items),
                       "모수 n=%d" % comp["n_good"]])
        else:
            r4.append(["분포(통과)", lab, "-", memo + " (n=%d)" % comp["n_good"]])
    # 분포(판별) : 구조 결함과 무관하게 판별에 성공한 클론만. 그룹마다 모수가 다릅니다.
    for key, group, label in CALLED_GROUPS:
        g = called["groups"][key]
        if g["tally"]:
            for k, v in g["tally"]:
                r4.append(["분포(판별)", label, "%s : %s" % (k, _pct(v, g["n"])), group])
            r4.append(["분포(판별)", label + " (종 수)", fmt_species(g),
                       "클론 수 분포는 동점 항목도 그대로 셉니다. 모수 n=%d" % g["n"]])
        else:
            r4.append(["분포(판별)", label, "-", "%s · 판별 성공 0 건" % group])
    cg = called["chain"]
    for k, v in cg["tally"]:
        r4.append(["분포(판별)", cg["label"], "%s : %s" % (k, _pct(v, cg["n"])),
                   "kappa/lambda 동점이 아닌 클론만"])
    if not cg["tally"]:
        r4.append(["분포(판별)", cg["label"], "-", "경쇄가 확정된 클론 없음"])

    # 커버리지 : 미관측 프라이머가 dropout 신호인지, 표본이 부족한 것인지
    alpha = RULES["COVERAGE_ALPHA"]
    names_of = lambda lst: ", ".join("%s(%s)" % (u["name"], "|".join(u["families"]) or "?")
                                     for u in lst)
    for cv in coverage or []:
        tag = "%s / %s" % (cv["group"], cv["chain"])
        if not cv["expected_n"]:
            # 배치에 지정되지 않은 버킷. 표본 부족과는 다른 상태입니다.
            r4.append(["커버리지", tag, "해당 배치 없음",
                       "이 배치에 지정되지 않아 나올 수 없는 프라이머 %d 종입니다. "
                       "미관측이 dropout 이 아니므로 검정력을 계산하지 않습니다."
                       % cv["total_n"]])
            continue
        r4.append(["커버리지", tag,
                   "%d / %d 관측 (n=%d)"
                   % (cv["observed_expected"], cv["expected_n"], cv["n"]),
                   "기대 대상 %d 종 (전체 %d) · 미관측 %d 종 · 균등 사용 가정에서 "
                   "특정 1 종이 안 나올 확률 (1-1/%d)^%d = %.3f %s %.2f%s. 실제 "
                   "germline 사용 빈도는 균등하지 않으므로 이 값은 하한입니다."
                   % (cv["expected_n"], cv["total_n"], len(cv["missing_expected"]),
                      cv["expected_n"], cv["n"], cv["power_p"],
                      ">=" if cv["underpowered"] else "<", alpha,
                      " · 표본 부족" if cv["underpowered"] else " · 검정력 충분")])
        if cv["missing_expected"]:
            r4.append(["커버리지", tag + " 미관측",
                       names_of(cv["missing_expected"]),
                       "기대 대상인데 한 번도 나오지 않았습니다. dropout 후보입니다."])
        if cv["missing_other"]:
            r4.append(["커버리지", tag + " 배치 무관",
                       names_of(cv["missing_other"]),
                       "배치 무관 항목은 이 배치에 지정되지 않아 나올 수 없었던 "
                       "것이며 dropout 이 아닙니다."])

    if comp["cdr3_lens"]:
        r4.append(["CDR3-H3", "길이 목록 (aa)",
                   ", ".join(str(x) for x in comp["cdr3_lens"]), ""])
        r4.append(["CDR3-H3", "중앙값 (aa)", comp["cdr3_median"], ""])
        r4.append(["CDR3-H3", "길이 분포 (통과)",
                   "최소 %d · 중앙 %s · 최대 %d"
                   % (min(comp["cdr3_lens"]), comp["cdr3_median"],
                      max(comp["cdr3_lens"])),
                   "모수 n=%d" % len(comp["cdr3_lens"])])
    cd = called["cdr3"]
    if cd["n"]:
        r4.append(["CDR3-H3", "길이 분포 (판별)",
                   "최소 %d · 중앙 %s · 최대 %d" % (cd["min"], cd["median"], cd["max"]),
                   "CDR3-H3 가 추출된 클론 전체. 모수 n=%d" % cd["n"]])
    r4.append(["CDR3-H3", "중복 서열 수", len(comp["cdr3_dup"]),
               ", ".join(comp["cdr3_dup"]) if comp["cdr3_dup"]
               else "동일 CDR3-H3 를 가진 클론 없음"])
    # 랜드마크 이상의 판독 품질. 판정에는 반영되지 않는 표시 전용 근거입니다.
    mark = RULES["LM_LOWQ_MARK"]
    lqs = [landmark_quality(r, cfg) for r in qc_results]
    with_anom = [q for q in lqs if q["anomalies"]]
    lowq = [q for q in with_anom
            if q["worst_minq"] is not None and q["worst_minq"] < mark]
    r4.append(["점검", "랜드마크 이상 중 저품질(Q<%d)" % mark,
               "%d / %d 건" % (len(lowq), len(with_anom)),
               "판정에는 반영되지 않습니다. 해당 갭이 실제 결실인지 판독 오류인지는 "
               "크로마토그램으로 확인하세요. 클론별 값은 01_판정요약의 "
               "'이상 랜드마크 최저Q' 열, 랜드마크별 값은 02_구조QC상세에 있습니다."])
    for k in _LM_KEYS:
        hit = [a for q in lqs for a in q["anomalies"] if a["lm"] == k]
        if not hit:
            continue
        low = [a for a in hit if a["minq"] is not None and a["minq"] < mark]
        qs = sorted(a["minq"] for a in hit if a["minq"] is not None)
        med = _fmt_med(median(qs)) if qs else "-"
        r4.append(["점검", "%s 이상" % k,
                   "%d 건 · 저품질 %d 건" % (len(hit), len(low)),
                   "최저 Q 중앙값 %s" % med])
    if comp["overlong_suspect"]:
        r4.append(["점검", "F1_For 불일치 앞쪽 편중",
                   ", ".join(comp["overlong_suspect"]),
                   "판별구간 앞 %d nt 안에만 불일치가 몰린 클론. "
                   "For-Over-long 프라이머 3' 말단 설계 영향 가능성"
                   % comp["overlong_zone"]])
    sheets.append({
        "title": "04_배치조성", "headers": h4, "rows": r4, "wrap": [3], "maxw": 60,
        "note": "구조QC 와 프라이머 판별을 모두 통과한 클론만 집계합니다."})

    sheets.append(_sheet_config(cfg, meta, "05_실행설정", primers))
    sheets.append(_sheet_glossary("06_용어설명", primers, cfg))

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


# --- 여러 모드가 함께 쓰는 시트 --------------------------------------------
def _sheet_struct_qc(qc_results, title):
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
    return {"title": title, "headers": h2, "rows": r2,
            "wrap": [len(h2) - 1], "maxw": 60,
            "note": "위치는 트리밍 후 서열 기준 1-based. 랜드마크 상태 표기"
                    "(OK / S#G# / GAP# / ABSENT / NA)는 용어설명 시트 참조."}


def _sheet_glossary(title, primers=None, cfg=None):
    return {"title": title, "headers": ["구분", "용어 / 컬럼", "설명"],
            "rows": glossary(primers, cfg), "wrap": [2], "maxw": 110,
            "note": "이 도구를 처음 보는 사람도 이 시트만으로 모든 컬럼과 판정을 "
                    "이해할 수 있도록 정리했습니다."}


def _sheet_config(cfg, meta, title, primers=None):
    h5 = ["구분", "항목", "값", "기본값", "기본값과 동일", "설명"]
    r5 = [["실행", "core 버전", cfg["core_version"], "", "", "기준 노트북 v" + cfg["nb_version"]],
          ["실행", "실행 환경", meta.get("runtime", ""), "", "", ""],
          ["실행", "실행 시각", meta.get("timestamp", ""), "", "", ""],
          ["실행", "param_hash", cfg["param_hash"], "", "",
           "판정 임계값 16개의 지문. 배치 간 비교 시 이 값이 같아야 동일 기준"],
          ["실행", "design_hash", cfg["design_hash"], "", "",
           "판정에 관여하는 설계 %d 개(배치 VH family·경쇄, fragment 1·2 For 구성)의 "
           "지문. param_hash 가 같아도 이 값이 다르면 서로 다른 기준으로 판정된 것"
           % len(JUDGMENT_DESIGN_KEYS)],
          ["실행", "프라이머 FASTA", meta.get("primer_file", ""), "", "",
           "%d 종" % meta.get("primer_n", 0)],
          ["실행", "입력 파일", ", ".join(meta.get("files", [])), "", "", ""]]
    amb = primer_ambiguity(primers, cfg) if primers else []
    if primers:
        s = ambiguity_summary(amb)
        r5.append(["실행", "모호성 쌍",
                   "%d 건 (certain %d · split %d)" % (s["total"], s["certain"], s["split"]),
                   "", "",
                   "같은 family %d · 다른 family %d. 판별구간이 겹쳐 한 read 가 두 "
                   "프라이머 모두와 허용 미스매치 안에서 맞을 수 있는 쌍입니다. "
                   "등급의 뜻은 용어설명 참조." % (s["same_family"], s["cross_family"])])
        for tier in ("certain", "split"):
            sub = [p for p in amb if p["tie"] == tier and not p["same_family"]]
            r5.append(["실행", "모호성 · family 다름 (%s)" % tier, "%d 건" % len(sub), "", "",
                       " · ".join(fmt_ambiguity_pair(p) for p in sub) or "없음"])
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
    return {"title": title, "headers": h5, "rows": r5, "wrap": [5], "maxw": 60,
            "note": "이 배치의 판정에 실제로 쓰인 값 전부입니다. "
                    "'기본값과 동일' 이 X 인 항목은 결과에 직접 영향을 줍니다."}


def insert_seq(r):
    if r["pos_notI"] >= 0 and r["pos_ascI"] > r["pos_notI"]:
        return r["seq"][r["pos_notI"]:r["pos_ascI"] + len(CONST["AscI"])]
    return ""


# =============================================================================
#  12-2. 음성 대조군
# =============================================================================
# 벡터만 ligation 한 대조군은 기대값이 정반대입니다. 일반 배치에서 PASS 는
# 성공이지만 대조군에서 완전한 scFv 가 나오면 오염 신호입니다.
# qc_one 은 손대지 않고, 그 결과를 대조군 관점으로 다시 읽기만 합니다.
#
# 분류는 "무엇이 관측됐나" 가 아니라 "벡터에 무엇이 들어갈 수 있나" 로 짰습니다.
# 스터퍼 보유(PARENTAL)는 지금 표본에 없어도 벡터 준비물을 새로 만들면 주된
# 유형이 될 수 있어 분기를 둡니다. 어디에도 맞지 않으면 억지로 끼우지 않고
# CHECK 로 남깁니다.
NEGCTRL_VERDICTS = [
    "CONCATEMER", "MIXED", "EMPTY_VECTOR", "PARENTAL", "PARENTAL?",
    "CONTAMINATED", "CONTAMINATED?", "PARTIAL_INSERT", "CARRYOVER", "CHECK",
]

# 화면 색 구분에만 쓰는 등급입니다. verdict 선택 우선순위는 negctrl_verdict 의
# 판정 순서로 고정되어 있고 이 등급을 쓰지 않습니다. FLAG_SEV 도 쓰지 않습니다.
# 빨강은 "scFv 크기의 인서트가 링커까지 갖춘 채 들어감" 에만 씁니다.
# CONTAMINATED? 의 물음표는 오염 여부가 아니라 프레임 이상을 가리킵니다.
NEGCTRL_LEVEL = {
    "CONTAMINATED": "fail",     # 완전한 scFv 가 들어감
    "CONTAMINATED?": "fail",    # scFv 가 들어갔고 프레임까지 밀림
    "EMPTY_VECTOR": "none",     # 대조군에서 기대되는 배경
    "CONCATEMER": "check", "MIXED": "check", "PARENTAL": "check",
    "PARENTAL?": "check", "PARTIAL_INSERT": "check",
    "CARRYOVER": "check", "CHECK": "check",
}


NEGCTRL_TABLE_HEADERS = ["clone", "판정", "근거", "인서트(bp)", "인서트 md5"]


def _neg_facts(r):
    ins = r["insert_bp"]
    return {
        "notI": "NotI " + ("검출" if r["pos_notI"] >= 0 else "미검출"),
        "ascI": "AscI " + ("검출" if r["pos_ascI"] >= 0 else "미검출"),
        "link": "링커 " + ("검출" if r["pos_link"] >= 0 else "미검출"),
        "ins": "인서트 " + ("%d bp" % ins if ins is not None else "계산 불가"),
        "stuf": "스터퍼 " + ("검출" if r["stuffer"] else "미검출"),
        # 프레임 값과 인서트 길이는 모든 reason 에 들어갑니다. 분류가 틀렸을 때
        # 근거만 보고도 알아챌 수 있어야 하기 때문입니다.
        "frame": ("프레임 %%3=%d (유지 조건 %d)" % (ins % 3, CONST["FRAME_MOD"])
                  if ins is not None else "프레임 판단 불가"),
    }


def negctrl_verdict(r, cfg):
    """대조군 클론 하나를 (verdict, reason) 으로 읽는다. 위에서부터 먼저 맞는 것."""
    f = _neg_facts(r)
    ins = r["insert_bp"]
    tail = [f["ins"], f["frame"]]

    def out(v, *head):
        return v, " · ".join(list(head) + tail)

    if max(r["n_notI"], r["n_ascI"], r["n_link"]) > 1:
        return out("CONCATEMER", "NotI %d / AscI %d / 링커 %d 회 검출 (각 1 회여야 함)"
                   % (r["n_notI"], r["n_ascI"], r["n_link"]))
    if r["mix_pct"] is not None and r["mix_pct"] > cfg["mix_pct"]:
        return out("MIXED", "2순위 피크 초과 위치 %.1f%% (임계 %.1f%%)"
                   % (r["mix_pct"], cfg["mix_pct"]))
    if r["pos_notI"] < 0 or r["pos_ascI"] < 0:
        return out("EMPTY_VECTOR", f["notI"], f["ascI"], f["link"])
    if r["stuffer"]:
        return out("PARENTAL", "스터퍼 서열 검출")
    if ins is not None and ins == CONST["STUFFER_INSERT_BP"]:
        return out("PARENTAL?", "인서트 길이가 스터퍼 보유 클론과 동일", f["stuf"])
    if r["pos_link"] >= 0 and ins is not None:
        rng = "범위 %d~%d" % (cfg["insert_min"], cfg["insert_max"])
        if cfg["insert_min"] <= ins <= cfg["insert_max"]:
            if ins % 3 == CONST["FRAME_MOD"]:
                return out("CONTAMINATED", f["link"], rng + " 내", "프레임 정상")
            return out("CONTAMINATED?", f["link"], rng + " 내", "프레임 이상")
        if ins < cfg["insert_min"]:
            return out("PARTIAL_INSERT", f["link"],
                       "인서트가 하한 %d bp 미만" % cfg["insert_min"])
    if ins is not None and r["pos_link"] < 0:
        return out("CARRYOVER", f["link"], f["stuf"])
    return out("CHECK", f["notI"], f["ascI"], f["link"], f["stuf"])


def insert_md5(r):
    """인서트 염기서열의 지문. 다른 모드 결과와 수기 대조하기 위한 것입니다."""
    s = insert_seq(r)
    return hashlib.md5(s.encode("ascii")).hexdigest() if s else ""


def negctrl_summary(qc_results, cfg):
    """대조군 집계. 비율이나 배경률은 내지 않습니다 (n 이 작아 정당화되지 않음)."""
    counts = dict((v, 0) for v in NEGCTRL_VERDICTS)
    clones = []
    for r in qc_results:
        v, why = negctrl_verdict(r, cfg)
        counts[v] = counts.get(v, 0) + 1
        clones.append({"id": r["id"], "verdict": v, "reason": why,
                       "insert_bp": r["insert_bp"], "insert_md5": insert_md5(r),
                       "level": NEGCTRL_LEVEL.get(v, "check")})
    return {"n_total": len(qc_results),
            "counts": [[v, counts[v]] for v in NEGCTRL_VERDICTS],
            "clones": clones,
            # 화면 표도 컬럼명을 여기서 받아 갑니다. index.html 이 따로 들고 있지
            # 않도록 01_대조군판정 시트와 같은 상수를 씁니다.
            "table": {"headers": list(NEGCTRL_TABLE_HEADERS),
                      "rows": [[c["id"], c["verdict"], c["reason"],
                                c["insert_bp"], c["insert_md5"]] for c in clones]},
            "fingerprints": sorted(set(c["insert_md5"] for c in clones if c["insert_md5"]))}


def build_negctrl_sheets(qc_results, cfg, meta, primers=None):
    """대조군 워크북. 프라이머 판별과 배치 조성은 인서트가 scFv 가 아니라 뺍니다."""
    neg = negctrl_summary(qc_results, cfg)
    by_id = dict((c["id"], c) for c in neg["clones"])

    h1 = list(NEGCTRL_TABLE_HEADERS) + \
        ["NotI위치", "링커위치", "AscI위치", "스터퍼", "프레임(%3)"]
    r1 = []
    for r in qc_results:
        c = by_id[r["id"]]
        ins = r["insert_bp"]
        r1.append([r["id"], c["verdict"], c["reason"], ins, c["insert_md5"],
                   r["pos_notI"] + 1 if r["pos_notI"] >= 0 else None,
                   r["pos_link"] + 1 if r["pos_link"] >= 0 else None,
                   r["pos_ascI"] + 1 if r["pos_ascI"] >= 0 else None,
                   _yn(r["stuffer"]), (ins % 3) if ins is not None else None])
    sheets = [{"title": "01_대조군판정", "headers": h1, "rows": r1,
               "wrap": [2], "maxw": 70,
               "note": "음성 대조군은 기대값이 반대입니다. 완전한 scFv 가 나오면 "
                       "오염입니다. 판정 코드의 뜻은 05_용어설명 참조."}]

    sheets.append(_sheet_struct_qc(qc_results, "02_구조QC상세"))

    h3 = ["구분", "항목", "값", "비고"]
    r3 = [["개요", "분석 클론 수", neg["n_total"], ""]]
    for v, n in neg["counts"]:
        r3.append(["유형", v, n, NEGCTRL_LEVEL.get(v, "")])
    for fp in neg["fingerprints"]:
        ids = [c["id"] for c in neg["clones"] if c["insert_md5"] == fp]
        r3.append(["인서트 지문", fp, ", ".join(ids),
                   "다른 모드 결과의 같은 지문과 대조하세요. 모드별로 파일이 "
                   "분리되므로 자동 대조는 하지 않습니다."])
    if not neg["fingerprints"]:
        r3.append(["인서트 지문", "-", "", "인서트가 계산된 클론이 없습니다"])
    sheets.append({"title": "03_대조군요약", "headers": h3, "rows": r3,
                   "wrap": [3], "maxw": 70,
                   "note": "관측 사실만 적습니다. 비율·배경률 추정과 해석은 "
                           "넣지 않았습니다."})

    sheets.append(_sheet_config(cfg, meta, "04_실행설정", primers))
    sheets.append(_sheet_glossary("05_용어설명", primers, cfg))

    h6 = ["clone", "판정", "인서트 길이(bp)", "인서트 염기서열 (NotI~AscI)", "인서트 md5"]
    r6 = [[r["id"], by_id[r["id"]]["verdict"], r["insert_bp"], insert_seq(r),
           by_id[r["id"]]["insert_md5"]] for r in qc_results]
    sheets.append({"title": "06_서열", "headers": h6, "rows": r6,
                   "wrap": [3], "maxw": 70,
                   "note": "IMGT 조회나 수기 대조용입니다. 인서트는 NotI 인식서열 "
                           "첫 염기부터 AscI 인식서열 마지막 염기까지입니다."})
    return sheets


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
    info = read_id_from_name(filename)
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
    # direction / dir_guessed / parsed 는 파일명 파싱 폐기와 함께 사라졌습니다.
    # 판독 방향은 qc_one 의 orient() 결과(analyze 의 qc 항목)에서 확인하세요.
    return {"id": r["id"], "filename": r["filename"], "clone": r["clone"],
            "batch": r["batch"], "date": r["date"], "primer": r["primer"],
            "raw_len": r["raw_len"], "mean_q": round(r["mean_q"], 1),
            "has_qual": r["has_qual"], "has_trace": r["has_trace"]}


def scan_inputs(files, primer_text=""):
    """폼을 채우기 위한 입력 요약. 판정은 하지 않는다."""
    reads, errors = _load_reads(files)
    primers, trims, pwarns = parse_primer_fasta(primer_text) if primer_text else ([], {}, [])
    batches = []          # 파일명에서 배치를 뽑지 않습니다. 폼에서 입력받습니다.
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
    # 배치명은 파일명이 아니라 폼에서 받습니다.
    label = str(meta.get("batch_label", "") or "")
    batches = [label] if label else []
    cfg = build_config(overrides, batches)
    if cfg["errors"]:
        return {"ok": False, "errors": cfg["errors"], "config": cfg}

    qc_results = [qc_one(r, cfg) for r in reads]
    if primers:
        amb = ambiguity_context(primers, cfg)
        calls = [call_one(r, primers, cfg, amb) for r in qc_results]
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
    coverage = (primer_coverage(calls, primers, cfg, batch_assigned(cfg))
                if primers else [])
    sheets = build_sheets(qc_results, calls, cfg, comp, meta, primers, coverage)
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
        "called": compose_called(calls, cfg),
        "coverage": coverage,
        # 대조군 모드에서만 채웁니다. 키는 항상 있어야 화면 쪽 계약이 흔들리지 않습니다.
        "negctrl": (negctrl_summary(qc_results, cfg)
                    if cfg["analysis_mode"] == MODE_NEGCTRL else None),
        "summary": {"headers": SUMMARY_HEADERS, "rows": summary_rows},
        "sheets": sheets,
        "fasta": build_fasta(qc_results, calls),
    }
