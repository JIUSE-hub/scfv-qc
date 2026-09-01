# pAIM1 scFv Library — Sanger Read QC (웹 버전)

`.ab1` 판독 파일에서 scFv 카세트의 조립·클로닝 정합성을 검사하고,
유래 프라이머를 역추적해 xlsx / csv / fasta 로 내보냅니다.

브라우저 안에서 전부 처리되며 **업로드한 파일은 어디로도 전송되지 않습니다.**

## 파일 구성

| 파일 | 역할 |
|---|---|
| `index.html` | 화면과 파일 입출력. **판정 로직·서열·임계값 없음** |
| `core.py` | 판정 로직 전부 |
| `xlsx_writer.py` | `.xlsx` / `.csv` 생성. 표준 라이브러리만 사용 |
| `verify.py` | 검증 스크립트. 구문·계약·누출·회귀 98 개 검사 |
| `CLAUDE.md` | 작업 규칙과 알려진 이슈 |
| `scFv_primers.fa` | 프라이머 목록. 보안상 저장소에 포함하지 않습니다. 사용 시 직접 업로드하세요 |

외부 파이썬 패키지 의존성 **0 개**. `hashlib, io, json, re, struct, unicodedata, zipfile` 만 사용합니다.

## GitHub Pages 배포

1. 저장소 루트에 `index.html` · `core.py` · `xlsx_writer.py` 를 올립니다.
2. 빈 파일 `.nojekyll` 을 만듭니다 (Jekyll 이 `.py` 를 건드리지 않게).
3. `Settings → Pages → Source` 를 `Deploy from a branch` / `main` / `/ (root)` 로 지정합니다.
4. `https://<사용자명>.github.io/<저장소명>/` 으로 접속합니다.

> `index.html` 을 더블클릭해 `file://` 로 여는 방식은 브라우저의 `fetch` 제약 때문에
> 동작하지 않을 수 있습니다. 정적 호스팅에서 열어 주세요.

## Pyodide 벤더링 (선택)

기본값은 jsDelivr CDN 에서 Pyodide 를 내려받습니다(최초 1회 10~20 MB, 이후 브라우저 캐시).
완전 오프라인으로 쓰려면 `pyodide-core-314.0.5.tar.bz2` 를 받아 `pyodide/` 폴더에 풀고,
`index.html` 의 아래 한 줄만 바꾸면 됩니다.

```js
const PYODIDE_BASE = "https://cdn.jsdelivr.net/pyodide/v314.0.5/full/";
// → const PYODIDE_BASE = "pyodide/";
```

우리는 표준 라이브러리만 쓰므로 전체 배포판(200 MB+)이 아니라 `pyodide-core` 로 충분합니다.

## 입력

### `.ab1` 파일

**파일명 규칙은 없습니다.** 확장자를 뗀 파일명이 그대로 클론 ID 가 됩니다.
같은 이름이 둘 이상이면 결과가 잘못 짝지어지므로 경고가 표시됩니다.

방향(정방향/역방향)은 파일명이 아니라 서열로 판정합니다.

배치명과 날짜는 설정 탭에서 입력하며, 결과 파일명과 `05_실행설정` 시트에 기록됩니다.

### 프라이머 FASTA

```
>For3-k-16 | group=F3_For | chain=kappa | family=IGKV4 | target=IGKV4-1 | fragment=K3 | dir=F | tm=52-54
TCAGGGGGCGGTGGATCCGACATCGTGATGACCCA
```

프라이머를 추가할 때는 이 형식으로 두 줄을 넣으면 됩니다. 코드 수정은 필요 없습니다.
판별구간(그룹 내 공통 5' 태그를 뺀 나머지)과 프라이머 간 모호성은 실행 시 자동으로 다시 계산됩니다.

## 분석 모드

| 모드 | 용도 | 지정 |
|---|---|---|
| **Assigned batch** | VH family 와 경쇄를 아는 배치를 검증 | 배치 카드마다 VH · 경쇄 · 파일 |
| **Library** | 라이브러리 전체. 프라이머로만 판정 | 없음 |
| **Negative control** | 벡터만 ligation 한 음성 대조군 | 없음 |

**Assigned batch** 는 배치 카드를 자유롭게 추가·삭제할 수 있고, 배치마다
`core.analyze` 를 따로 호출한 뒤 결과를 합칩니다. 지정한 family 와 다르게
판정되면 `WRONG_FAMILY`, 프라이머 모호성으로 설명되는 경우는 `AMBIG_FAMILY?` 가
붙습니다. 커버리지도 배치에 지정된 프라이머만 기대 대상으로 세므로,
지정하지 않은 family 나 경쇄의 프라이머가 dropout 으로 잘못 보고되지 않습니다.

**Negative control** 은 판정 어휘가 다릅니다. `PASS` 가 성공이 아니라
"인서트가 들어갔다"는 오염 신호입니다.

```
EMPTY_VECTOR    NotI 또는 AscI 소실. 클로닝 자리가 파괴된 재결합
PARENTAL        스터퍼 검출. 미절단 또는 단일절단 벡터
PARENTAL?       인서트 길이만 스터퍼와 같고 서열은 미검출
CONTAMINATED    링커 검출 + 인서트 범위 내 + 프레임 정상. ligation 단계 오염
CONTAMINATED?   위와 같으나 프레임 이상
PARTIAL_INSERT  링커 검출 + 인서트가 하한 미만
CARRYOVER       인서트는 있으나 링커 없음. 벡터 준비물 오염 의심
CONCATEMER · MIXED · CHECK
```

모드는 `design_hash` 에 기록되므로 결과 파일만 보고도 어느 모드로 판정했는지 알 수 있습니다.

## 결과 파일

| 파일 | 내용 |
|---|---|
| `*_scFvQC.xlsx` | 마스터 문서 |
| `*_scFvQC_summary.csv` | 판정요약 시트 (UTF-8 BOM) |
| `*_scFvQC_inserts.fa` | 통과 클론의 인서트 염기서열 |

xlsx 시트 구성은 모드에 따라 다릅니다.

```
Assigned batch / Library
  01_판정요약 · 02_구조QC상세 · 03_프라이머판별 · 04_배치조성
  05_실행설정 · 06_용어설명 · 07_서열

Negative control
  01_대조군판정 · 02_구조QC상세 · 03_대조군요약
  04_실행설정 · 05_용어설명 · 06_서열
```

용어설명 시트에 모든 컬럼·지표·판정코드의 정의가 정리되어 있어, 결과 파일만
받은 사람도 해석할 수 있습니다.

### 결과 파일 공유 시 주의

`03_프라이머판별`(프라이머 이름·좌표)과 `07_서열`(인서트 염기서열)을 합치면
프라이머 판별구간을 복원할 수 있습니다. 외부에 공유할 때는 두 시트 중 하나를
삭제하거나 csv 만 전달하세요.

## 재현성

결과 파일에 두 가지 지문이 기록됩니다.

| | 대상 | 언제 같아야 하나 |
|---|---|---|
| `param_hash` | 판정 임계값 16 개 | 두 배치를 비교하려면 같아야 함 |
| `design_hash` | 판정에 관여하는 설계 5 개 | 배치마다 달라도 정상 |

`design_hash` 대상은 `analysis_mode` · `batch_vh_family` · `batch_chain` ·
`f1_for_mode` · `f2_for_mode` 입니다. 같은 `param_hash` 라도 이 값이 다르면
서로 다른 기준으로 판정된 것입니다.

코드에 고정된 값은 두 계층입니다.

```
CONST   서열과 프레임 규칙
RULES   알고리즘 규칙 (정렬 점수 · 탐색 범위 · 플래그 심각도)
```

둘 다 UI 로 바꿀 수 없고 `param_hash` 에도 들어가지 않으며, `05_실행설정` 시트와
화면 하단에 전체 목록이 기록됩니다. 판정 임계값 16 개만 UI 로 조절하며 이 값들이
`param_hash` 를 이룹니다.

## 검증

```bash
python verify.py
```

98 개 검사가 돌고 하나라도 실패하면 종료코드 1 을 반환합니다.

| 구분 | 내용 |
|---|---|
| `[A]` | 구문 — AST · JS · 포맷 문자열 · 미정의 이름 |
| `[B]` | 계약 — index.html 이 참조하는 키가 core.py 에 실재하는가 |
| `[C]` | 누출 — index.html 에 서열 · 임계값 · 라벨이 하드코딩되지 않았는가 |
| `[D]` | 회귀 — 실측 4 클론 |
| `[E]` | 단위 — 합성 서열로 실측이 못 덮는 분기 검증 |
| `[F]` | GLUE — 배치별 호출과 병합 |
| `[G]` | 회귀 — 음성 대조군 3 클론 |
| `[H]` | 회귀 — 49 클론 배치 통계 |

`testdata/` 는 `.gitignore` 로 제외되므로 검증을 돌리려면 `.ab1` 파일이 따로
필요합니다. 파일이 없으면 `[D][G][H]` 는 건너뜁니다.

검사를 새로 추가하거나 예외 규칙을 넣을 때는 **고장 주입으로 그 검사가 실제로
실패를 잡는지 확인합니다.** 통과만으로는 검사가 작동한다는 증거가 되지 않습니다
(자세한 내용은 `CLAUDE.md`).

## 노트북과의 관계

`core.py` 는 Colab 노트북 `scFv_sequencing_QC.ipynb` v1.0 의 Cell 1/3/4/5 판정
로직에서 출발했습니다. 추출 당시의 차이는 세 가지였고 판정 결과에는 영향이
없었습니다.

- 전역 `CFG` / `PRIMERS` → 함수 인자 `cfg` / `primers`
- Biopython `SeqIO(abi)` / `Seq.translate` → 표준 라이브러리 구현
  (실측 `.ab1` 4 개 완전 일치, IUPAC 3375 코돈 전수 대조 일치)
- 화면 출력 제거, 결과를 JSON 직렬화 가능한 값으로 반환

이후 `core.py` 는 노트북과 별개로 발전했습니다 — 분석 모드 3 종, `RULES` 계층,
프라이머 모호성 계산, 배치 지정 기반 커버리지 등이 추가되었습니다. 변경 이력은
`CLAUDE.md` 의 「알려진 이슈」에 있습니다.

**판정 임계값은 여전히 노트북과 같습니다** — `param_hash` 가 `3473927a` 로
유지되며, 동일 입력에서 구조 QC 와 CDR3-H3 중앙값도 일치합니다.
