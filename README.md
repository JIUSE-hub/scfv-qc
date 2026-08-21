# pAIM1 scFv Library — Sanger Read QC (웹 버전)

`.ab1` 판독 파일에서 scFv 카세트의 조립·클로닝 정합성을 검사하고,
유래 프라이머를 역추적해 xlsx / csv / fasta 로 내보냅니다.

브라우저 안에서 전부 처리되며 **업로드한 파일은 어디로도 전송되지 않습니다.**

## 파일 구성

| 파일 | 역할 |
|---|---|
| `index.html` | 화면과 파일 입출력. **판정 로직·서열·임계값 없음** |
| `core.py` | 판정 로직 전부. 노트북 `scFv_sequencing_QC.ipynb` v1.0 에서 추출 |
| `xlsx_writer.py` | `.xlsx` / `.csv` 생성. 표준 라이브러리만 사용 |
| `scFv_primers.fa` | 프라이머 목록. 보안상 저장소에 포함하지 않습니다. 사용 시 직접 업로드하세요 |

외부 파이썬 패키지 의존성 **0 개**. `hashlib, io, json, re, struct, unicodedata, zipfile` 만 사용합니다.

## GitHub Pages 배포

1. 이 폴더의 4개 파일을 저장소 루트에 올립니다.
2. `Settings → Pages → Source` 를 `Deploy from a branch` / `main` / `/ (root)` 로 지정합니다.
3. `https://<사용자명>.github.io/<저장소명>/` 으로 접속합니다.

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

## 노트북과의 관계

`core.py` 는 노트북 Cell 1/3/4/5 의 판정 로직을 그대로 옮긴 것입니다.
차이는 세 가지뿐이며 판정 결과에는 영향이 없습니다.

- 전역 `CFG` / `PRIMERS` → 함수 인자 `cfg` / `primers`
- Biopython `SeqIO(abi)` / `Seq.translate` → 표준 라이브러리 구현
  (실측 `.ab1` 4개 완전 일치, IUPAC 3375 코돈 전수 대조 일치)
- 화면 출력 제거, 결과를 JSON 직렬화 가능한 값으로 반환

동일 입력에서 노트북과 `param_hash`, 구조 QC, CDR3-H3 중앙값이 일치함을 확인했습니다.

## 결과 파일

| 파일 | 내용 |
|---|---|
| `*_scFvQC.xlsx` | 7 시트 — 판정요약 / 구조QC상세 / 프라이머판별 / 배치조성 / 실행설정 / 용어설명 / 서열 |
| `*_scFvQC_summary.csv` | 판정요약 시트 (UTF-8 BOM) |
| `*_scFvQC_inserts.fa` | 통과 클론의 인서트 염기서열 |

`06_용어설명` 시트에 모든 컬럼·지표·판정코드의 정의가 정리되어 있어,
결과 파일만 받은 사람도 해석할 수 있습니다.

## 입력 규격

```
YYMMDD_배치명_c클론번호_프라이머명.ab1
예)  260819_VH6-VK_c01_pAIM1-seq-For.ab1
```

프라이머 FASTA 헤더 형식:

```
>For3-k-16 | group=F3_For | chain=kappa | family=IGKV4 | target=IGKV4-1 | fragment=K3 | dir=F | tm=52-54
TCAGGGGGCGGTGGATCCGACATCGTGATGACCCA
```

프라이머를 추가할 때는 이 형식으로 두 줄을 넣으면 됩니다. 코드 수정은 필요 없습니다.
