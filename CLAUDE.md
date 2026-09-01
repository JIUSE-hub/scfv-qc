# 작업 규칙

## 프로젝트
pAIM1 scFv 라이브러리의 Sanger 판독(.ab1) QC 웹도구. GitHub Pages 배포.
core.py 는 Colab 노트북 scFv_sequencing_QC.ipynb v1.0 의 판정 로직을 옮긴 것이다.

## 원칙
- 판정에 쓰이는 수치는 CFG_DEFAULTS 를 통해 UI 로 노출한다.
- index.html 에 서열·임계값·용어설명을 넣지 않는다.
  화면의 모든 과학적 내용은 실행 시 core.py 에서 읽어온다.
- 판정 결과가 노트북과 달라지는 변경은 사전에 보고하고 승인을 받는다.

## 수정 방식
- 파일 전체 재작성 금지. 바꿀 부분만 고치고 diff 를 보여준다.
- 수정 후 반드시 verify.py 를 실행하고 결과를 보고한다.
- 의도적 변경은 아무리 사소해도 전부 밝힌다. 문구 하나, 정렬 순서 하나까지.
- 판정 로직이 바뀌면 CORE_VERSION 을 올린다.
  표시·문구만 바뀌면 소수점, 판정이 달라지면 정수 자리를 올린다.
- 검사를 새로 추가하거나 예외 규칙을 넣을 때는, 그 검사가 실제로
  실패를 잡는지 고장 주입으로 확인한다. 통과만으로는 검사가
  작동한다는 증거가 되지 않는다.

## 금지
- testdata/ 의 .ab1 파일을 커밋하지 않는다.
- 프라이머 FASTA(scFv_primers.fa)를 저장소에 올리지 않는다.
- verify.py 통과 없이 커밋하지 않는다.

## 알려진 이슈

1. [버그] · 처리 완료 (core 1.0 / verify 1.1)
   index.html 의 badge() 가 check_landmark 의 WARN 상태값
   "S{sub}G{gap}" 을 처리하지 않아 FAIL 과 같은 빨간 배지로 표시된다.
   b-warn 으로 나와야 한다.

2. [문서] · 처리 완료 (core 1.1 · RULES 계층 신설)
   core.py 모듈 docstring 이 "코드에 고정된 것은 서열과 프레임 규칙뿐"
   이라고 하지만, 아래도 판정에 직접 작용한다. CONST_DOC 또는 별도 목록으로
   05_실행설정·06_용어설명 에 노출해야 한다.
     FLAG_SEV(verdict 우선순위) · _FR4_MOTIF(CDR3-H3 절단) ·
     _AL_MATCH/_AL_MIS/_AL_GAP/_AL_FLANK(정렬) ·
     _EXO_TRIM(CDR3 경계) · _REPEAT_MAX_PERIOD · _OVERLONG_ZONE
   _FR4_MOTIF 는 서열이므로 CONST 로 옮기는 것이 적절하다.

3. [일관성] · 처리 완료 (core 1.2)
   DESIGN_DEFAULTS 는 14개인데 DESIGN_DOC 은 12개다.
   rna_bone_marrow / rna_peripheral 이 빠져 있어 index.html 이 라벨과
   기본값을 자체 하드코딩하고 있다. DESIGN_DOC 에 두 줄을 추가하고
   index.html 의 하드코딩을 걷어낸다.

4. [재현성] · 처리 완료 (core 1.3 · design_hash 신설)
   param_hash 가 임계값 16개만 덮는다. 판정에 직접 관여하는
   batch_vh_family / batch_chain / f1_for_mode / f2_for_mode 는 빠져 있어,
   같은 param_hash 로 서로 다른 기준의 판정이 나올 수 있다.
   design_hash 를 별도로 추가하는 안을 검토한다.

5. [취약성] · 처리 완료 (verify 1.0 부터 [B] 에 상시 포함)
   index.html 의 COMPACT/NUMCOL/MONOCOL/BADGECOL 에
   SUMMARY_HEADERS 의 컬럼명이 문자열로 복제되어 있다.
   core 에서 컬럼명을 바꾸면 에러 없이 열이 조용히 사라진다.
   verify.py 가 이 대조를 반드시 포함해야 한다.

6. [미사용] · 처리 완료 (core 2.0 · PARENTAL? 신설)
   CONST["STUFFER_INSERT_BP"] = 386 은 판정에 쓰이지 않고
   문서 표기에만 나온다. 스터퍼 클론 판정에 실제로 쓸지, 아니면
   CONST 에서 빼고 용어설명으로 옮길지 결정한다.

7. [설계] verdict 선택이 FLAG_SEV 심각도가 아니라 flags 의 append 순서로
   결정된다. verdict = [f for f in flags if sev(f)==3][0] 이므로 같은
   심각도 안에서는 코드에서 먼저 append 된 플래그가 이긴다.
   실측 예 : 인서트 386 bp 클론은 PARENTAL? 를 심각도 3 으로 올려도
   verdict 가 LINKER_DEL 이다. 랜드마크 루프가 인서트 계산보다 앞이라서다.
   06_용어설명의 FLAG_SEV 설명도 "심각도 3 중 첫 번째" 라고만 하고
   그 순서가 무엇인지 밝히지 않는다. 기존 클론의 verdict 가 바뀔 수 있어
   결정 전 영향 범위 확인이 필요하다.

8. [진행 중] 대량 처리 전환 — 파일명 파싱 폐기, 분석 모드 분리,
   배치 단위 지정. 3 단계로 나눠 처리한다.

12. [해석] F1_For 프라이머 중 family 가 다른 certain 쌍이 8 건 있다.
    For-1-1b/For-1-3b(k=1) · For-1-1b/For-1-5(k=1) · For-1-3b/For-1-5(k=1) 등
    VH1 · VH3 · VH5 가 서로 구분되지 않는 조합이 있고, primer_max_mismatch=2
    에서 일상적으로 동점이 발생할 수 있다. split 까지 포함하면 F1_For 에서
    family 가 다른 쌍이 15 건이다.

    영향 : VH family 판정, WRONG_FAMILY 플래그, 04_배치조성의 family 분포.
    실측 예 — 260901 배치에서 VH6 배치 4 클론이 VH4|VH6 으로 판정됐다.
    이는 For-1-4b/For-1-6(k=4) · For-1-4c/For-1-6(k=3) 의 split 동점이며
    모호성으로 설명된다.

    반면 같은 배치의 VH5-VK_1 이 VH4(For-1-4c) 로 판정되어 WRONG_FAMILY 가
    붙은 것은 성격이 다르다. 판별구간 22 nt 중 21 곳이 For-1-4c 와 일치하고
    (mm1, Δ2, 불일치 1 번 위치 하나), For-1-5 는 mm4 로 멀다. VH4|VH5 조합은
    certain·split 어느 목록에도 없다(For-1-4b/For-1-5 k=6, For-1-4c/For-1-5 k=5).
    즉 모호성이 아니라 진짜 교차오염이거나 판별구간 밖의 문제이며,
    검토 방향 (c)로는 걸러지지 않는다. 별도로 다뤄야 한다.

    검토 방향 : (a) 동점·근소차 판정에 Δ 와 함께 모호성 등급을 표시한다
    (b) family 가 다른 certain 쌍이 관여한 판정은 별도 표시한다
    (c) WRONG_FAMILY 를 붙이기 전에 해당 쌍이 모호성 목록에 있는지 본다

13. [문서] 06_용어설명의 Δ · 모호성 표기 · CDR3 경계 · Over 프라이머
    네 항목이 프라이머 이름을 예시로 하드코딩하고 있다. 이슈 11 대상은
    아니었으나 같은 성격의 문제다.

14. [사각지대] · 처리 완료 (core 3.0 · AMBIG_CALL? 신설)
    배치 지정 family 가 판정 집합에 포함되지만 단독이 아닌
    클론에는 아무 플래그도 붙지 않는다. WRONG_FAMILY / AMBIG_FAMILY? 는
    'batch_vh_family not in fam_matched' 일 때만 실행되기 때문이다.

    실측 — 260901 VH6 배치의 4 클론이 VH4|VH6 동점(For-1-4b/For-1-6,
    split, k4)인데 배치 지정 VH6 가 집합 안에 있어 분기에 닿지 않는다.
    유일한 신호는 03_프라이머판별의 모호성 열이며, 01_판정요약만 보면
    family 가 확정되지 않았다는 사실을 알 수 없다.

    처리 : AMBIG_CALL? (심각도 2) 를 신설했다. 판정 후보가 2 개 이상이고
    family 가 여러 종이면 배치 지정 여부와 무관하게 붙는다. 260901 에서
    7 건 · 6 클론(VH6-VK_10 이 vh·vj 양쪽 해당)이며, 전부 이미 FAIL 이라
    통과 클론은 줄지 않았다(통과 1 건은 VH5-VK_5 이고 네 그룹 모두 후보 1).
    그래서 별도 필드가 아니라 플래그로 두었다.

15. [불일치] · 처리 완료 (core 3.0 · QC2 상태가 OK 일 때만 위치 보정)
    랜드마크 상태와 구조 플래그가 서로 다른 집합을 만든다.
    check_landmark 는 허용 치환 안에서 링커를 찾아 QC2 를 OK 로 판정하지만,
    구조 문법은 s.find(CONST["QC2"]) 로 완전 일치만 보므로 pos_link = -1 이
    되어 NO_LINKER 가 붙는다.

    실측 — 260901 에서 NO_LINKER 19 건 중 2 건(VH6-VK_1 · VH6-VK_3)이
    QC2 상태 OK 이면서 링커 위치 미검출이다.

    영향 : 01_판정요약에서 QC2 열은 OK 인데 구조QC 는 NO_LINKER 로 나와
    읽는 사람이 모순으로 느낀다. d1·d2 도 계산되지 않는다.

    처리 : s.find 가 실패하고 QC2 상태가 OK 일 때만 랜드마크가 찾은 위치를
    쓴다. GAP·ABSENT·S#G# 는 랜드마크가 깨진 것이라 보정하지 않는다.
    NotI·AscI 에는 같은 보정을 할 수 없다 — QC1 은 pos_notI 가, QC3·QC4 는
    pos_ascI 가 있어야 covered 로 검사되므로 그 자리가 비면 랜드마크도 NA 가
    되어 되살릴 근거가 없다. 즉 이 어긋남은 구조상 QC2 에서만 생긴다.
    실측 2 건 모두 d1 이 범위 안이고(360 · 384) F3_For 판별이 되살아나
    NO_VL 이 사라졌으며, verdict 는 NO_LINKER 에서 NO_ASCI 로 바뀌었다.
