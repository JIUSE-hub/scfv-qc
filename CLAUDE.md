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

## 금지
- testdata/ 의 .ab1 파일을 커밋하지 않는다.
- 프라이머 FASTA(scFv_primers.fa)를 저장소에 올리지 않는다.
- verify.py 통과 없이 커밋하지 않는다.

## 알려진 이슈 (verify.py 구축 후 순서대로 처리)

1. [버그] index.html 의 badge() 가 check_landmark 의 WARN 상태값
   "S{sub}G{gap}" 을 처리하지 않아 FAIL 과 같은 빨간 배지로 표시된다.
   b-warn 으로 나와야 한다.

2. [문서] core.py 모듈 docstring 이 "코드에 고정된 것은 서열과 프레임 규칙뿐"
   이라고 하지만, 아래도 판정에 직접 작용한다. CONST_DOC 또는 별도 목록으로
   05_실행설정·06_용어설명 에 노출해야 한다.
     FLAG_SEV(verdict 우선순위) · _FR4_MOTIF(CDR3-H3 절단) ·
     _AL_MATCH/_AL_MIS/_AL_GAP/_AL_FLANK(정렬) ·
     _EXO_TRIM(CDR3 경계) · _REPEAT_MAX_PERIOD · _OVERLONG_ZONE
   _FR4_MOTIF 는 서열이므로 CONST 로 옮기는 것이 적절하다.

3. [일관성] DESIGN_DEFAULTS 는 14개인데 DESIGN_DOC 은 12개다.
   rna_bone_marrow / rna_peripheral 이 빠져 있어 index.html 이 라벨과
   기본값을 자체 하드코딩하고 있다. DESIGN_DOC 에 두 줄을 추가하고
   index.html 의 하드코딩을 걷어낸다.

4. [재현성] param_hash 가 임계값 16개만 덮는다. 판정에 직접 관여하는
   batch_vh_family / batch_chain / f1_for_mode / f2_for_mode 는 빠져 있어,
   같은 param_hash 로 서로 다른 기준의 판정이 나올 수 있다.
   design_hash 를 별도로 추가하는 안을 검토한다.

5. [취약성] index.html 의 COMPACT/NUMCOL/MONOCOL/BADGECOL 에
   SUMMARY_HEADERS 의 컬럼명이 문자열로 복제되어 있다.
   core 에서 컬럼명을 바꾸면 에러 없이 열이 조용히 사라진다.
   verify.py 가 이 대조를 반드시 포함해야 한다.

6. [미사용] CONST["STUFFER_INSERT_BP"] = 386 은 판정에 쓰이지 않고
   문서 표기에만 나온다. 스터퍼 클론 판정에 실제로 쓸지, 아니면
   CONST 에서 빼고 용어설명으로 옮길지 결정한다.
