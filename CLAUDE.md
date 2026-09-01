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

13. [문서] · 처리 완료 (core 3.2 · 용어설명 예시 런타임 선택)
    06_용어설명의 Δ · 모호성 표기 · CDR3 경계 · Over 프라이머
    네 항목이 프라이머 이름을 예시로 하드코딩하고 있었다. 이슈 11 대상은
    아니었으나 같은 성격의 문제다.

    처리 : 다섯 항목을 (나) 런타임 선택 / (다) 유지 로 나눴다. 새로 넣은
    [B] 검사가 다섯 번째 항목(이슈 12 때 추가한 'VH4|VH6 같은 표기의 두 가지
    뜻')도 같은 부류로 잡아내 함께 처리했고, 제목이 특정 이름에 묶여 있어
    'family 칸의 | 가 갖는 두 가지 뜻' 으로 바꿨다.

      Δ (차순위 간격)      (나) F2_Rev 에서 비호환 최소 쌍 · 이름순 tie-break
      모호성 표기          (나) 채점 그룹 · certain · 같은 family 인 쌍
      CDR3 경계            (나) 다른 family 를 축퇴 공간으로 완전히 포함하는 쌍
      family 칸의 |        (나) (a) 다표적 프라이머 이름순 첫 개
                                (b) 다른 family 인 채점 그룹 쌍 중 비호환 최소
      Over 프라이머        (다) 이름이 데이터가 아니라 설계상 고정된 역할이다.
                                GLOSSARY_NAME_OK 예외로 두되 그 이름이 실제로
                                FASTA 에 있는지 함께 검사한다.

    조건에 맞는 쌍이 없으면 예시 없이 개념만 서술한다. 선택 규칙과 폴백은
    [E] 단위시험이, 이름이 리터럴로 되박히는 것은 [B] 가 잡는다.

    'CDR3 경계' 의 '완전 포함' 주장은 검증됐다. 참이며, 원문이 든 두 쌍보다
    많다 — F1_Rev 4 쌍(Rev-1-3 ⊇ Rev-1-1 · Rev-1-4 · Rev-1-6, Rev-1-1 ⊇
    Rev-1-6, 21 nt 전체) · F2_For 4 쌍(For-2-3 ⊇ For-2-1 · For-2-4 · For-2-6,
    For-2-1 ⊇ For-2-6, 19 nt 전체) 으로 그룹당 4 쌍씩 모두 8 쌍이다.

    하드코딩 예시가 낡는다는 실증 : 모호성 표기의 예시 'IGKV1 (For3-k-5|6)'
    은 현재 프라이머 FASTA 에 없는 이름이었다. 지난 세트의 잔재가 문구에
    남아 아무도 모르는 채 계속 표시되고 있었다. 이 항목을 (나) 로 돌린
    이유가 편의가 아니라 정확성인 근거다.

    Δ 항목에서 "Δ 가 항상 1 입니다" 단정을 뺐다. 이 세트에서만 참인 값이고,
    프라이머가 추가되면 같은 그룹의 Δ 가 2 이상이 될 수 있다. 예시를 런타임에
    고르면서 단정만 남기면 문장 안에서 서로 어긋나므로 "이번 프라이머
    세트에서는 … 그 그룹의 Δ 가 작습니다" 로 바꿨다.

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
