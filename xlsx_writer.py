# -*- coding: utf-8 -*-
"""
xlsx_writer.py — 표준 라이브러리만으로 .xlsx 를 쓴다
==============================================================================
.xlsx 는 OOXML(ECMA-376) 문서 몇 장을 zip 으로 묶은 것이고, zipfile 과
문자열 조립은 모두 표준 라이브러리에 있습니다. 그래서 openpyxl 없이도
동일한 결과물을 만들 수 있고, 브라우저(Pyodide)에서도 그대로 돕니다.

core.build_sheets() 가 만든 시트 데이터를 그대로 받습니다.
    sheets = [{"title", "headers", "rows", "wrap", "maxw", "note"}, ...]
    data   = write_xlsx(sheets)      # bytes

서식은 노트북(openpyxl) 판과 동일하게 맞췄습니다.
    안내문   Arial 9 이탤릭 회색
    헤더     Arial 10 볼드 흰색 / 남색 배경 / 가운데 정렬 / 줄바꿈 / 테두리
    본문     Arial 10 / 위쪽 정렬 / 테두리 / wrap 지정 열만 줄바꿈
    틀 고정  헤더 행 아래
    열 너비  내용 표시폭 기준 자동 (한글은 2칸으로 계산)

동일 입력이면 항상 동일한 바이트가 나오도록 zip 타임스탬프를 고정했습니다.
"""

import io
import re
import unicodedata
import zipfile

WRITER_VERSION = "1.0"

# 스타일 인덱스 (styles.xml 의 cellXfs 순서와 일치)
_S_DEFAULT = 0
_S_NOTE = 1
_S_HEADER = 2
_S_BODY = 3
_S_BODY_WRAP = 4

_FIXED_TIME = (1980, 1, 1, 0, 0, 0)      # zip 타임스탬프 고정 (재현성)

_NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_NS_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
_CT_SHEET = ("application/vnd.openxmlformats-officedocument"
             ".spreadsheetml.worksheet+xml")
_CT_BOOK = ("application/vnd.openxmlformats-officedocument"
            ".spreadsheetml.sheet.main+xml")
_CT_STYLES = ("application/vnd.openxmlformats-officedocument"
              ".spreadsheetml.styles+xml")

# XML 1.0 이 허용하지 않는 제어문자
_BAD_CTRL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")


# --- 저수준 유틸 --------------------------------------------------------------
def _esc(text):
    """XML 텍스트 이스케이프."""
    t = _BAD_CTRL.sub("", str(text))
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _esc_attr(text):
    return _esc(text).replace('"', "&quot;")


def col_letter(n):
    """1 -> A, 26 -> Z, 27 -> AA"""
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def display_width(v):
    """동아시아 문자를 2칸으로 세는 표시폭."""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1
               for c in str(v))


_BAD_TITLE = re.compile(r"[\[\]:*?/\\]")


def safe_title(name, used):
    """엑셀 시트명 규칙: 31자 이내, []:*?/\\ 금지, 중복 불가."""
    t = _BAD_TITLE.sub("_", str(name)).strip().strip("'")
    if not t:
        t = "Sheet"
    t = t[:31]
    base, i = t, 2
    while t in used:
        suffix = "_" + str(i)
        t = base[:31 - len(suffix)] + suffix
        i += 1
    used.add(t)
    return t


# --- 셀 --------------------------------------------------------------------
def _cell_xml(ref, style, value):
    if value is None or value == "":
        return '<c r="%s" s="%d"/>' % (ref, style)
    if isinstance(value, bool):
        return '<c r="%s" s="%d" t="b"><v>%d</v></c>' % (ref, style, 1 if value else 0)
    if isinstance(value, int):
        return '<c r="%s" s="%d"><v>%d</v></c>' % (ref, style, value)
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return '<c r="%s" s="%d" t="inlineStr"><is><t>%s</t></is></c>' % (
                ref, style, _esc(value))
        return '<c r="%s" s="%d"><v>%r</v></c>' % (ref, style, value)
    s = str(value)
    space = ' xml:space="preserve"' if (s != s.strip()) else ""
    return '<c r="%s" s="%d" t="inlineStr"><is><t%s>%s</t></is></c>' % (
        ref, style, space, _esc(s))


# --- 시트 --------------------------------------------------------------------
def _sheet_xml(sheet):
    headers = list(sheet.get("headers") or [])
    rows = list(sheet.get("rows") or [])
    wrap = set(sheet.get("wrap") or [])
    maxw = int(sheet.get("maxw") or 60)
    note = sheet.get("note")

    ncol = max([len(headers)] + [len(r) for r in rows]) if (headers or rows) else 1
    top = 3 if note else 1
    nrow = top + len(rows)

    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="%s" xmlns:r="%s">' % (_NS_MAIN, _NS_R),
           '<dimension ref="A1:%s%d"/>' % (col_letter(ncol), max(nrow, 1)),
           '<sheetViews><sheetView workbookViewId="0">',
           '<pane ySplit="%d" topLeftCell="A%d" activePane="bottomLeft" state="frozen"/>'
           % (top, top + 1),
           '</sheetView></sheetViews>',
           '<sheetFormatPr defaultRowHeight="15"/>']

    # 열 너비 : 헤더와 모든 셀의 표시폭 최대값
    cols = []
    for j in range(1, ncol + 1):
        w = display_width(headers[j - 1]) if j - 1 < len(headers) else 0
        for r in rows:
            if j - 1 < len(r) and r[j - 1] is not None:
                w = max(w, display_width(r[j - 1]))
        w = min(max(w + 2, 9), maxw)
        cols.append('<col min="%d" max="%d" width="%d" customWidth="1"/>' % (j, j, w))
    if cols:
        out.append("<cols>" + "".join(cols) + "</cols>")

    out.append("<sheetData>")
    if note:
        cells = [_cell_xml("A1", _S_NOTE, note)]
        for j in range(2, ncol + 1):
            cells.append('<c r="%s1" s="%d"/>' % (col_letter(j), _S_DEFAULT))
        out.append('<row r="1">' + "".join(cells) + "</row>")
        out.append('<row r="2"/>')
    if headers:
        cells = [_cell_xml(col_letter(j) + str(top), _S_HEADER,
                           headers[j - 1] if j - 1 < len(headers) else None)
                 for j in range(1, ncol + 1)]
        out.append('<row r="%d">' % top + "".join(cells) + "</row>")
    for i, row in enumerate(rows):
        rn = top + 1 + i
        cells = []
        for j in range(1, ncol + 1):
            v = row[j - 1] if j - 1 < len(row) else None
            st = _S_BODY_WRAP if (j - 1) in wrap else _S_BODY
            cells.append(_cell_xml(col_letter(j) + str(rn), st, v))
        out.append('<row r="%d">' % rn + "".join(cells) + "</row>")
    out.append("</sheetData></worksheet>")
    return "".join(out)


# --- 고정 파트 ---------------------------------------------------------------
_STYLES_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<styleSheet xmlns="%s">'
    '<fonts count="3">'
    '<font><sz val="10"/><color theme="1"/><name val="Arial"/><family val="2"/></font>'
    '<font><i/><sz val="9"/><color rgb="FF555555"/><name val="Arial"/><family val="2"/></font>'
    '<font><b/><sz val="10"/><color rgb="FFFFFFFF"/><name val="Arial"/><family val="2"/></font>'
    '</fonts>'
    '<fills count="3">'
    '<fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill>'
    '<fill><patternFill patternType="solid">'
    '<fgColor rgb="FF1F4E79"/><bgColor indexed="64"/></patternFill></fill>'
    '</fills>'
    '<borders count="2">'
    '<border><left/><right/><top/><bottom/><diagonal/></border>'
    '<border>'
    '<left style="thin"><color rgb="FFBFBFBF"/></left>'
    '<right style="thin"><color rgb="FFBFBFBF"/></right>'
    '<top style="thin"><color rgb="FFBFBFBF"/></top>'
    '<bottom style="thin"><color rgb="FFBFBFBF"/></bottom>'
    '<diagonal/></border>'
    '</borders>'
    '<cellStyleXfs count="1">'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>'
    '</cellStyleXfs>'
    '<cellXfs count="5">'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="2" fillId="2" borderId="1" xfId="0"'
    ' applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1">'
    '<alignment horizontal="center" vertical="center" wrapText="1"/></xf>'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0"'
    ' applyFont="1" applyBorder="1" applyAlignment="1">'
    '<alignment vertical="top"/></xf>'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0"'
    ' applyFont="1" applyBorder="1" applyAlignment="1">'
    '<alignment vertical="top" wrapText="1"/></xf>'
    '</cellXfs>'
    '<cellStyles count="1">'
    '<cellStyle name="Normal" xfId="0" builtinId="0"/>'
    '</cellStyles>'
    '</styleSheet>'
) % _NS_MAIN

_ROOT_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="%s">'
    '<Relationship Id="rId1" Type="%s/officeDocument" Target="xl/workbook.xml"/>'
    '</Relationships>'
) % (_NS_PKG_REL, _NS_R.rsplit("/relationships", 1)[0] + "/relationships")


def _workbook_xml(titles):
    sheets = "".join(
        '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (_esc_attr(t), i, i)
        for i, t in enumerate(titles, start=1))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="%s" xmlns:r="%s"><sheets>%s</sheets></workbook>'
            % (_NS_MAIN, _NS_R, sheets))


def _workbook_rels_xml(n):
    rels = "".join(
        '<Relationship Id="rId%d" Type="%s/worksheet" Target="worksheets/sheet%d.xml"/>'
        % (i, _NS_R, i) for i in range(1, n + 1))
    rels += ('<Relationship Id="rId%d" Type="%s/styles" Target="styles.xml"/>'
             % (n + 1, _NS_R))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="%s">%s</Relationships>' % (_NS_PKG_REL, rels))


def _content_types_xml(n):
    over = "".join(
        '<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="%s"/>'
        % (i, _CT_SHEET) for i in range(1, n + 1))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-'
            'package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="%s"/>'
            '%s'
            '<Override PartName="/xl/styles.xml" ContentType="%s"/>'
            '</Types>' % (_CT_BOOK, over, _CT_STYLES))


# --- 진입점 ------------------------------------------------------------------
def write_xlsx(sheets):
    """시트 데이터 리스트를 받아 .xlsx 바이트를 반환한다."""
    if not sheets:
        raise ValueError("시트가 없습니다")
    used = set()
    titles = [safe_title(s.get("title", "Sheet"), used) for s in sheets]

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        def put(name, text):
            info = zipfile.ZipInfo(name, date_time=_FIXED_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            z.writestr(info, text.encode("utf-8"))

        put("[Content_Types].xml", _content_types_xml(len(sheets)))
        put("_rels/.rels", _ROOT_RELS)
        put("xl/workbook.xml", _workbook_xml(titles))
        put("xl/_rels/workbook.xml.rels", _workbook_rels_xml(len(sheets)))
        put("xl/styles.xml", _STYLES_XML)
        for i, sh in enumerate(sheets, start=1):
            put("xl/worksheets/sheet%d.xml" % i, _sheet_xml(sh))
    return buf.getvalue()


def build_csv(headers, rows, bom=True):
    """엑셀에서 한글이 깨지지 않도록 UTF-8 BOM 을 붙인 CSV 바이트."""
    def field(v):
        if v is None:
            return ""
        s = str(v)
        if any(ch in s for ch in ',"\n\r'):
            return '"' + s.replace('"', '""') + '"'
        return s

    lines = [",".join(field(h) for h in headers)]
    for r in rows:
        lines.append(",".join(field(v) for v in r))
    text = "\r\n".join(lines) + "\r\n"
    return (u"\ufeff" + text).encode("utf-8") if bom else text.encode("utf-8")
