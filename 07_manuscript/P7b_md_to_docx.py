#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P7b_md_to_docx.py — 将嵌入了图/表的手稿 MD 转换为 DOCX（含真实图片与表格）
  - 标题层级: # -> Title, ## -> Heading 1, ### -> Heading 2
  - 图片行: ![alt](relpath) -> 居中嵌入（宽 16cm，纵向图 14cm）
  - 管道表格: -> Table Grid 样式，表头加粗+底纹，9pt
  - 行内格式: **粗体**、[[n]](url) PubMed 超链接（蓝色下划线）
  - *注：...* / *Figure ...* -> 斜体灰色小字
  - 参考文献行 [n] ... -> 9pt
"""
import os
import re
import sys

import docx
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.opc.constants import RELATIONSHIP_TYPE

BASE = r"D:\projects\OS_matrix_dualKO"
_argv = sys.argv[1:]
MD = _argv[0] if len(_argv) >= 1 else os.path.join(BASE, "manuscript", "OS_manuscript_v1_cn.md")
OUT = _argv[1] if len(_argv) >= 2 else os.path.join(BASE, "manuscript", "OS_manuscript_v1_cn_embedded.docx")
RANGE = _argv[2] if len(_argv) >= 3 else "full"   # full | main | supp
MD_DIR = os.path.dirname(MD)

_all_lines = open(MD, encoding="utf-8").read().splitlines()
if RANGE == "main":
    # 截断到补充材料节之前，但保留参考文献列表
    _cut = next(i for i, L in enumerate(_all_lines) if L.startswith("## 补充材料图注"))
    _ref = next(i for i, L in enumerate(_all_lines) if L.startswith("## 参考文献"))
    lines = _all_lines[:_cut] + [""] + _all_lines[_ref:]
elif RANGE == "supp":
    _start = next(i for i, L in enumerate(_all_lines) if L.startswith("## 补充材料图注"))
    try:  # 补充材料不含正文参考文献
        _end = next(i for i, L in enumerate(_all_lines) if L.startswith("## 参考文献"))
    except StopIteration:
        _end = len(_all_lines)
    lines = _all_lines[_start:_end]
else:
    lines = _all_lines

doc = Document()

# ---------------------------------------------------------------- 页面与默认样式
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)   # A4
sec.top_margin = sec.bottom_margin = Cm(2.5)
sec.left_margin = sec.right_margin = Cm(2.2)

style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(11)
style.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
pf = style.paragraph_format
pf.space_after = Pt(6)
pf.line_spacing = 1.3

for hname, sz in [("Heading 1", 15), ("Heading 2", 13), ("Heading 3", 12)]:
    hs = doc.styles[hname]
    hs.font.name = "Times New Roman"
    hs.font.size = Pt(sz)
    hs.font.bold = True
    hs.font.color.rgb = RGBColor(0, 0, 0)
    hs.element.get_or_add_rPr()
    rf = hs.element.rPr.get_or_add_rFonts()
    rf.set(qn("w:eastAsia"), "黑体")


def add_hyperlink(paragraph, text, url, size=None):
    """在段落中加入真超链接（蓝字下划线）。"""
    r_id = paragraph.part.relate_to(url, RELATIONSHIP_TYPE.HYPERLINK, is_external=True)
    hl = OxmlElement("w:hyperlink")
    hl.set(qn("r:id"), r_id)
    run = OxmlElement("w:r")
    rPr = OxmlElement("w:rPr")
    color = OxmlElement("w:color"); color.set(qn("w:val"), "0563C1"); rPr.append(color)
    u = OxmlElement("w:u"); u.set(qn("w:val"), "single"); rPr.append(u)
    if size:
        szel = OxmlElement("w:sz"); szel.set(qn("w:val"), str(int(size * 2))); rPr.append(szel)
    run.append(rPr)
    wt = OxmlElement("w:t")
    wt.text = text
    wt.set(qn("xml:space"), "preserve")
    run.append(wt)
    hl.append(run)
    paragraph._p.append(hl)


CITE_RE = re.compile(r"\[\[(\d+)\]\]\((https?://[^)]+)\)")
BOLD_RE = re.compile(r"\*\*(.+?)\*\*")


def add_rich_paragraph(text, style=None, size=None, italic=False, color=None):
    """解析 **粗体** 与 [[n]](url) 的富文本段落。"""
    p = doc.add_paragraph(style=style)
    # 先按粗体切分，再在每段内切分引用链接
    pos = 0
    for bm in BOLD_RE.finditer(text):
        _emit_plain(p, text[pos:bm.start()], size, italic, color)
        r = p.add_run(bm.group(1))
        r.bold = True
        if size: r.font.size = Pt(size)
        if italic: r.italic = True
        if color: r.font.color.rgb = color
        pos = bm.end()
    _emit_plain(p, text[pos:], size, italic, color)
    return p


def _emit_plain(p, chunk, size, italic, color):
    pos = 0
    for m in CITE_RE.finditer(chunk):
        if m.start() > pos:
            _run(p, chunk[pos:m.start()], size, italic, color)
        r = p.add_run("[" + m.group(1) + "]")
        if size: r.font.size = Pt(size)
        if italic: r.italic = True
        if color: r.font.color.rgb = color
        add_hyperlink(p, "[" + m.group(1) + "]", m.group(2), size)
        # 移除上一个纯文本 run（链接文本已由 hyperlink 承担）
        r._r.getparent().remove(r._r)
        pos = m.end()
    if pos < len(chunk):
        _run(p, chunk[pos:], size, italic, color)


def _run(p, txt, size, italic, color):
    if not txt:
        return
    r = p.add_run(txt)
    if size: r.font.size = Pt(size)
    if italic: r.italic = True
    if color: r.font.color.rgb = color
    return r


def add_image(path_rel, alt):
    path = os.path.normpath(os.path.join(MD_DIR, path_rel))
    if not os.path.exists(path):
        print("[WARN] 图片缺失:", path)
        return
    # 纵向图（高>宽）用窄宽度，避免超页
    from PIL import Image
    with Image.open(path) as im:
        w, h = im.size
    width = Cm(13.5) if h > w else Cm(16.0)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(path, width=width)


def _set_cell_borders(cell, edges):
    """edges: {边名: 线宽(1/8 pt)}，宽 0 表示无框线"""
    tcPr = cell._tc.get_or_add_tcPr()
    for el in tcPr.findall(qn("w:tcBorders")):
        tcPr.remove(el)
    tcB = OxmlElement("w:tcBorders")
    for edge, sz in edges.items():
        e = OxmlElement("w:" + edge)
        e.set(qn("w:val"), "single" if sz else "none")
        e.set(qn("w:sz"), str(sz or 0))
        e.set(qn("w:space"), "0")
        e.set(qn("w:color"), "000000")
        tcB.append(e)
    tcPr.append(tcB)


def _make_three_line(tb):
    """三线表：顶线 1.5pt + 表头下线 0.75pt + 底线 1.5pt；无竖线、无内部横线"""
    tblPr = tb._tbl.tblPr
    for el in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(el)
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement("w:" + edge)
        e.set(qn("w:val"), "none")
        e.set(qn("w:sz"), "0")
        e.set(qn("w:space"), "0")
        e.set(qn("w:color"), "auto")
        borders.append(e)
    tblPr.append(borders)
    for cell in tb.rows[0].cells:
        _set_cell_borders(cell, {"top": 12, "bottom": 6})
    for cell in tb.rows[-1].cells:
        _set_cell_borders(cell, {"bottom": 12})


def add_table(rows):
    """rows: 已解析的单元格二维列表（含表头）→ 标准三线表"""
    n_cols = max(len(r) for r in rows)
    tb = doc.add_table(rows=len(rows), cols=n_cols)
    tb.style = "Normal Table"
    tb.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        for j in range(n_cols):
            cell = tb.cell(i, j)
            txt = row[j] if j < len(row) else ""
            cell.text = ""
            cp = cell.paragraphs[0]
            _emit_cell_rich(cp, txt, header=(i == 0))
            for r in cp.runs:
                r.font.size = Pt(9)
                if i == 0:
                    r.bold = True
    _make_three_line(tb)
    doc.add_paragraph()


def _emit_cell_rich(p, txt, header=False):
    pos = 0
    for bm in BOLD_RE.finditer(txt):
        if bm.start() > pos:
            _cell_plain(p, txt[pos:bm.start()], header)
        r = p.add_run(bm.group(1)); r.bold = True
        pos = bm.end()
    if pos < len(txt):
        _cell_plain(p, txt[pos:], header)


def _cell_plain(p, chunk, header):
    pos = 0
    for m in CITE_RE.finditer(chunk):
        if m.start() > pos:
            p.add_run(chunk[pos:m.start()])
        add_hyperlink(p, "[" + m.group(1) + "]", m.group(2), size=9)
        pos = m.end()
    if pos < len(chunk):
        p.add_run(chunk[pos:])


# ---------------------------------------------------------------- 逐行解析
i = 0
n = len(lines)
while i < n:
    line = lines[i]

    if not line.strip():
        i += 1
        continue

    # 图片行
    mimg = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", line)
    if mimg:
        add_image(mimg.group(2), mimg.group(1))
        i += 1
        continue

    # 表格块
    if line.lstrip().startswith("|"):
        rows = []
        while i < n and lines[i].lstrip().startswith("|"):
            raw = lines[i].strip().strip("|")
            cells = [c.strip() for c in raw.split("|")]
            if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):  # 跳过分隔行
                rows.append(cells)
            i += 1
        if rows:
            add_table(rows)
        continue

    # 标题
    if line.startswith("# "):
        p = doc.add_heading(line[2:].strip(), level=0)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        i += 1
        continue
    if line.startswith("## "):
        doc.add_heading(line[3:].strip(), level=1)
        i += 1
        continue
    if line.startswith("### "):
        doc.add_heading(line[4:].strip(), level=2)
        i += 1
        continue

    # 斜体注释行（*注：...* 或 *Figure ...* 或 *Compartment...*）
    stripped = line.strip()
    if stripped.startswith("*") and stripped.endswith("*") and not stripped.startswith("**"):
        p = doc.add_paragraph()
        r = p.add_run(stripped.strip("*"))
        r.italic = True
        r.font.size = Pt(9)
        r.font.color.rgb = RGBColor(0x59, 0x59, 0x59)
        i += 1
        continue

    # 参考文献行
    if re.match(r"^\[\d+\]\s", stripped):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(3)
        r = p.add_run(stripped)
        r.font.size = Pt(9)
        i += 1
        continue

    # 普通段落
    add_rich_paragraph(stripped)
    i += 1

doc.save(OUT)
print("[write] ->", OUT)
print("size:", os.path.getsize(OUT), "bytes")
