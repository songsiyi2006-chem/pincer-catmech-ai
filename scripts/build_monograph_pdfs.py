"""Render scientific Markdown into paginated PDFs with typeset mathematics.

Requirements: reportlab, matplotlib, pillow, jieba. Windows font locations are
defaults only; --font-directory accepts another directory with equivalent TTFs.
PDFs are derived presentation artifacts; Markdown remains the text authority.
"""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re

import jieba
import matplotlib
matplotlib.use("Agg")
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph,
    Spacer, PageBreak, Table, TableStyle, Image, Flowable, KeepTogether)
from reportlab.platypus.tableofcontents import TableOfContents

REPO = Path(__file__).resolve().parents[1]
NAVY = colors.HexColor("#17324D")
TEAL = colors.HexColor("#087F8C")
GRAY = colors.HexColor("#54616E")


class Matrix(Flowable):
    def __init__(self, expression):
        super().__init__()
        content = re.search(r"\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", expression, re.S).group(1)
        self.rows = [[x.strip() for x in row.strip().split("&")] for row in content.split("\\\\") if row.strip()]
        self.width, self.height = 440, 18 * len(self.rows) + 16

    def draw(self):
        c = self.canv
        c.setFillColor(NAVY)
        c.setFont("Times", 12)
        c.drawString(40, self.height / 2, "S =")
        x0, y0 = 94, self.height - 18
        for i, row in enumerate(self.rows):
            for j, entry in enumerate(row):
                c.drawCentredString(x0 + j * 29, y0 - i * 18, entry.replace("-", "−"))
        c.setStrokeColor(NAVY)
        for x, sign in ((76, 1), (258, -1)):
            c.line(x, 8, x, self.height - 4)
            c.line(x, 8, x + 6 * sign, 8)
            c.line(x, self.height - 4, x + 6 * sign, self.height - 4)
        c.drawString(290, self.height / 2, "dc/dt = S r")


class Document(BaseDocTemplate):
    def __init__(self, output, language, **kwargs):
        super().__init__(str(output), pagesize=A4, leftMargin=48, rightMargin=48,
                         topMargin=62, bottomMargin=52, **kwargs)
        self.language = language
        self.section_counter = 0
        self.addPageTemplates(PageTemplate(id="paper", frames=[Frame(
            self.leftMargin, self.bottomMargin, self.width, self.height,
            leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)], onPageEnd=self.page_header))

    def page_header(self, canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(1)
        canvas.line(48, A4[1] - 40, A4[0] - 48, A4[1] - 40)
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(GRAY)
        canvas.drawString(48, A4[1] - 30, "PINCER-CATMECH-AI  |  COMPUTATIONAL RESEARCH RECORD")
        canvas.drawRightString(A4[0] - 48, 31, str(doc.page))
        canvas.drawString(48, 31, "GFN2-xTB / ALPB toluene  |  Model evidence; no experimental validation")
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph) and flowable.style.name == "Section":
            text = flowable.getPlainText()
            if not re.match(r"^\d+\.", text):
                return
            key = "section_" + hashlib.sha256(text.encode()).hexdigest()[:12]
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, 0)
            self.notify("TOCEntry", (0, text, self.page, key))


def styles(language):
    family = "Song" if language == "ZH" else "Times"
    body = ParagraphStyle("Body", fontName=family, fontSize=10.3,
        leading=17 if language == "ZH" else 15, spaceAfter=8, textColor=NAVY,
        wordWrap=None, splitLongWords=True, autoLeading="max",
        allowWidows=0, allowOrphans=0)
    return {"body": body,
        "title": ParagraphStyle("Title", parent=body, fontName="Hei" if language == "ZH" else "Times-Bold",
            fontSize=26, leading=36, spaceAfter=24),
        "section": ParagraphStyle("Section", parent=body, fontName="Hei" if language == "ZH" else "Times-Bold",
            fontSize=17, leading=24, spaceAfter=16, keepWithNext=True),
        "subheading": ParagraphStyle("Subheading", parent=body, fontName="Hei" if language == "ZH" else "Times-Bold",
            fontSize=12, leading=18, spaceBefore=12, spaceAfter=8, keepWithNext=True),
        "small": ParagraphStyle("Small", parent=body, fontSize=8, leading=12, spaceAfter=4),
        "table": ParagraphStyle("Table", parent=body, fontSize=7.4, leading=11, spaceAfter=0),
        "caption": ParagraphStyle("Caption", parent=body, fontSize=8.3, leading=12, textColor=GRAY)}


def math_image(expression, cache, inline=False):
    expression = expression.replace("\n", " ").replace("≈", "\\approx ")
    expression = expression.replace("\\rm ", "\\mathrm ").replace("\\text{", "\\mathrm{")
    expression = re.sub(r"\\(mathbf|mathrm|mathit)\s+([A-Za-z]+)", r"\\\1{\2}", expression)
    expression = re.sub(r"\\(dot|hat|tilde|widetilde)\s+([A-Za-z])", r"\\\1{\2}", expression)
    expression = re.sub(r"\\frac(\d)(\d)", r"\\frac{\1}{\2}", expression)
    key = hashlib.sha256((str(inline) + expression).encode()).hexdigest()
    path = cache / (key + ".png")
    if not path.exists():
        mathtext.math_to_image("$" + expression.strip() + "$", str(path), dpi=240,
                               color="#17324D", prop=FontProperties(size=10.3 if inline else 12))
    with PILImage.open(path) as bitmap:
        width, height = bitmap.size
    width, height = width * 72 / 240, height * 72 / 240
    limit = 410 if inline else 489
    scale = min(1.0, limit / max(width, 1))
    return path, width * scale, height * scale


def escaped_prose(text):
    """Retain scientific signs absent from SimSun using an embedded fallback."""
    primary = pdfmetrics.getFont("Song").face.charToGlyph
    fallback = pdfmetrics.getFont("Times").face.charToGlyph
    pieces = []
    for character in text:
        encoded = html.escape(character)
        if ord(character) >= 128 and not primary.get(ord(character), 0):
            if not fallback.get(ord(character), 0):
                raise ValueError(f"No embedded font glyph for U+{ord(character):04X}")
            encoded = '<font name="Times">' + encoded + '</font>'
        pieces.append(encoded)
    return "".join(pieces)


def inline_markup(text, cache):
    items = []
    def token(value):
        key = f"ZZTOKEN{len(items)}ZZ"
        items.append(value)
        return key
    def equation(match):
        path, width, height = math_image(match.group(1), cache, inline=True)
        return token(f'<img src="{html.escape(str(path), quote=True)}" width="{width:.2f}" height="{height:.2f}" valign="-2"/>')
    text = re.sub(r"\$([^$]+)\$", equation, text)
    def link(match):
        return token('<link href="' + html.escape(match.group(2), quote=True) + '" color="#087F8C">' + escaped_prose(match.group(1)) + '</link>')
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
    text = escaped_prose(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(chr(96) + "([^" + chr(96) + "]+)" + chr(96), r'<font color="#087F8C">\1</font>', text)
    for index, value in enumerate(items):
        text = text.replace(f"ZZTOKEN{index}ZZ", value)
    return text


def render(source, output, cache):
    language = "ZH" if source.stem.endswith("ZH") else "EN"
    source_text = source.read_text(encoding="utf-8")
    lines = source_text.splitlines()
    sty = styles(language)
    title = lines[0].lstrip("# ")
    story = [Spacer(1, 80), Paragraph(escaped_prose(title), sty["title"]), Spacer(1, 18),
        Paragraph("GFN2-xTB · ALPB / Toluene · 383.15 K", sty["subheading"]),
        Paragraph("24 designed catalysts | 4 metals × 3 backbones × 2 substituents", sty["body"]),
        Spacer(1, 26), Paragraph("可追溯的计算研究记录" if language == "ZH" else "An auditable computational research record", sty["body"]),
        Paragraph("模型结果、计算失败与尚缺证据分别记录。" if language == "ZH" else "Calculated results, rejected searches and missing evidence are documented separately.", sty["body"]),
        Spacer(1, 30), Paragraph("https://github.com/songsiyi2006-chem/pincer-catmech-ai", sty["small"]),
        PageBreak(), Paragraph("目录" if language == "ZH" else "Contents", sty["section"])]
    toc = TableOfContents()
    toc.levelStyles = [ParagraphStyle("TOC", parent=sty["body"], fontSize=9.2, leading=15, spaceBefore=3)]
    story.extend([toc])
    index, paragraphs, formula_count = 1, [], 0
    def flush():
        if paragraphs:
            story.append(Paragraph(inline_markup(" ".join(paragraphs), cache), sty["body"]))
            paragraphs.clear()
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            flush(); index += 1; continue
        if line.startswith("## "):
            flush(); story.extend([PageBreak(), Paragraph(escaped_prose(line[3:]), sty["section"])])
            index += 1; continue
        if line.startswith("### "):
            flush(); story.append(Paragraph(escaped_prose(line[4:]), sty["subheading"])); index += 1; continue
        if line.startswith("$$"):
            flush()
            expression = line[2:]
            while "$$" not in expression:
                index += 1
                if index >= len(lines):
                    raise ValueError("Unclosed display equation")
                expression += "\n" + lines[index]
            expression = expression.split("$$")[0]
            if "\\begin{pmatrix}" in expression:
                equation = Matrix(expression)
            else:
                path, width, height = math_image(expression, cache)
                equation = Image(str(path), width=width, height=height)
                equation.hAlign = "CENTER"
            story.extend([Spacer(1, 5), equation, Spacer(1, 11)])
            formula_count += 1; index += 1; continue
        if line.startswith("|"):
            flush(); rows = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                row = [item.strip() for item in lines[index].strip().strip("|").split("|")]
                if not all(re.fullmatch(r"[-: ]+", item) for item in row):
                    rows.append(row)
                index += 1
            count = max(map(len, rows))
            weights = [max(5, min(30, max(len(re.sub(r"[*$]", "", row[col]))
                for row in rows if col < len(row)))) for col in range(count)]
            widths = [499 * weight / sum(weights) for weight in weights]
            table = Table([[Paragraph(inline_markup(cell, cache), sty["table"]) for cell in row] for row in rows],
                          colWidths=widths, repeatRows=1, hAlign="LEFT")
            table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#DDECF0")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7F9")]),
                ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), .7, TEAL),
                ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
            story.extend([table, Spacer(1, 12)]); continue
        image_match = re.fullmatch(r"!\[([^\]]*)\]\(([^)]+)\)", line)
        if image_match:
            flush(); image_path = (source.parent / image_match.group(2)).resolve()
            if image_path.suffix == ".svg":
                image_path = image_path.with_suffix(".png")
            with PILImage.open(image_path) as bitmap:
                width, height = bitmap.size
            story.append(KeepTogether([Image(str(image_path), width=489, height=489 * height / width),
                          Paragraph(escaped_prose(image_match.group(1)), sty["caption"])])); index += 1; continue
        if line.startswith("- "):
            flush(); story.append(Paragraph(inline_markup(line[2:], cache), sty["body"], bulletText="•"))
            index += 1; continue
        paragraphs.append(line); index += 1
    flush()
    doc = Document(output, language, title=title, author="Pincer-CatMech-AI computational campaign")
    doc.multiBuild(story)
    # Chinese words are explicitly segmented, not equated with whitespace tokens.
    prose = re.sub(r"\$\$.*?\$\$|\$[^$]+\$", " ", source_text, flags=re.S)
    words = [w for w in jieba.cut(prose) if re.search(r"[A-Za-z\u4e00-\u9fff]", w)] if language == "ZH" else re.findall(r"\b[A-Za-z]+(?:[-'][A-Za-z]+)*\b", prose)
    return {"source": str(source), "pdf": str(output), "language": language,
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "pdf_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            "word_count": len(words), "word_count_method": "jieba 0.42.1 segmented lexical tokens" if language == "ZH" else "English lexical word regex excluding equations",
            "CJK_characters": len(re.findall(r"[\u4e00-\u9fff]", source_text)),
            "display_equations": formula_count, "pages": doc.page,
            "unicode_prose_fallback": "Embedded Times for scientific Unicode signs not present in SimSun; missing fallback glyphs are fatal",
            "inline_math_layout": "autoLeading=max; tall activity and flux products are separate display equations",
            "figure_caption_layout": "KeepTogether group"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--font-directory", type=Path, default=Path("C:/Windows/Fonts"))
    args = parser.parse_args()
    fonts = args.font_directory
    for name, file in [("Times", "times.ttf"), ("Times-Bold", "timesbd.ttf"), ("Times-Italic", "timesi.ttf"),
                        ("Times-BoldItalic", "timesbi.ttf"), ("Song", "simsun.ttc"), ("Hei", "simhei.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(fonts / file), subfontIndex=0))
    pdfmetrics.registerFontFamily("Times", normal="Times", bold="Times-Bold", italic="Times-Italic", boldItalic="Times-BoldItalic")
    pdfmetrics.registerFontFamily("Song", normal="Song", bold="Hei", italic="Song", boldItalic="Hei")
    cache = REPO.parent / "pdf-math-cache"
    cache.mkdir(exist_ok=True)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    results = []
    for language in ("EN", "ZH"):
        source = REPO / "docs" / f"MONOGRAPH_BORROWING_HYDROGEN_{language}.md"
        results.append(render(source, args.output_directory / (source.stem + ".pdf"), cache))
    (args.output_directory / "PDF_BUILD_AUDIT.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
