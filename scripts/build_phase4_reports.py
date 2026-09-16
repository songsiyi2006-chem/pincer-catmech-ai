"""Render the bilingual Phase4 report; Markdown and native data are authoritative."""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, PageBreak,
                               Table, TableStyle, Image, KeepTogether, CondPageBreak)
from PIL import Image as PILImage
from build_monograph_pdfs import inline_markup, escaped_prose, math_image

REPO = Path(__file__).resolve().parents[1]
NAVY = colors.HexColor('#18334C')
TEAL = colors.HexColor('#087F8C')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render(source, output, cache):
    zh = source.stem.endswith('ZH')
    family = 'Song' if zh else 'Times'
    body = ParagraphStyle('body', fontName=family, fontSize=10.4,
                          leading=16.6 if zh else 14.6, textColor=NAVY, spaceAfter=8,
                          # ReportLab's CJK splitter fails on inline math-image
                          # fragments. The regular long-word splitter handles
                          # both Unicode prose and those images, as in the
                          # existing bilingual monograph renderer.
                          autoLeading='max', wordWrap=None, splitLongWords=True,
                          allowWidows=0, allowOrphans=0)
    heading = ParagraphStyle('heading', parent=body, fontName='Hei' if zh else 'Times-Bold',
                             fontSize=16, leading=21, spaceAfter=12, keepWithNext=True)
    subheading = ParagraphStyle('subheading', parent=heading, fontSize=12, leading=17, spaceAfter=9)
    title = ParagraphStyle('title', parent=heading, fontSize=25, leading=34, spaceAfter=24)
    table_style = ParagraphStyle('table', parent=body, fontSize=9.1,
                                leading=14.5 if zh else 12, spaceAfter=0)
    small = ParagraphStyle('small', parent=body, fontSize=9, leading=13)
    lines = source.read_text(encoding='utf-8').splitlines()
    if '<!-- LOCAL_PILOT_RESULTS -->' in '\n'.join(lines) or '<!-- PUBLIC_STRUCTURE_RESULTS -->' in '\n'.join(lines):
        raise ValueError('Unresolved result placeholders')
    cover_title = ('从机理原型到可证伪研究<br/>Pincer-CatMech-AI<br/>第四阶段技术报告' if zh else
                   'From a mechanism prototype<br/>to falsifiable research<br/>Pincer-CatMech-AI<br/>Phase 4 technical report')
    story = [Spacer(1, 65), Paragraph(cover_title, title),
             Paragraph('PINCER-CATMECH-AI / PHASE 4 / 16 SEPTEMBER 2026', small),
             Spacer(1, 20), Paragraph('创新假说、数理模型与实验桥接' if zh else
                 'Innovation hypotheses, mathematical models and experimental bridges', heading),
             Paragraph('以 CNS 级原创性与严谨性规划研究；当前成果未达到主刊投稿标准。' if zh else
                 'A programme with CNS-level ambitions; not a claim of flagship-journal readiness.', body),
             Spacer(1, 15), Paragraph('公开结构 · 真实试算 · 失败记录 · 可核验输入' if zh else
                 'Public structures / native calculations / retained failures / auditable inputs', body),
             Spacer(1, 36), Paragraph('https://github.com/songsiyi2006-chem/pincer-catmech-ai', small), PageBreak()]
    paragraphs = []
    formulas = 0
    def flush():
        if paragraphs:
            prose = ' '.join(paragraphs)
            block = Paragraph(inline_markup(prose, cache), body)
            if prose.rstrip().endswith((':', '：')):
                block.keepWithNext = True
            story.append(block)
            paragraphs.clear()
    i = 1
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            flush(); i += 1; continue
        if line.startswith('## '):
            flush()
            story.extend([CondPageBreak(120), Spacer(1, 16), Paragraph(escaped_prose(line[3:]), heading)])
            i += 1; continue
        if line.startswith('### '):
            flush()
            story.extend([CondPageBreak(90), Spacer(1, 10), Paragraph(escaped_prose(line[4:]), subheading)])
            i += 1; continue
        if line.startswith('$$'):
            flush()
            expression = line[2:]
            while '$$' not in expression:
                i += 1
                expression += '\n' + lines[i]
            expression = expression.split('$$')[0]
            path, width, height = math_image(expression, cache)
            figure = Image(str(path), width=width, height=height)
            figure.hAlign = 'CENTER'
            lead_space = Spacer(1, 3)
            lead_space.keepWithNext = True
            story.extend([lead_space, figure, Spacer(1, 10)])
            formulas += 1
            i += 1; continue
        if line.startswith('|'):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                row = [s.strip() for s in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[-: ]+', s) for s in row):
                    rows.append(row)
                i += 1
            n = len(rows[0])
            weights = [max(9, min(38, max(len(re.sub(r'[*$]', '', r[j])) for r in rows))) for j in range(n)]
            widths = [499*x/sum(weights) for x in weights]
            table = Table([[Paragraph(inline_markup(s, cache), table_style) for s in r] for r in rows],
                          colWidths=widths, repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DBEBEF')),
                ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F3F6F8')]),
                ('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.8,TEAL),
                ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),
                ('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
            story.extend([table, Spacer(1,12)])
            continue
        match = re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if match:
            flush()
            path = (source.parent / match.group(2)).resolve()
            with PILImage.open(path) as im:
                w,h = im.size
            story.append(KeepTogether([Image(str(path),width=490,height=490*h/w),
                                      Paragraph(escaped_prose(match.group(1)),small)]))
            i += 1; continue
        if line.startswith('- '):
            flush(); story.append(Paragraph(inline_markup(line[2:],cache),body,bulletText='•'))
            i += 1; continue
        if line.startswith('<!--'):
            raise ValueError('Unresolved markup: '+line)
        paragraphs.append(line)
        i += 1
    flush()
    pages = []
    def page(canvas, doc):
        pages.append(doc.page)
        canvas.saveState()
        canvas.setFont('Helvetica',8)
        canvas.setFillColor(NAVY)
        canvas.drawString(48,A4[1]-31,'PINCER-CATMECH-AI  |  PHASE 4 RESEARCH RECORD')
        canvas.setStrokeColor(TEAL)
        canvas.line(48,A4[1]-40,A4[0]-48,A4[1]-40)
        canvas.drawString(48,29,'Model diagnostics and research plan | No new wet-lab or HPC validation')
        canvas.drawRightString(A4[0]-48,29,str(doc.page))
        canvas.restoreState()
    document = SimpleDocTemplate(str(output),pagesize=A4,leftMargin=48,rightMargin=48,
        topMargin=62,bottomMargin=50,title=lines[0].lstrip('# '),author='Pincer-CatMech-AI research record')
    document.build(story,onFirstPage=page,onLaterPages=page)
    return dict(source=source.relative_to(REPO).as_posix(),source_sha256=sha(source),
                pdf=output.name,pdf_sha256=sha(output),pages=max(pages),display_equations=formulas)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args()
    for name,file in [('Times','times.ttf'),('Times-Bold','timesbd.ttf'),('Times-Italic','timesi.ttf'),
                      ('Times-BoldItalic','timesbi.ttf'),('Song','simsun.ttc'),('Hei','simhei.ttf')]:
        pdfmetrics.registerFont(TTFont(name,str(Path('C:/Windows/Fonts')/file),subfontIndex=0))
    pdfmetrics.registerFontFamily('Times',normal='Times',bold='Times-Bold',italic='Times-Italic',boldItalic='Times-BoldItalic')
    pdfmetrics.registerFontFamily('Song',normal='Song',bold='Hei',italic='Song',boldItalic='Hei')
    cache=REPO.parent/'phase4'/'pdf-math-cache'
    cache.mkdir(parents=True,exist_ok=True)
    args.output_directory.mkdir(parents=True,exist_ok=True)
    result=[render(REPO/'docs'/f'PHASE4_RESEARCH_REPORT_{lang}.md',
                   args.output_directory/f'PHASE4_RESEARCH_REPORT_{lang}.pdf',cache) for lang in ('ZH','EN')]
    (args.output_directory/'PDF_BUILD_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
