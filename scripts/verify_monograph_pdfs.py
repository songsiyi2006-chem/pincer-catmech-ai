"""Render every PDF page and record text bounds, page counts and image hashes."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json

from PIL import Image, ImageOps, ImageDraw
import pypdfium2 as pdfium
from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reports = []
    for pdf in args.directory.glob("MONOGRAPH*.pdf"):
        document = pdfium.PdfDocument(pdf)
        reader = PdfReader(pdf)
        folder = args.output / pdf.stem
        folder.mkdir(exist_ok=True)
        pages, outside = [], []
        for index in range(len(document)):
            page = document[index]
            bitmap = page.render(scale=1.25)
            image = bitmap.to_pil()
            path = folder / f"page_{index+1:03d}.png"
            image.save(path)
            text_page = page.get_textpage()
            text = text_page.get_text_range()
            width, height = page.get_size()
            violations = []
            for character in range(text_page.count_chars()):
                left, bottom, right, top = text_page.get_charbox(character)
                if left < -1 or bottom < -1 or right > width+1 or top > height+1:
                    violations.append(character)
            outside.extend([[index + 1, item] for item in violations])
            pages.append({"page": index + 1, "characters": len(text), "outside_page_characters": len(violations),
                          "render": str(path), "render_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            text_page.close(); page.close(); bitmap.close()
        for start in range(0, len(pages), 6):
            sheet = Image.new("RGB", (1490, 3230), "#dbe5ec")
            draw = ImageDraw.Draw(sheet)
            for cell, item in enumerate(pages[start:start+6]):
                image = Image.open(item["render"]).convert("RGB")
                image.thumbnail((710, 1030))
                x, y = (cell % 2) * 745 + 17, (cell // 2) * 1075 + 28
                sheet.paste(image, (x, y))
                draw.text((x, y-18), f"PAGE {item['page']}", fill="black")
            sheet.save(folder / f"contact_{start//6+1:02d}.png")
        report = {"pdf": str(pdf), "sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
                  "render_runtime": {name: importlib.metadata.version(name) for name in ("pypdf", "pypdfium2", "Pillow")},
                  "pdfium_pages": len(document), "pypdf_pages": len(reader.pages),
                  "all_pages_rendered": len(pages) == len(reader.pages),
                  "outside_page_text_count": len(outside), "pages": pages}
        (folder / "extracted_text.txt").write_text("\n\f\n".join(p.extract_text() for p in reader.pages), encoding="utf-8")
        document.close()
        reports.append(report)
    (args.output / "PDF_QA.json").write_text(json.dumps(reports, indent=2), encoding="utf-8")
    print(json.dumps([{k:v for k,v in r.items() if k != "pages"} for r in reports], indent=2))


if __name__ == "__main__":
    main()
