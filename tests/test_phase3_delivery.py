"""Local export fixtures; never push, run git archive, or touch real outputs."""
import copy
import importlib.util
from pathlib import Path
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("phase3_delivery", ROOT / "scripts" / "export_phase3_delivery.py")
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


@pytest.fixture
def pdf_receipts(tmp_path):
    """Blank PDF/hash fixtures test receipt binding, not visual review quality."""
    PdfWriter = pytest.importorskip("pypdf", reason="PDF export checks use the declared documents extra").PdfWriter
    repo = tmp_path / "repo"
    build, qa, visual = [], [], {"documents": []}
    for name, relative in EXPORT.EXPECTED_PDFS.items():
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        count = 35 if name in EXPORT.PAPER_NAMES else 1
        writer = PdfWriter()
        for _ in range(count):
            writer.add_blank_page(width=200, height=300)
        writer.write(path)
        digest = EXPORT.sha(path)
        qa.append({"pdf": relative, "sha256": digest, "pdfium_pages": count, "pypdf_pages": count,
            "all_pages_rendered": True, "outside_page_text_count": 0,
            "pages": [{"page": i+1, "outside_page_characters": 0, "render_sha256": "0"*64} for i in range(count)]})
        visual["documents"].append({"name": name, "sha256": digest, "all_pages_visually_reviewed": True})
        if name in EXPORT.PAPER_NAMES:
            source = path.with_suffix(".md")
            source.write_text("Synthetic source fixture\n", encoding="utf-8")
            build.append({"pdf": relative, "source": str(source), "language": EXPORT.PAPER_NAMES[name],
                "source_sha256": EXPORT.sha(source), "pdf_sha256": digest, "pages": count})
    return repo, build, qa, visual


def test_all_six_pdfs_are_bound_to_actual_pages_and_final_receipts(pdf_receipts):
    verified = EXPORT.validate_pdf_evidence(*pdf_receipts)
    assert set(verified) == set(EXPORT.EXPECTED_PDFS)
    assert len(verified) == 6


def test_paper_repository_source_remains_pdf_when_figures_are_bound(pdf_receipts):
    repo, build, qa, visual = pdf_receipts
    source = repo / "docs/MAGNUM_OPUS_CATMECH_EN.md"
    figure = repo / "examples/plots/phase3_micro_solvation.png"
    figure.write_bytes(b"fixture image bytes")
    source.write_text("![Figure](../examples/plots/phase3_micro_solvation.png)\n", encoding="utf-8")
    paper = next(x for x in build if x["language"] == "EN")
    paper["source_sha256"] = EXPORT.sha(source)
    paper["figure_source_sha256"] = {"examples/plots/phase3_micro_solvation.png": EXPORT.sha(figure)}
    verified = EXPORT.validate_pdf_evidence(repo, build, qa, visual)
    assert verified["MAGNUM_OPUS_CATMECH_EN.pdf"]["repository_source"] == "docs/MAGNUM_OPUS_CATMECH_EN.pdf"


@pytest.mark.parametrize("fault", ["missing", "changed"])
def test_embedded_figure_must_match_the_pdf_build_inputs(pdf_receipts, fault):
    repo, build, qa, visual = pdf_receipts
    source = repo / "docs/MAGNUM_OPUS_CATMECH_EN.md"
    figure = repo / "examples/plots/phase3_micro_solvation.png"
    figure.write_bytes(b"fixture image bytes")
    source.write_text("![Figure](../examples/plots/phase3_micro_solvation.png)\n", encoding="utf-8")
    paper = next(x for x in build if x["language"] == "EN")
    paper["source_sha256"] = EXPORT.sha(source)
    if fault == "changed":
        paper["figure_source_sha256"] = {"examples/plots/phase3_micro_solvation.png": EXPORT.sha(figure)}
        figure.write_bytes(b"updated after paper build")
    with pytest.raises(RuntimeError, match="figure"):
        EXPORT.validate_pdf_evidence(repo, build, qa, visual)


@pytest.mark.parametrize("missing", ["build", "plot_qa", "plot_visual"])
def test_missing_required_paper_or_plot_cannot_pass(pdf_receipts, missing):
    repo, build, qa, visual = copy.deepcopy(pdf_receipts)
    if missing == "build":
        build.pop()
    elif missing == "plot_qa":
        qa.pop()
    else:
        visual["documents"].pop()
    with pytest.raises(RuntimeError, match="membership mismatch"):
        EXPORT.validate_pdf_evidence(repo, build, qa, visual)


@pytest.mark.parametrize("fault", ["duplicate", "stale_source", "stale_plot", "incomplete_pages", "unreviewed_plot"])
def test_stale_or_incomplete_pdf_evidence_is_rejected(pdf_receipts, fault):
    repo, build, qa, visual = copy.deepcopy(pdf_receipts)
    if fault == "duplicate":
        qa.append(copy.deepcopy(qa[0]))
    elif fault == "stale_source":
        (repo / "docs/MAGNUM_OPUS_CATMECH_EN.md").write_text("Changed after rendering", encoding="utf-8")
    elif fault == "stale_plot":
        qa[-1]["sha256"] = "f"*64
    elif fault == "incomplete_pages":
        qa[0]["pages"].pop()
    else:
        visual["documents"][-1]["all_pages_visually_reviewed"] = False
    with pytest.raises(RuntimeError):
        EXPORT.validate_pdf_evidence(repo, build, qa, visual)


@pytest.mark.parametrize("algorithm", ["sha1", "sha256"])
def test_git_blob_membership_handles_unicode_quotes_and_tabs(tmp_path, algorithm):
    names = ["docs/中文 paper.md", 'docs/literal"quote.txt', "docs/tab\tname.txt"]
    payloads = [b"one\n", b"two\n", b"three\n"]
    listing = "".join(f"100644 blob {EXPORT.blob_digest(data, algorithm)}\t{name}\0" for name, data in zip(names, payloads))
    tree = EXPORT.parse_git_tree(listing)
    assert set(tree) == {"pincer-catmech-ai/"+name for name in names}
    archive = tmp_path / "fixture.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        for name, data in zip(names, payloads):
            bundle.writestr("pincer-catmech-ai/"+name, data)
    EXPORT.verify_git_archive(archive, tree, algorithm)
    with zipfile.ZipFile(archive, "a") as bundle:
        bundle.writestr("unexpected.txt", b"extra")
    with pytest.raises(RuntimeError, match="membership"):
        EXPORT.verify_git_archive(archive, tree, algorithm)


def test_modified_zip_payload_fails_git_blob_verification(tmp_path):
    name = "pincer-catmech-ai/a.txt"
    archive = tmp_path / "fixture.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr(name, b"modified")
    with pytest.raises(RuntimeError, match="Git blob"):
        EXPORT.verify_git_archive(archive, {name: EXPORT.blob_digest(b"committed")})


def test_only_markdown_figure_images_are_rewritten():
    image = "![plot](../examples/plots/phase3_micro_solvation.png)"
    prose = "Path ../examples/plots/phase3_micro_solvation.png stays.\n"
    ordinary = "[source](../examples/plots/phase3_micro_solvation.png)\n"
    source = image + "\n" + prose + ordinary + "`" + image + "`\n```md\n" + image + "\n```\n"
    adapted = EXPORT.adapt_figure_image_links(source, {"phase3_micro_solvation.png"})
    assert adapted == source.replace(image, "![plot](phase3_micro_solvation.png)", 1)
    with pytest.raises(RuntimeError, match="missing"):
        EXPORT.adapt_figure_image_links(image, set())


def test_prior_release_is_rejected_and_unchanged(tmp_path):
    prior = tmp_path / "old-release"
    prior.mkdir()
    sentinel = prior / "PHASE3_DELIVERY.json"
    sentinel.write_bytes(b"prior receipt")
    with pytest.raises(FileExistsError, match="preserve"):
        EXPORT.require_fresh_output(prior, tmp_path / "repo")
    assert sentinel.read_bytes() == b"prior receipt"
    with pytest.raises(ValueError, match="outside"):
        EXPORT.require_fresh_output(tmp_path / "repo" / "new-output", tmp_path / "repo")


def test_source_snapshot_must_cover_every_current_code_file(tmp_path):
    path = tmp_path / "scripts" / "exporter.py"
    path.parent.mkdir()
    path.write_text("# tested source\n", encoding="utf-8")
    good = {"files": {"scripts/exporter.py": EXPORT.sha(path)}}
    EXPORT.validate_source_snapshot(tmp_path, good)
    (path.parent / "later.py").write_text("# new untested source\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="omits"):
        EXPORT.validate_source_snapshot(tmp_path, good)
