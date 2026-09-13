"""Export a verified, already-pushed Phase 3 commit without replacing prior releases."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import uuid
import xml.etree.ElementTree as ET
import zipfile

REPO = Path(__file__).resolve().parents[1]
PAPER_NAMES = {f"MAGNUM_OPUS_CATMECH_{language}.pdf": language for language in ("EN", "ZH")}
PLOT_STEMS = ("phase3_micro_solvation", "phase3_solver_verification", "phase3_egnn_auxiliary", "phase3_spin_diagnostics")
EXPECTED_PDFS = {**{name: "docs/" + name for name in PAPER_NAMES},
    **{stem + ".pdf": "examples/plots/" + stem + ".pdf" for stem in PLOT_STEMS}}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args):
    return subprocess.check_output(["git", "-C", str(REPO), *args], text=True, encoding="utf-8").strip()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def basename(value):
    return PurePosixPath(str(value).replace("\\", "/")).name


def indexed_records(records, key, expected, label):
    if not isinstance(records, list):
        raise RuntimeError(label + " must be a list")
    indexed = {}
    for record in records:
        name = basename(record[key])
        if name in indexed:
            raise RuntimeError(f"Duplicate {label} record: {name}")
        indexed[name] = record
    if set(indexed) != set(expected):
        raise RuntimeError(f"{label} membership mismatch; missing={sorted(set(expected)-set(indexed))}, extra={sorted(set(indexed)-set(expected))}")
    return indexed


def validate_pdf_evidence(repo, build, qa, visual):
    """Bind exactly two current papers and four current plot PDFs to final QA."""
    from pypdf import PdfReader
    papers = indexed_records(build, "pdf", PAPER_NAMES, "paper build")
    qa_records = indexed_records(qa, "pdf", EXPECTED_PDFS, "PDF QA")
    reviews = indexed_records(visual["documents"], "name", EXPECTED_PDFS, "visual review")
    verified = {}
    for name, relative in EXPECTED_PDFS.items():
        path = repo / relative
        digest = sha(path)
        count = len(PdfReader(path).pages)
        record, reviewed = qa_records[name], reviews[name]
        pages = record["pages"]
        if count < 1 or record["pdfium_pages"] != count or record["pypdf_pages"] != count or len(pages) != count:
            raise RuntimeError("PDF page-count binding failed: " + name)
        if [page["page"] for page in pages] != list(range(1, count+1)):
            raise RuntimeError("PDF page coverage is incomplete: " + name)
        if digest != record["sha256"] or digest != reviewed["sha256"]:
            raise RuntimeError("PDF hash binding failed: " + name)
        if record["all_pages_rendered"] is not True or record["outside_page_text_count"] != 0 or reviewed["all_pages_visually_reviewed"] is not True:
            raise RuntimeError("PDF review is incomplete: " + name)
        if any(page["outside_page_characters"] != 0 or not re.fullmatch(r"[0-9a-fA-F]{64}", page["render_sha256"]) for page in pages):
            raise RuntimeError("PDF page render evidence is incomplete: " + name)
        if name in papers:
            paper = papers[name]
            source = repo / "docs" / (Path(name).stem + ".md")
            if paper["language"] != PAPER_NAMES[name] or basename(paper["source"]) != source.name:
                raise RuntimeError("Paper source/language identity differs: " + name)
            if sha(source) != paper["source_sha256"] or digest != paper["pdf_sha256"] or count != paper["pages"] or count < 35:
                raise RuntimeError("Paper source/PDF/page binding failed: " + name)
            image_references = re.findall(r"(?m)^!\[[^\]]*\]\(([^)]+)\)\s*$", source.read_text(encoding="utf-8"))
            expected_figures = set()
            for reference in image_references:
                figure = (source.parent / reference).resolve()
                if figure.suffix == ".svg":
                    figure = figure.with_suffix(".png")
                if not figure.is_relative_to(repo.resolve()):
                    raise RuntimeError("Paper figure lies outside repository")
                expected_figures.add(figure.relative_to(repo.resolve()).as_posix())
            bound_figures = paper.get("figure_source_sha256", {})
            if set(bound_figures) != expected_figures:
                raise RuntimeError("Paper figure input binding is incomplete: " + name)
            for figure_relative, checksum in bound_figures.items():
                if sha(repo / figure_relative) != checksum:
                    raise RuntimeError("Paper figure changed after PDF build: " + figure_relative)
        verified[name] = {"repository_source": relative, "sha256": digest, "pages": count,
            "all_pages_rendered": True, "all_pages_visually_reviewed": True}
    return verified


def validate_source_snapshot(repo, snapshot):
    files = snapshot.get("files")
    if not isinstance(files, dict) or not files:
        raise RuntimeError("Tested source snapshot must be nonempty")
    required = {path.relative_to(repo).as_posix() for folder in ("src", "scripts", "tests")
        for path in (repo / folder).rglob("*.py") if path.is_file()}
    if (repo / "pyproject.toml").exists():
        required.add("pyproject.toml")
    normalized = {str(name).replace("\\", "/"): expected for name, expected in files.items()}
    if len(normalized) != len(files):
        raise RuntimeError("Tested source snapshot has duplicate normalized paths")
    if not required.issubset(normalized):
        raise RuntimeError("Tested source snapshot omits current code: " + ", ".join(sorted(required-set(normalized))))
    for relative, expected in normalized.items():
        path = (repo / relative).resolve()
        if not path.is_relative_to(repo.resolve()) or sha(path) != expected:
            raise RuntimeError("Source changed after full tests: " + relative)


def parse_git_tree(raw, prefix="pincer-catmech-ai/"):
    """Parse -z output, retaining Unicode, spaces, tabs and literal quotes."""
    tree = {}
    for entry in raw.split("\0"):
        if not entry:
            continue
        metadata, path = entry.split("\t", 1)
        _, kind, object_id = metadata.split()
        if kind != "blob" or prefix + path in tree:
            raise RuntimeError("Unexpected non-blob or duplicate Git tree entry")
        tree[prefix + path] = object_id
    return tree


def blob_digest(payload, algorithm="sha1"):
    if algorithm not in {"sha1", "sha256"}:
        raise RuntimeError("Unsupported Git object format")
    return hashlib.new(algorithm, f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def verify_git_archive(path, tree, algorithm="sha1"):
    with zipfile.ZipFile(path) as bundle:
        members = [item for item in bundle.infolist() if not item.is_dir()]
        names = [item.filename for item in members]
        if len(names) != len(tree) or len(set(names)) != len(names) or set(names) != set(tree):
            raise RuntimeError("Export archive membership differs from committed tree")
        for member in members:
            digest = hashlib.new(algorithm, f"blob {member.file_size}\0".encode("ascii"))
            with bundle.open(member) as stream:
                for block in iter(lambda: stream.read(1024*1024), b""):
                    digest.update(block)
            if digest.hexdigest() != tree[member.filename]:
                raise RuntimeError("Archive differs from Git blob: " + member.filename)


def adapt_figure_image_links(text, available_names):
    """Adapt only inline figure-image destinations, excluding code spans/fences."""
    pattern = re.compile(r'(!\[[^\]\n]*\]\(\s*<?)\.\./examples/plots/(?P<name>phase3_[A-Za-z0-9_.-]+\.(?:png|svg|pdf))(?=[>\s)])')
    def replace(match):
        if match["name"] not in available_names:
            raise RuntimeError("Markdown image is missing from exported figures: " + match["name"])
        return match[1] + match["name"]
    output, fence = [], None
    for line in text.splitlines(keepends=True):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if marker:
            token = marker[1]
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            output.append(line)
        elif fence is not None:
            output.append(line)
        else:
            position, pieces = 0, []
            for code in re.finditer(r"(`+)(.*?)\1", line):
                pieces.extend((pattern.sub(replace, line[position:code.start()]), code[0]))
                position = code.end()
            pieces.append(pattern.sub(replace, line[position:]))
            output.append("".join(pieces))
    return "".join(output)


def require_fresh_output(output, repo):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("Output already exists; choose a new release directory to preserve previous outputs")
    if output.is_relative_to(Path(repo).resolve()):
        raise ValueError("Release output must be outside the committed repository")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    destination = require_fresh_output(args.output, REPO)
    if git("status", "--porcelain"):
        raise RuntimeError("Export requires a clean, committed checkout")
    if git("branch", "--show-current") != "main":
        raise RuntimeError("Export requires main")
    head = git("rev-parse", "HEAD")
    remote_record = git("ls-remote", "origin", "refs/heads/main").split()
    if len(remote_record) != 2 or remote_record[1] != "refs/heads/main":
        raise RuntimeError("Remote main could not be verified")
    remote = remote_record[0]
    if head != remote:
        raise RuntimeError("Local and remote main must match")
    base = REPO / "data/phase3"
    status = read(base / "PHASE3_STATUS.json")
    if not status["stages"]["spin"]["run_complete"]:
        raise RuntimeError("Spin matrix remains incomplete")
    build = read(REPO / "docs/PHASE3_PDF_BUILD_AUDIT.json")
    qa = read(base / "verification/PDF_QA.json")
    visual = read(base / "verification/PDF_VISUAL_REVIEW.json")
    source_snapshot = read(base / "verification/tested_source_manifest.json")
    validate_source_snapshot(REPO, source_snapshot)
    junit = ET.parse(base / "verification/pytest.xml").getroot()
    suites = [suite for suite in junit.iter("testsuite") if not list(suite.iter("testsuite"))[1:]]
    counts = {key: sum(int(s.get(key, "0")) for s in suites)
              for key in ("tests", "failures", "errors", "skipped")}
    if counts["tests"] <= counts["skipped"] or counts["failures"] or counts["errors"]:
        raise RuntimeError("Full pytest evidence does not pass")
    pdf_evidence = validate_pdf_evidence(REPO, build, qa, visual)
    tree = parse_git_tree(git("ls-tree", "-r", "-z", "--full-tree", head))
    algorithm = git("rev-parse", "--show-object-format")
    if algorithm not in {"sha1", "sha256"}:
        raise RuntimeError("Unsupported Git object format")
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = destination.parent / ("." + destination.name + ".staging-" + uuid.uuid4().hex)
    output.mkdir(exist_ok=False)
    selected = {}

    def copy(source, name=None):
        target = output / (name or source.name)
        if target.exists() or target.name in selected:
            raise RuntimeError("Duplicate output artifact: " + target.name)
        payload = source.read_bytes()
        relative = source.relative_to(REPO).as_posix()
        if blob_digest(payload, algorithm) != tree.get("pincer-catmech-ai/" + relative):
            raise RuntimeError("Export source differs from its committed Git blob: " + relative)
        target.write_bytes(payload)
        selected[target.name] = {"sha256": sha(target), "bytes": target.stat().st_size,
                                 "repository_source": relative, "source_git_blob": tree["pincer-catmech-ai/" + relative]}

    for language in ("EN", "ZH"):
        for extension in ("md", "pdf"):
            source = REPO / f"docs/MAGNUM_OPUS_CATMECH_{language}.{extension}"
            copy(source)
            if extension == "md":
                target = output / source.name
                available = {stem + suffix for stem in PLOT_STEMS for suffix in (".png", ".svg", ".pdf")}
                text = adapt_figure_image_links(target.read_text(encoding="utf-8"), available)
                target.write_text(text, encoding="utf-8", newline="\n")
                selected[target.name].update(sha256=sha(target), bytes=target.stat().st_size,
                                             transformation="Figure links adapted to output siblings; repository Markdown retained in ZIP")
    for stem in PLOT_STEMS:
        for suffix in (".png", ".svg", ".pdf"):
            copy(REPO / "examples/plots" / (stem + suffix))
    exports = {
        "PHASE3_STATUS.json": base / "PHASE3_STATUS.json",
        "PHASE3_PDF_BUILD_AUDIT.json": REPO / "docs/PHASE3_PDF_BUILD_AUDIT.json",
        "PHASE3_PDF_QA.json": base / "verification/PDF_QA.json",
        "PHASE3_PDF_VISUAL_REVIEW.json": base / "verification/PDF_VISUAL_REVIEW.json",
        "PHASE3_pytest.xml": base / "verification/pytest.xml",
        "PHASE3_TESTED_SOURCE_MANIFEST.json": base / "verification/tested_source_manifest.json",
        "PHASE3_SPIN_EVIDENCE_PUBLICATION.json": base / "spin/evidence_publication.json",
        "phase3_readiness.csv": base / "catalyst_temperature_base_readiness.csv",
        "phase3_solvation_association.csv": base / "solvation/association_energies.csv",
        "phase3_solvation_thermochemistry.csv": base / "solvation/association_thermochemistry.csv",
        "phase3_software_grid.csv": base / "kinetics/software_verification/grid_summary.csv",
        "phase3_software_controls.csv": base / "kinetics/software_verification/baseline_control_coefficients.csv",
        "PHASE3_PROTON_WIRE_REPRODUCIBILITY.md": REPO / "docs/PROTON_WIRE_REPRODUCIBILITY.md",
    }
    spin_value = status["stages"]["spin"]["active_matrix_summary"].replace("\\", "/")
    spin_relative = spin_value if spin_value.startswith("data/") else "data/" + spin_value.split("/data/", 1)[1]
    if not (REPO / spin_relative).resolve().is_relative_to(base.resolve()):
        raise RuntimeError("Spin matrix summary is outside Phase 3 data")
    for filename in ("spin_states.csv", "vertical_gaps.csv"):
        exports["phase3_" + filename] = (REPO / spin_relative).parent / filename
    for name, source in exports.items():
        copy(source, name)
    archive = output / "pincer-catmech-ai-phase3.zip"
    subprocess.run(["git", "-C", str(REPO), "archive", "--format=zip", "--prefix=pincer-catmech-ai/",
                    "--output=" + str(archive), head], check=True)
    # Bind every exported member directly to its committed Git blob. This also
    # reads each member fully, invoking ZIP CRC validation without extraction.
    verify_git_archive(archive, tree, algorithm)
    selected[archive.name] = {"sha256": sha(archive), "bytes": archive.stat().st_size,
                              "git_commit": head, "git_blob_members_verified": len(tree)}
    readme = output / "PHASE3_README_ZH.md"
    readme.write_text(
        "# Phase 3 交付说明\n\n"
        f"已提交并核对远程 main：`{head}`。完整仓库与原生计算证据位于 `pincer-catmech-ai-phase3.zip`。\n\n"
        f"完整回归测试：{counts['tests']} 项，失败 {counts['failures']}，错误 {counts['errors']}。"
        "两份白皮书及四组图已按最终文件哈希核验。\n\n"
        "真实结果包括微溶剂簇优化与频率、自旋态计算尝试、脱水路径尝试和辅助构象能差 EGNN 训练。"
        "自旋质量、MECP、TS 与完整势垒网络分别验收，缺失结果没有填零。"
        "动力学网格与控制系数文件名称含 software，使用任意速率检验算法，不是催化性能预测。"
        "EGNN 未优于零能差基线，势垒与自旋能隙输出头保持未训练。\n\n"
        "原生 ZIP、完整 CSV、失败记录、脚本和测试均在完整仓库包中。"
        "Psi4 的大型 DF-SCF 积分与 DIIS 中间缓存按清单保留在原计算机，未放入 Git/ZIP；"
        "其路径、大小与 SHA256 在 PHASE3_SPIN_EVIDENCE_PUBLICATION.json 中逐项记录。"
        "运行前复用既有科学环境；Windows xTB/Psi4 的 DLL 路径应按项目说明设置。"
        "详细计算条件和数量以 `PHASE3_STATUS.json` 为准，文件完整性以 `PHASE3_DELIVERY.json` 为准。\n",
        encoding="utf-8", newline="\n")
    selected[readme.name] = {"sha256": sha(readme), "bytes": readme.stat().st_size}
    receipt = {"schema": "phase3_delivery_v1", "repository": git("remote", "get-url", "origin"),
        "local_commit": head, "remote_main_verified": remote, "working_tree_clean": True,
        "pytest": counts, "paper_pages": {x["language"]: x["pages"] for x in build},
        "all_archive_members_match_committed_git_blobs": True, "git_object_format": algorithm,
        "verified_pdfs": pdf_evidence, "files": selected,
        "scientific_status": status["scientific_status"],
        "physical_18_step_kinetic_predictions": status["accepted_physical_18_step_kinetic_predictions"],
        "accepted_proton_wire_TS": status["stages"]["proton_wire"]["accepted_TS_count"],
        "spin_cache_publication_policy": read(base / "spin/evidence_publication.json")["policy"],
        "prior_outputs_preserved": True}
    (output / "PHASE3_DELIVERY.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if git("status", "--porcelain") or git("rev-parse", "HEAD") != head or git("ls-remote", "origin", "refs/heads/main").split() != [head, "refs/heads/main"]:
        raise RuntimeError("Repository or remote changed during export; staged output was not published")
    require_fresh_output(destination, REPO)
    output.rename(destination)
    print(json.dumps({"commit": head, "archive_bytes": selected[archive.name]["bytes"],
                      "archive_members": len(tree), "pytest": counts, "output": str(destination)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
