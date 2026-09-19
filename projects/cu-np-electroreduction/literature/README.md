# Literature audit / 文献审计

Audit cutoff: **2026-09-20**. This directory contains original research notes, bibliographic facts, source locations, access records and hashes. It does not contain third-party article bodies, SI PDFs, downloaded figures or failed-request HTML.

本目录仅提交原创整理、书目事实、来源与访问审计；第三方正文、补充材料、图片和403响应页面保留在本地审计存储，不在公共仓库中发布。

## Read in this order / 阅读顺序

1. [evidence_matrix.md](evidence_matrix.md): twelve original-study records, scientific overlap, limitations and potential contribution.
2. [evidence_matrix.json](evidence_matrix.json): the same curated records, evidence labels, relevance classes, access boundaries and author-role provenance.
3. [cu_evidence.json](cu_evidence.json) and [CU_AUDIT_ZH_EN.md](CU_AUDIT_ZH_EN.md): detailed primary Cu–N–P audit, reaction-figure extraction and missing reproduction inputs.
4. [design_references.json](design_references.json): supporting mechanistic/process references. This is not an additional disjoint nearest-neighbor count.
5. [search_log.json](search_log.json): search scope and later works, including records excluded because they are reviews or only conceptual references.

## Counts and evidence levels / 数量与证据层次

The matrix contains one primary Cu–N–P/quinazolinone study, three other organic-electroreduction comparators, four interface/speciation/method neighbors and four mandatory conceptual anchors. Eight relevant original studies are obtained only with the stated method-inclusive definition. Direct organic-electroreduction studies number four, leaving a shortfall of four under that stricter interpretation. The exact reaction has one verified original study; eight identical-reaction papers are not a requirement. Thermal Mn–Rh and aqueous CO₂ studies are never presented as direct validation of the DMF/TEOA reaction.

`LITERATURE FACT` means a source reports a claim. Publisher abstracts, original full articles, SI and an author-written institutional account have different evidentiary scope. An institutional figure supplies only visible labels and connectivity; it is not a replacement for the original SI. Crossref, PubMed and OpenAlex supply bibliography or discovery leads, not proof that a PDF was obtained. `MODEL INFERENCE`, `WORKING HYPOTHESIS` and `PLANNED CALCULATION` remain distinct from executed calculations and independent experiments.

作者信息与角色均标有来源。能核验通讯或共同贡献的才列为已核验；其他角色明确为未知。发表时机构、作者排序或通讯身份不证明当前合作、仪器、服务器或未发表数据权限。

## Access and provenance / 访问及追溯

- A0 Cu, A2 Mn and A3 Fe: publisher abstract/metadata are accessible; original body/SI remain unavailable after the recorded attempts. No inaccessible SI conditions or coordinates are asserted.
- A1 Nature chlorine anchor: public SI obtained; main article access remained limited to the publisher preview/metadata.
- A4 Mn–Rh and N5 Pd hydrogen-bond study: original full text obtained through the official EuropePMC `fullTextXML` API; SI PDF members obtained from the official `supplementaryFiles` ZIP API after publisher/PMC page failures.
- N2 dynamic Cu: original publisher HTML and SI obtained. N6 pulsed Cu and N7 2026 furfural: original publisher HTML obtained; their SI was not downloaded in this audit.
- N1 Cu/P furfural, N3 cluster theory and N4 Cu–N interface theory: primary abstract or publisher-indexed content only for the claims marked in the matrix; full numerical model details remain unaudited.

[raw_source_index.json](raw_source_index.json) maps successful local source files to source URLs, byte sizes and SHA256. For repository APIs, the ZIP member name is included. [cu_source_hashes.json](cu_source_hashes.json) indexes the separate Cu audit; unresolved exact asset/query URLs are explicitly identified rather than invented. Basenames are local audit identifiers, not files bundled in this repository. Original summaries edited during integration are excluded from raw-source hashes.

The following logs preserve outcomes without response bodies or machine-specific paths:

- [download_access_log.json](download_access_log.json)
- [supplement_access_log.json](supplement_access_log.json)
- [europepmc_access_log.json](europepmc_access_log.json)
- [europepmc_supplement_log.json](europepmc_supplement_log.json)
- [subsequent_access_log.json](subsequent_access_log.json)

A later successful official repository retrieval does not erase an earlier failed attempt. Download success also does not establish full scientific review: selected methods and relevant figures were inspected; there was no exhaustive review of all pages of every SI.

## License and local-only source storage / 许可与本地原文

Public access is not the same as permission to republish a complete work. Wiley, Nature-anchor and institutional assets had no blanket redistribution permission verified here. Accessible sources carry differing licenses: Mn–Rh CC BY-NC 3.0, dynamic/pulsed Cu CC BY 4.0, Pd study CC BY-NC 4.0, and 2026 furfural CC BY-NC-ND 4.0 according to their source notices. This package avoids redistributing any complete third-party source. URLs and hashes allow independently authorized retrieval; no credentials or paid access are included.

## Regeneration / 再生成

From this directory, run:

```console
python render_matrix.py
```

The standard-library script reads the curated JSON, checks DOI uniqueness and required evidence/role fields, and regenerates the Markdown. It works relative to its own file and contains no workstation or scratch-directory dependency. It does not download full source files or execute chemistry calculations.

The molecular pilot wording is deliberately **fixed-nuclear-geometry equilibrium-ALPB charging diagnostic**. Nuclear coordinates stay fixed while equilibrium implicit solvent responds to charge; this is not a physical vertical electron affinity, an electrode potential, or an electron-transfer barrier.
