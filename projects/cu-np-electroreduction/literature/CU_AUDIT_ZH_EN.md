# Cu 文献核验 / Cu literature audit

核验日期 / audit date: 2026-09-20。机器可读事实及逐项来源见 `cu_evidence.json`。本文件是本次审查的原创笔记；下载的机构图片、网页和文献元数据仅留作本地证据，不默认具有公开再分发许可。

## 1. 核验结论 / verified result

**LITERATURE FACT.** DOI `10.1002/anie.202505085` 的准确题名、12位作者、三位通讯作者及三种出版日期已由 [Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202505085)、[Crossref](https://api.crossref.org/works/10.1002/anie.202505085)、[PubMed](https://pubmed.ncbi.nlm.nih.gov/40107943/) 交叉核验。通讯为 Wen-Hao Li、Hai-Tao Tang、Dingsheng Wang；Ying-Ming Pan 为共同作者。accepted-online 2025-03-19、VoR-online 2025-03-27、issue 2025-05-26 并非彼此冲突。

**EXECUTED RESULT.** Publisher metadata and abstract were accessible. Main PDF and SI download attempts returned HTTP 403; full HTML redirected to the abstract. The publisher lists the SI filename `anie202505085-sup-0001-SuppMat.pdf` (8.9 MB), but its contents have **not** been read. Consequently, no SI-derived DFT protocol, Cartesian coordinates, isolated yields, charge passed or potential conversion is asserted.

**LITERATURE FACT.** [东北大学通讯作者署名报道](https://neunews.neu.edu.cn/info/1341/941741.htm)提供可公开读取的反应图。图中 1a 的 N3-甲基连接得到确认。PubChem CID 17088 对应 C9H8N2O、SMILES `CN1C=NC2=CC=CC=C2C1=O`。该连接式是从公开二维图重建，不是原文优化坐标。图中还可识别芳腈偶联伙伴、DMF、TBAPF6、TEOA、10 mA 与 C(+)|Cu-N3P1@NC(-)。天然醇先被制成4-氰基苯甲酸酯；图示前处理为 DCC/DMAP/THF。

**MODEL INFERENCE.** This information is sufficient for a molecular substrate pilot in DMF, but insufficient for reconstructing the published catalytic model. The provisional project remains scientifically relevant because stable-site adsorption, changing surface populations and interfacial transfer can yield similar gross conversion trends. Their distinction requires measurements and calculations that are absent from the accessible evidence, rather than a larger static adsorption screen.

## 2. 不应声称的内容 / claims not established

- **MODEL INFERENCE:** Source accessibility is not the same as openness or a redistribution license. A search-index OpenAlex page called the article bronze OA, whereas the fresh API contained no PDF URL or OA location. The actual failed download governs this audit.
- **MODEL INFERENCE:** The N3-methyl substrate must not be replaced silently by unsubstituted quinazolinone. Such a replacement introduces an N-H site and changes possible protonation and reduction chemistry.
- **MODEL INFERENCE:** The reaction contains an aromatic-nitrile partner. An isolated two-electron hydrogenation network would omit the central coupling and cleavage chemistry. No accepted product atom map or electron number is available yet.
- **MODEL INFERENCE:** The small source figure suggests an opened anthranilamide derivative, but this visual interpretation has not been accepted as a definitive product SMILES. High-resolution product schemes and characterization must precede net reaction balancing.
- **MODEL INFERENCE:** CuN3P1 is a coordination label, not a unique coordinate file. A hand-built graphene defect may be a proposed model, but cannot be called a published structure or a baseline reproduction.
- **REQUIRES EXPERIMENT:** Five repeated batches and reported ex-situ structural retention cannot alone establish invariant operando populations, absence of transient clusters or absence of dissolved active copper.
- **MODEL INFERENCE:** No industrial lifetime, current density, energy per mass or production cost follows from these statements. Full-cell voltage, electron balance and acceptable-product amount remain essential.

## 3. 最小计算与决策 / minimum calculations and decisions

| Calculation / 计算 | Decision / 决策 | Boundary / 边界 |
| --- | --- | --- |
| **PLANNED CALCULATION:** Optimize neutral 1a, then vibrationally check it | Is the generated molecular input a stationary minimum? | Geometry generation is not literature reproduction |
| **PLANNED CALCULATION:** Fixed-nuclear-geometry equilibrium-ALPB charging diagnostic at the neutral geometry | Is the first reduction qualitatively solvent-sensitive? | Constant-charge molecular diagnostic with equilibrium solvent response; not a physical vertical electron affinity, electrode potential or ET barrier |
| **PLANNED CALCULATION:** Relax radical anion and vibrationally check it | Does reduction materially change the C=N / C-N geometry or cause spontaneous rearrangement? | A local minimum does not establish the surface reaction pathway |
| **PLANNED CALCULATION:** Compare DMF solvation, gas-phase reference and a second affordable electronic-structure choice | Is the interpreted effect larger than model sensitivity? | Gas/solvent differences cannot be interpreted as catalytic P effects |
| **PLANNED CALCULATION, GATED:** Restore published CuN4/CuN3P1 models and selected adsorption values | Can the baseline be reproduced before new working-state screening? | Unexecuted until coordinates/settings and resource permissions are available |

**MODEL INFERENCE.** If the molecular radical anion cannot be stabilized consistently, this is a pilot diagnosis that informs model redesign. It is not a falsification of the published electrocatalytic mechanism. Likewise, a small reduction-induced bond-length change cannot validate ring-opening selectivity.

## 4. 后续文献审计 / later-work audit

**LITERATURE FACT.** A bounded exact-title/DOI search and an OpenAlex cited-by query cut off at 2026-09-20 were performed. Selected later records were checked against publisher-deposited Crossref metadata. The raw cited-by response is retained only in local audit storage as `citing_openalex.json`, indexed in `cu_source_hashes.json`; citation databases can have incomplete coverage and inconsistent dates.

**LITERATURE FACT.** [Chem review, DOI 10.1016/j.chempr.2025.102838](https://doi.org/10.1016/j.chempr.2025.102838) is a 2026 review, not an original experiment for the required nearest-neighbor count. [JACS DOI 10.1021/jacs.6c07265](https://doi.org/10.1021/jacs.6c07265), published 2026-06-16 according to Crossref, concerns a zinc single-atom cathode/iodide relay for pyridine functionalization. It is a conceptual later reference and does not directly validate Cu/quinazolinone speciation.

**MODEL INFERENCE.** No independently verified same-system follow-up emerged from this bounded search. This is a documented search outcome, not proof that no such work exists. Do not claim priority on that basis.

## 5. 最短待补证据 / shortest evidence request

**REQUIRES EXPERIMENT / LITERATURE ACCESS.** Obtain the original article and SI lawfully through existing institutional access or an author-provided copy; request CuN4/CuN3P1 coordinates and numerical adsorption outputs if not included. Extract product mapping, net balance, concentrations, time/temperature, electrode dimensions/loading, potential/reference information, FE, leaching controls and model settings. No external message has been sent.

**MODEL INFERENCE.** The immediate stop gate applies to published-catalyst reproduction and quantitative catalytic predictions. It does not prevent transparent local execution of a molecular preflight, tested data handling, conditional experimental discrimination or preparation of restartable production templates.
