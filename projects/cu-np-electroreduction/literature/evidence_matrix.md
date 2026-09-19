# Literature evidence and novelty audit / 文献证据与创新边界审计

Audit date / 检索截止：2026-09-20。

`LITERATURE FACT` records what a source reports; it does not mean independent reproduction. `MODEL INFERENCE` denotes this audit’s interpretation. `WORKING HYPOTHESIS` denotes an unproven contribution. / 文献事实不等于独立复现；推断与待检验假说均单独标示。

## Count and relevance / 数量与相关性

| Class | Count | Entries |
|---|---:|---|
| Primary Cu–N–P/quinazolinone system | 1 | A0 |
| Other organic-electroreduction comparators | 3 | N1, N5, N7 |
| Interface/speciation/method neighbors | 4 | N2, N3, N4, N6 |
| Mandatory conceptual anchors | 4 | A1–A4 |

The matrix contains 12 original studies. Eight are relevant under an explicit primary-plus-organic-plus-methodological definition; four conceptual anchors are excluded from that count. Only four directly study organic electroreduction, a shortfall of four if the eight-study threshold is interpreted to require that narrower scope. One study concerns the exact primary reaction; no additional same-system original study was verified. The contract does not require eight studies on an identical reaction.

共12项原创研究；仅在明确包含界面、位点和方法近邻的定义下，核心及近邻合计8项。直接研究有机电还原的只有4项；若将“8项直接相关”限定于有机电还原，则尚缺4项。同一Cu–N–P/喹唑啉酮体系仅1项，没有将热催化或CO₂研究称作同反应研究。

## Evidence matrix / 证据矩阵

| ID / source | System/electrolyte: LITERATURE FACT | Claim/evidence: LITERATURE FACT | Interface: LITERATURE FACT | Limits/overlap: MODEL INFERENCE | Possible contribution: WORKING HYPOTHESIS |
|---|---|---|---|---|---|
| A0 [10.1002/anie.202505085](https://doi.org/10.1002/anie.202505085) | Cu-N-P@NC, assigned CuN3P1; 3-methylquinazolin-4(3H)-one and aromatic nitrile coupling/ring opening. Author-institution figure labels TBAPF6, TEOA, DMF, 10 mA, carbon anode and Cu-N3P1@NC cathode. Concentrations and cathode potential unknown. | Authors attribute P promotion to metal-support electronic effects and substrate adsorption/desorption; abstract reports gram-scale chemistry and 11 natural-product modifications. Publisher abstract/metadata plus author-written NEU account and inspected reaction figure; no raw computational or experimental reanalysis. | Exact coordinates, charge/spin, functional, solvent and potential treatment are unavailable because the original body/SI were not obtained. | Cannot yet reproduce catalyst energetics or balance the net electron/proton reaction. A low-resolution product depiction is insufficient for an accepted atom map. Repeating CuN3P1 versus CuN4 static adsorption logic would mainly reproduce the published rationale. | Discriminate coordination-only adsorption from potential-dependent populations and donor/solvent competition, then test a frozen prospective condition set. |
| A1 [10.1038/s41586-023-05886-z](https://doi.org/10.1038/s41586-023-05886-z) | Amide organocatalyst; chloride oxidation; SI electrochemical tests use 5 M NaCl, pH 2. Not a Cu SAC or organic electroreduction. | Authors attribute improved chlorine evolution to reversible CO2 binding and radical generation. Abstract reports 10 kA m^-2, 99.6% selectivity, 89 mV overpotential. Electrochemistry, operando ATR-SEIRAS, molecular characterization and calculated redox thermodynamics. | SI: Gaussian09; B3LYP/6-31G(d) geometries; water IEF-PCM; B3LYP-D3(BJ)/def2-TZVPP single points plus thermal corrections. Redox potentials from thermodynamic cycles; not a constant-potential electrode model. | Aqueous CER performance cannot establish nonaqueous Cu-organic process performance. SI S69 notes calculated water/pH7 conditions versus pH2 experiments. Reproducing CO2-enabled redox-mediator logic is conceptual reproduction. | Test whether potential-dependent catalyst populations, rather than a single molecular redox shift, control the target Cu reaction. |
| A2 [10.1002/anie.202315032](https://doi.org/10.1002/anie.202315032) | Mn-SA@NC; silane-to-silanol oxidation. Electrolyte, solvent, electrode and applied-current conditions unavailable in inspected abstract. | Authors propose an OER-derived MOOH intermediate mediates silane oxidation; abstract reports 600 ppm loading, TON 9132 and a 10 mmol flow demonstration. Publisher abstract describes synthesis, catalytic scope and flow demonstration; mechanistic claim remains author-reported at abstract evidence level. | Detailed computational charge, solvent, potential and ion treatments not accessible; no assumptions filled in. | Cannot independently assess mechanistic controls or model reproduction without original body/SI. Transferring the OER-intermediate concept alone supplies no new evidence for reductive Cu chemistry. | Use adsorbate-state competition as a question; do not transplant an oxidizing MOOH pathway to the Cu cathode. |
| A3 [10.1002/anie.202404295](https://doi.org/10.1002/anie.202404295) | Fe-SA@NC; organic anodic oxidation, including cyclopropyl-amide ring opening. Exact electrolyte and electrolysis conditions not accessible. | Authors describe heterogeneous redox mediation, C-C/C-X construction, complex-molecule tolerance, gram scale and flow electrochemistry. Publisher abstract verifies the reported reaction family and scale claims; individual mechanistic experiments not inspected. | No accessible detailed charge, potential-reference, solvation or ion model; not reconstructed from abstract. | Abstract cannot separate electron mediation from adsorption, support or leached-species effects. A heterogeneous mediator for ring opening is already established conceptually; the target must add a discriminating explanation. | Separate static coordination, changing site populations and interfacial transfer using the target reductive system. |
| A4 [10.1039/D4SC08658A](https://doi.org/10.1039/D4SC08658A) | Mn1-Rh1@O-TiC; cumene to acetophenone. SI: 40 mg catalyst, 2.9 mmol cumene, O2 balloon, 120 C, 400 rpm, 12-18 h; no added reaction solvent specified. Thermal, not electrochemical. | Authors propose Mn-promoted oxidation and Rh-assisted dehydration/relay oxidation; main text reports 99% conversion and 94% isolated yield. XAS/STEM, intermediates, EPR, hot filtration, ICP, cycling and calculated intermediate free energies. | SI: spin-polarized VASP/PAW/PBE-D3, 450 eV, 15 A vacuum, 1e-5 eV and 0.05 eV/A tolerances; thermal corrections at 298.15 K. No electrode-potential model. | Methods inspected do not establish validated transition states. Intermediate-energy differences must not become kinetic barriers. Experimental 393 K versus thermochemical corrections at 298 K requires care. Relay-catalysis explanation is not novel, and this reaction cannot validate cathodic Cu kinetics. | Apply falsification and explicit competing-site controls without adding Mn/Rh to the core project. |
| N1 [10.1039/D0GC03999C](https://doi.org/10.1039/D0GC03999C) | Single-atom Cu with phosphorus dopants; furfural reduction at aqueous pH5. Electrolyte salt not verified from accessible abstract. | Authors associate P with sequential reduction toward methylfuran and atomically dispersed Cu with HER suppression; FE >90% furfuryl alcohol at -0.75 V and 60% methylfuran at -0.90 V vs RHE. Primary publisher abstract, authors, DOI, pages and date verified. | Potential reference and bulk pH explicit; computational interface treatment unavailable in inspected material. | Phosphorus doping does not by itself prove the same CuN3P1 local motif. Furfural and aqueous media differ from the target. P-enhanced Cu-SAC organic electroreduction is prior art; introduction of P alone is not novelty. | Establish a target-specific potential/speciation crossover with prospective predictions rather than repeat P promotion. |
| N2 [10.1038/s41467-023-40970-y](https://doi.org/10.1038/s41467-023-40970-y) | N-coordinated single-atom Cu; CO2 reduction in 0.1 M KHCO3; reported potentials vs RHE. | Potential-driven Cu-N weakening creates low-coordination states; Cu-N changes partly reverse, whereas Cu-Cu formation is irreversible for this material. Operando/in situ TEM, XAS, spectroscopy, product distributions and electronic-structure calculations. | Methods: Gaussian09 B3LYP-D3(BJ); H/C/N/O 6-31G(d,p), Cu LANL2DZ; no constant-potential electrode ensemble established in inspected computational details. | Does not show that every Cu-N-C support, solvent or adsorbate reconstructs identically; ex situ recovery is not a universal state descriptor. Claiming that Cu coordination changes under bias would reproduce known behavior. | Determine whether P shifts the speciation crossover in the actual organic electrolyte and affects product release. |
| N3 [10.1021/jacs.4c05669](https://doi.org/10.1021/jacs.4c05669) | Cu sites on carbon nitride, CO2/CO electroreduction; solvent and ionic composition not independently extracted. | Authors calculate that isolated Cu is unstable and cannot accomplish the considered multicarbon pathways; leaching and clustering create alternative activity. Verified DOI and primary abstract; this is computational evidence, not independent experimental confirmation. | Publisher-indexed text describes constant-potential/hybrid-solvent calculations. Detailed algorithm and settings remain unverified because full text/SI access failed. | Requirement for multiple Cu centers in CO-CO chemistry cannot be transferred to quinazolinone ring opening. Generic statements that Cu clusters may be active are prior art. | Test bounded Cu-cluster/leached-species controls for this organic substrate without assuming their dominance. |
| N4 [10.1021/acs.nanolett.5c01245](https://doi.org/10.1021/acs.nanolett.5c01245) | Cu-N-C at an electrochemical interface; proton transfer and metal-N stability, not a quinazolinone reaction. | Authors link potential-sensitive Cu-N antibonding occupancy and proton-induced electronic reordering to metal leaching. Original computational study; title, authors and journal bibliographic record independently verified. | Hybrid-solvation constant-potential simulations claimed in primary abstract. SI contents list E-U fits and charge-potential tables; numerical settings not inspected. | Potential fit quality, sampled solvent configurations and transfer to nonaqueous organic electrolyte remain unaudited. Combining Cu-N speciation, proton transfer and electrode potential is already a published methodological direction. | Measure the decision value of applying a calibrated version to Cu-N-P and substrate/product occupancy. |
| N5 [10.1126/sciadv.adu1602](https://doi.org/10.1126/sciadv.adu1602) | Pd1-CuOx on Cu foam; chloramphenicol hydrodechlorination in 0.2 M Na2SO4, -1.2 V vs Hg/HgO; carbonyl/nitrile tests also reported. | Authors connect water organization and proton transfer to hydrogenation. CAP dechlorination ratio 99% is distinct from its 17.5% hydrodechlorination FE. Operando Raman/FTIR, EPR, H/D KIE, hydrogen-bond perturbation and AIMD. | VASP PBE-D3/PAW 400 eV; 68 explicit waters; 300 K AIMD and slow-growth constraints. Methods describe adding 1 electron via NELECT and mapping potentials. | A fixed added electron and potential mapping do not alone verify constant-potential feedback; authors use that label. Aqueous Pd-water behavior cannot establish Cu-N-P nonaqueous behavior. Hydrogen-bond control of organic ECH is prior art. | Discriminate donor-mediated transfer from active-site populations using matched target-specific donor/activity experiments. |
| N6 [10.1038/s41467-024-50379-w](https://doi.org/10.1038/s41467-024-50379-w) | Cu-N-C; CO2 reduction in CO2-saturated 0.1 M KHCO3; alternating cathodic/anodic pulses. | Operando tracking shows reversible interconversion of single sites, clusters and nanoparticles with pulse-dependent product distributions. Operando quick XAS, electrolysis and a semiquantitative kinetic model. | Experimentally controlled potentials; Ag/AgCl calibrated against RHE. Four-process population model, not an atomistic constant-potential DFT calculation. | Local pH, adsorbates and capacitive-current corrections complicate pulse interpretation; reversibility differs from N2 because material/protocol differ. Pulse-controlled Cu ensembles and kinetic population models are prior art. | Use hysteresis and post-pulse persistence to distinguish target catalyst changes from instantaneous electron-transfer effects. |
| N7 [10.1038/s41467-026-74394-1](https://doi.org/10.1038/s41467-026-74394-1) | Mixed-valence Cu clusters/nanoparticles in N-doped carbon; furfural hydrogenolysis. Standard tests: 25 mM substrate, 30 mL 0.5 M H2SO4 with 20 vol% MeCN (reported pH1.0), -0.65 V vs RHE. | Authors associate Cu valence and intermediate coverage with PCET-selective methylfuran formation versus HAT-related alcohol formation. Operando XAS/FTIR/Raman, product-selective KIE, adsorption controls and DFT. Published 6 July 2026; version of record 14 August 2026. | VASP PBE-D3/PAW 450 eV, 20 A vacuum, 1e-5 eV and 0.05 eV/A tolerances. Inspected methods do not specify a grand-canonical electrode or explicit-solvent ensemble. | Cluster morphology, aqueous acid and furfural differ from the target; DFT alone cannot establish the interfacial kinetic mechanism. Cu-valence/coverage and PCET-versus-HAT arguments in organic reduction are already represented in subsequent work. | Determine whether the intact Cu-N-P model remains adequate against experimentally bounded cluster alternatives. |

## Source access, locators and author roles / 访问范围、证据定位与作者角色

### A0: Phosphorus-Doped Single Atom Copper Catalyst as a Redox Mediator in the Cathodic Reduction of Quinazolinones

- Classification: `primary_same_system`.
- Article: NOT OBTAINED: full HTML redirects to abstract; publisher PDF request returned HTTP 403; publisher abstract and author-institution account inspected
- Supporting information: NOT OBTAINED: publisher exposes file name and size, but download returned HTTP 403
- Evidence locations: Publisher Abstract, author information and publication history; NEU author account dated 2025-04-01; reaction figure panel a; Detailed extraction and missing fields: cu_evidence.json
- Sources: [1](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202505085); [2](https://neunews.neu.edu.cn/info/1341/941741.htm); [3](https://neunews.neu.edu.cn/__local/B/1C/67/245ABE2C221B0DDDFB1A80F654F_DB14C13F_2DF5D.png); [4](https://api.crossref.org/works/10.1002/anie.202505085); [5](https://pubmed.ncbi.nlm.nih.gov/40107943/)
- Authors in order: Xin-Yu Wang; Wan-Jie Wei; Si-Yu Zhou; Yong-Zhou Pan; Jiarui Yang; Tao Gan; Zechao Zhuang; Wen-Hao Li; Xia Zhang; Ying-Ming Pan; Hai-Tao Tang; Dingsheng Wang.
- Corresponding authors: Wen-Hao Li; Hai-Tao Tang; Dingsheng Wang.
- Role verification: Published GXNU affiliation and corresponding status of Hai-Tao Tang verified. Individual computation/experiment contributions unavailable in accessible article view. No current collaboration or access inferred.
- Role source: [original source](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202505085). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### A1: CO2-mediated organocatalytic chlorine evolution under industrial conditions

- Classification: `mandatory_conceptual_anchor`.
- Article: publisher abstract, author information and data statement read; full article unavailable
- Supporting information: downloaded public publisher PDF; selected methods pages inspected
- Evidence locations: publisher Abstract and Author information; SI S8 computational method; S16-S17 electrochemistry; S69 thermodynamic interpretation
- Sources: [1](https://www.nature.com/articles/s41586-023-05886-z); [2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-023-05886-z/MediaObjects/41586_2023_5886_MOESM1_ESM.pdf)
- Authors in order: Jiarui Yang; Wen-Hao Li; Hai-Tao Tang; Ying-Ming Pan; Dingsheng Wang; Yadong Li.
- Corresponding authors: Dingsheng Wang; Yadong Li.
- Role verification: Hai-Tao Tang and Ying-Ming Pan affiliated with GXNU; catalyst synthesis contributions. Jiarui Yang performed theoretical calculations.
- Role source: [original source](https://www.nature.com/articles/s41586-023-05886-z). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.
- Equal contribution explicitly reported: Jiarui Yang; Wen-Hao Li; Hai-Tao Tang.

### A2: Single-Atom Manganese-Catalyzed Oxygen Evolution Drives the Electrochemical Oxidation of Silane to Silanol

- Classification: `mandatory_conceptual_anchor`.
- Article: publisher abstract and author roles accessible; full-text link redirects to abstract
- Supporting information: publisher listed filename verified; download HTTP403
- Evidence locations: publisher Abstract, Author information, Supporting Information filename
- Sources: [1](https://onlinelibrary.wiley.com/doi/10.1002/anie.202315032)
- Authors in order: Hai-Tao Tang; He-Yang Zhou; Ying-Ming Pan; Jia-Lan Zhang; Fei-Hu Cui; Wen-Hao Li; Dingsheng Wang.
- Corresponding authors: Hai-Tao Tang; Wen-Hao Li; Dingsheng Wang.
- Role verification: First five authors affiliated with GXNU. Detailed individual research contributions unavailable in abstract view.
- Role source: [original source](https://onlinelibrary.wiley.com/doi/10.1002/anie.202315032). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### A3: Single-Atom Iron Catalyst as an Advanced Redox Mediator for Anodic Oxidation of Organic Electrosynthesis

- Classification: `mandatory_conceptual_anchor`.
- Article: publisher abstract and author roles read; full-text link redirects to abstract
- Supporting information: filename verified; HTTP403 download
- Evidence locations: publisher Abstract and Author information
- Sources: [1](https://onlinelibrary.wiley.com/doi/10.1002/anie.202404295); [2](https://pubmed.ncbi.nlm.nih.gov/38649323/)
- Authors in order: Xin-Yu Wang; Yong-Zhou Pan; Jiarui Yang; Wen-Hao Li; Tao Gan; Ying-Ming Pan; Hai-Tao Tang; Dingsheng Wang.
- Corresponding authors: Wen-Hao Li; Hai-Tao Tang; Dingsheng Wang.
- Role verification: GXNU affiliations for Xin-Yu Wang, Yong-Zhou Pan, Ying-Ming Pan and Hai-Tao Tang. Individual computation/experiment contributions unavailable in abstract view.
- Role source: [original source](https://onlinelibrary.wiley.com/doi/10.1002/anie.202404295). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.
- Equal contribution explicitly reported: Xin-Yu Wang; Yong-Zhou Pan.

### A4: A Mn-Rh dual single-atom catalyst for inducing C-C cleavage: relay catalysis reversing chemoselectivity in C-H oxidation

- Classification: `mandatory_conceptual_anchor`.
- Article: complete official publisher HTML/EuropePMC XML read
- Supporting information: complete 71-page PDF downloaded via official EuropePMC API; relevant methods inspected
- Evidence locations: main Fig5-Fig6 and Author contributions; SI pp4-5 computational method; p8 general procedure
- Sources: [1](https://pubs.rsc.org/en/content/articlehtml/2025/sc/d4sc08658a?page=search); [2](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11934264/fullTextXML); [3](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11934264/supplementaryFiles)
- Authors in order: Chang-Jie Yang; Yu-Da Huang; Yu-Yuan Zhang; Yong-Zhou Pan; Jiarui Yang; Ying-Ming Pan; Tao Gan; Hai-Tao Tang; Xia Zhang; Wen-Hao Li; Dingsheng Wang.
- Corresponding authors: Hai-Tao Tang; Wen-Hao Li; Dingsheng Wang.
- Role verification: GXNU affiliations for Chang-Jie Yang, Yu-Da Huang, Yu-Yuan Zhang, Ying-Ming Pan and Hai-Tao Tang. Ying-Ming Pan performed theoretical calculation; Li/Tang/Wang conceived, designed, planned synthesis and acquired funding.
- Role source: [original source](https://pubs.rsc.org/en/content/articlehtml/2025/sc/d4sc08658a?page=search). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N1: Selective electrochemical hydrogenation of furfural to 2-methylfuran over a single atom Cu catalyst under mild pH conditions

- Classification: `direct_organic_electroreduction_comparator`.
- Article: publisher abstract available; full-text fetch HTTP403
- Supporting information: listed publicly on publisher; content not obtained
- Evidence locations: publisher Abstract and Article information
- Sources: [1](https://pubs.rsc.org/en/content/articlelanding/2021/gc/d0gc03999c/unauth)
- Authors in order: Peng Zhou; Yu Chen; Peng Luan; Xiaolong Zhang; Ziliang Yuan; Si-Xuan Guo; Qinfen Gu; Bernt Johannessen; Mamun Mollah; Alan L. Chaffee; David R. Turner; Jie Zhang.
- Corresponding authors: Qinfen Gu; Jie Zhang.
- Role verification: Corresponding status verified on publisher landing page; individual experimental/computational roles not inspected.
- Role source: [original source](https://pubs.rsc.org/en/content/articlelanding/2021/gc/d0gc03999c/unauth). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N2: Activating dynamic atomic-configuration for single-site electrocatalyst in electrochemical CO2 reduction

- Classification: `methodological_or_speciation_neighbor`.
- Article: complete publisher HTML downloaded and methods read
- Supporting information: 98-page public publisher PDF downloaded; selected pages inspected
- Evidence locations: main Results and Computational details
- Sources: [1](https://www.nature.com/articles/s41467-023-40970-y); [2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-40970-y/MediaObjects/41467_2023_40970_MOESM1_ESM.pdf)
- Authors in order: Chia-Shuo Hsu; Jiali Wang; You-Chiuan Chu; Jui-Hsien Chen; Chia-Ying Chien; Kuo-Hsin Lin; Li Duan Tsai; Hsiao-Chien Chen; Yen-Fa Liao; Nozomu Hiraoka; Yuan-Chung Cheng; Hao Ming Chen.
- Corresponding authors: Yuan-Chung Cheng; Hao Ming Chen.
- Role verification: Publisher Contributions identifies Jui-Hsien Chen and Yuan-Chung Cheng with computational investigation; Hao Ming Chen directed the project.
- Role source: [original source](https://www.nature.com/articles/s41467-023-40970-y). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N3: Not One, Not Two, But at Least Three: Activity Origin of Copper Single-Atom Catalysts toward CO2/CO Electroreduction to C2+ Products

- Classification: `methodological_or_speciation_neighbor`.
- Article: publisher-indexed abstract accessible; direct publisher fetch HTTP403
- Supporting information: not retrieved
- Evidence locations: publisher abstract; complete model settings not verified
- Sources: [1](https://pubs.acs.org/doi/abs/10.1021/jacs.4c05669); [2](https://pubmed.ncbi.nlm.nih.gov/38804682/)
- Authors in order: Juan Zhang; Yu Wang; Yafei Li.
- Corresponding authors: NOT VERIFIED in accessible source.
- Role verification: Author list verified bibliographically; corresponding and individual contribution roles not verified.
- Role source: [original source](https://pubs.acs.org/doi/abs/10.1021/jacs.4c05669). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N4: Dynamic Structural Evolution of Single-Atom Catalysts at the Catalyst-Electrolyte Interface: Insights from Electrochemical Coupled Field

- Classification: `methodological_or_speciation_neighbor`.
- Article: primary abstract and ACS issue table read; direct article HTTP403
- Supporting information: publisher description available, content not retrieved
- Evidence locations: abstract and Supporting Information description only
- Sources: [1](https://pubmed.ncbi.nlm.nih.gov/40190161/); [2](https://pubs.acs.org/doi/10.1021/acs.nanolett.5c01245); [3](https://pubs.acs.org/nalefd/issue/25/15)
- Authors in order: Xiaotao Zhang; Jiao Chen; Hongyan Wang; Yongliang Tang; Yuan Ping Feng; Yuanzheng Chen; Zhongfang Chen.
- Corresponding authors: NOT VERIFIED in accessible source.
- Role verification: Author list verified in ACS issue table; corresponding and individual contributions not verified.
- Role source: [original source](https://pubs.acs.org/nalefd/issue/25/15). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N5: Pd single atoms guided proton transfer along an interfacial hydrogen bond network for efficient electrochemical hydrogenation

- Classification: `direct_organic_electroreduction_comparator`.
- Article: complete original article read as official EuropePMC XML
- Supporting information: public PDF downloaded via EuropePMC API; not exhaustively audited
- Evidence locations: Results Electrochemical hydrogenation performance; Materials and Methods Electrochemical measurements and Dynamics simulations
- Sources: [1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12333693/); [2](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12333693/fullTextXML); [3](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12333693/supplementaryFiles)
- Authors in order: Rui Zhao; Qi Wang; Yancai Yao; Ruizhao Wang; Long Zhao; Zhiwei Hu; Cheng-Wei Kao; Ting-Shan Chan; Wenhuai Li; Qian Zheng; Jiaxian Wang; Xingyue Zou; Kaiyuan Wang; Jie Dai; Xiang-Kui Gu; Lizhi Zhang.
- Corresponding authors: Yancai Yao; Jie Dai; Xiang-Kui Gu; Lizhi Zhang.
- Role verification: Author roles from official EuropePMC fullTextXML author metadata; no inferred roles.
- Role source: [original source](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12333693/fullTextXML). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.
- Equal contribution explicitly reported: Rui Zhao; Qi Wang.

### N6: Reversible metal cluster formation on Nitrogen-doped carbon controlling electrocatalyst particle size with subnanometer accuracy

- Classification: `methodological_or_speciation_neighbor`.
- Article: complete publisher HTML read
- Supporting information: not downloaded in this audit
- Evidence locations: main Size dependency, Catalytic selectivity, Discussion and Methods
- Sources: [1](https://www.nature.com/articles/s41467-024-50379-w); [2](https://pubmed.ncbi.nlm.nih.gov/39030207/)
- Authors in order: Janis Timoshenko; Clara Rettenmaier; Dorottya Hursan; Martina Ruscher; Eduardo Ortega; Antonia Herzog; Timon Wagner; Arno Bergmann; Uta Hejral; Aram Yoon; Andrea Martini; Eric Liberra; Mariana Cecilio de Oliveira Monteiro; Beatriz Roldan Cuenya.
- Corresponding authors: Janis Timoshenko; Beatriz Roldan Cuenya.
- Role verification: Publisher Contributions: Timoshenko/Cuenya designed study; Timoshenko led XAS analysis, Rettenmaier catalytic tests, Hursan synthesis. Names with diacritics are normalized here.
- Role source: [original source](https://www.nature.com/articles/s41467-024-50379-w). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

### N7: Breaking the activity-selectivity dilemma in direct electrocatalytic hydrogenolysis of furfural into furanic biofuel

- Classification: `direct_organic_electroreduction_comparator`.
- Article: complete official publisher HTML downloaded and methods read
- Supporting information: not downloaded in this audit
- Evidence locations: main Fig3 conditions; Fig4 KIE; Theoretical calculations; publication dates
- Sources: [1](https://www.nature.com/articles/s41467-026-74394-1)
- Authors in order: Xiuzheng Zhuang; Ke Liu; Wenjie Zhang; Gehao Chen; Qi Zhang; Lungang Chen; Xinghua Zhang; Longlong Ma.
- Corresponding authors: Xinghua Zhang; Longlong Ma.
- Role verification: Publisher Contributions: Longlong Ma designed/supervised; Xiuzheng Zhuang conducted experimental work and drafted manuscript.
- Role source: [original source](https://www.nature.com/articles/s41467-026-74394-1). Affiliation and author roles describe the publication; no present collaboration, facility access or available data are inferred.

## Novelty boundary / 创新边界

1. **MODEL INFERENCE:** N1 establishes a Cu/P organic-electroreduction precedent. Adding phosphorus or redrawing a static adsorption diagram is insufficient novelty.
2. **MODEL INFERENCE:** N2–N4 and N6 already address dynamic Cu sites, proton-associated instability or clusters. Their existence motivates controls but does not prove that these states occur in DMF/TEOA.
3. **MODEL INFERENCE:** N2 and N6 report different reconstruction reversibility. Support, electrolyte and potential protocol must be examined before transferring either result.
4. **MODEL INFERENCE:** N5’s dechlorination ratio and FE are different observables. Its fixed added-charge description does not independently establish constant-potential feedback.
5. **WORKING HYPOTHESIS:** A useful target-specific contribution would distinguish coordination-only adsorption from changing site populations or interfacial transfer and predict a condition where the Cu–N–P advantage weakens. No validated catalyst ranking is supplied by this audit.

Reproduction, application to a new system, mechanistic discovery, predictive methodology and prospective experimental validation are separate achievements. DFT labels do not supply independent experimental validation.

## Retrieval and search / 检索与原文追溯

See [search_log.json](search_log.json) for query scope and the bounded follow-up search. The 2026 original N7 predates the audit cutoff; it is not a same-system follow-up. Reviews and unrelated citing works are not counted. Absence from a bounded search is not proof of absence.

See [README.md](README.md) for access levels, source manifests and local-only third-party files. No inaccessible SI parameters have been reconstructed from abstracts.
