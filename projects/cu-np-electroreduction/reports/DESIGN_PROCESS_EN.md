# Cu–N–P organic electroreduction: mechanism discrimination and process assessment

**Status note:** This supporting design contains planned branches. The executed-status report `TECHNICAL_REPORT_EN.md` and frozen `../config/prospective_holdout.json` are authoritative for the delivered experiment and holdout definitions.

Audit date: 2026-09-20. This is a research-design report, not a completed catalytic mechanism. Evidence labels identify the status of substantive statements. A literature fact means what a source reports, not independent experimental replication.

## 1. Charter and feasibility decision

**[WORKING HYPOTHESIS]** Retain the CuN3P1/CuN4 catalyst family provisionally. Restrict the chemistry to quinazolinone–aryl-nitrile reductive coupling/ring opening. The central question is whether potential-dependent competitive occupation is necessary to explain the advantage of CuN3P1 over CuN4, and under which conditions that advantage decreases or disappears. Study reversible occupation/protonation first; activate reconstruction/leaching branches only when supported.

**[LITERATURE FACT]** The publisher reports CuN3P1, gram-scale chemistry and modification of 11 natural products. [Primary article record](https://doi.org/10.1002/anie.202505085) An author-institution scheme identifies an aryl-nitrile partner and TBAPF6/TEOA/DMF at 10 mA with a carbon anode; natural alcohols are first esterified with 4-cyanobenzoic acid using DCC/DMAP. [Author-institution report](https://neunews.neu.edu.cn/info/1341/941741.htm)

**[MODEL INFERENCE]** Both quinazolinone Q and nitrile B must therefore enter the model. Q-only adsorption calculations cannot establish which partner first accepts an electron. The public scheme is not a complete procedure. Concentrations, electrode area, time, working potential, exact product atom mapping, electron stoichiometry and published catalyst coordinates remain gates until verified from sufficient original material.

**[WORKING HYPOTHESIS]** The application is small-batch diversification of natural-alcohol-derived fine chemicals for discovery samples. It is a platform opportunity, not an established pharmaceutical manufacturing process. The decision metric is acceptable isolated product per electricity/separation burden. Core catalyst controls are CuN4 and CuN3P1; NC, NP@NC and bare electrodes are necessary experimental controls. Dissolved copper and clusters are alternative explanations, not additional screening programmes.

**[WORKING HYPOTHESIS]** Success requires rejection or constraint of one mechanistic alternative, a robust prospective trend, or a justified negative result that changes a design decision. Finding local minima or increasing calculation count is insufficient.

## 2. Novelty boundary

**[LITERATURE FACT]** Original Cu electrocatalysis studies in a different reaction demonstrate potential-driven coordination changes and pulse-dependent populations of atoms, clusters and particles. [2023 study](https://doi.org/10.1038/s41467-023-40970-y); [2024 study](https://doi.org/10.1038/s41467-024-50379-w).

**[MODEL INFERENCE]** These motivate a working-state audit but do not establish reconstruction in DMF/TEOA organic synthesis. Their aqueous potentials, RHE values and CO chemistry cannot be transferred. Recalculating published adsorption is reproduction; applying established ensemble theory is application; excluding a fixed-site interpretation with independent observables may establish mechanism; predictions on frozen holdouts establish prospective performance only after experiment. No priority claim is justified. An unmentioned control in accessible material is not proof that the original SI lacks it.

## 3. Competing hypotheses

**[WORKING HYPOTHESIS]** H1–H3 may coexist. The study must identify necessary contributions rather than force one universal winner. No direct contradictory result for these hypotheses in the target system has yet been independently verified.

| Hypothesis | Existing support and gaps | Required assumptions | Discriminating calculation | Discriminating experiment | Substantial weakening observation |
|---|---|---|---|---|---|
| H1: static coordination | The article proposes coordination/adsorption effects; operational persistence is not established in the accessible material | Comparable site numbers and stable coordination; transport controlled | Identical balanced Q/B adsorption and confirmed key-step comparisons on both sites; model/method sensitivity | Matched-potential initial kinetics with working-state characterization | Original coordination is a minor inactive population; fixed-site model misses a reproducible ranking reversal beyond its uncertainty |
| H2: dynamic occupation/speciation | External Cu studies make it plausible; no target-system population measurement verified | States are kinetically accessible; fast equilibrium distinguished from slow reconstruction | Balanced Q/B/product/DMF/TEOA/H-state free energies, then only relevant interconversion barriers | Potential steps and recovery, time-resolved rate/selectivity, operando spectra if available | Small bounded population changes cannot explain the rate, and ensemble predictions fail while a fixed-site model succeeds |
| H3: interfacial transfer | Electrode mediation is the proposed framework; target KIE and sampled barriers not verified | Proton source/activity and charge/potential treatment known; TEOA role established | Compare Q and B reduction response; balanced HA/A− transfer and targeted solvent configurations; ET/PCET theory only after checking assumptions | Donor variation and H/D substitution with matched medium | No discriminating response after occupation/transport controls weakens a specified rate-controlling H-transfer model, not every possible H3 mechanism |
| A1: leached Cu/clusters | External precedent; no verified positive target-system evidence | Mobile species survive or redeposit; detection limits constrain catalytic contribution | Only evidence-triggered, balanced dissolution/recoordination and minimal cluster checks | ICP time series, electrode/solution cross-transfer, dissolved-Cu spikes, operando Cu–Cu response | Quantitative upper bound is insufficient to account for activity; a negative filtration result alone is inadequate |
| A2: support contribution | Requires blank controls; doping may co-change wetting/porosity | Comparable supports and films | Only the most relevant support-site comparison if required | NC, NP@NC and bare electrodes with matched preparation | Support controls cannot produce relevant product within detection limits and metal-state perturbations explain rates |
| A3: product inhibition | A testable consequence of adsorption/desorption models, not demonstrated here | Product is stable and quantifiable; spikes do not change medium | Competitive product/reactant adsorption and release | Product spikes plus independent Q/B orders at low conversion | No inhibition where the constrained model predicts strong blocking |
| A4: transport/electrode architecture | Porous electrodes and constant current make this necessary to exclude | Comparable hydrodynamics and effective films | Transport/depletion bounds, separately from DFT | Stirring/flow and thickness/loading series; product partial current | Chemistry remains insensitive over a meaningful transport perturbation; one stirring setting is insufficient |

## 4. Minimal calculation sequence and decisions

**[PLANNED CALCULATION]** Execution status comes from individual calculation manifests, not this table.

1. **Identity gate:** verify Q, B, product and co-reactant connectivity/atom mapping; establish TEOA function, both half-reactions and electron count z. Without this, do not calculate target FE or absolute process energy.
2. **Molecular pilot:** verified Q and a simple aryl nitrile, limited conformers and fixed-charge neutral/reduced states with spin checks. Low-cost optimization plus selected DFT points measures feasibility and whether electron-response ordering is method-sensitive; it does not provide an interfacial ET barrier.
3. **Published baseline:** reconstruct both actual catalyst models and the explicitly documented adsorption comparison. An invented finite fragment is not reproduction of an inaccessible periodic model. Freeze tolerance before reviewing results; a proposed 0.10 eV adsorption-difference tolerance is an audit choice, not a claim of DFT accuracy.
4. **Small ensemble:** two sites times bare/Q-bound/B-bound/one relevant DMF-or-TEOA-bound/one chemically justified H or ligand-protonated state. Verify formal balance and viable spins. Select the few states capable of changing the decision.
5. **Sensitivity:** change model boundary, an electronic-structure choice and a few solvent configurations only for close competing differences. If uncertainty exceeds the interpreted effect, stop ranking.
6. **Minimal mechanism:** add only a necessary chemical, release or interconversion step that separates the remaining alternatives. Validate transition states with vibrations and connectivity. Stop if accessible observations cannot discriminate.

**[PLANNED CALCULATION]** Comparing total energies of different N/P compositions cannot rank catalyst stability. Compare balanced adsorption reactions on each scaffold. For composition-changing states use ΔΩ=ΔG−Σνkμk with explicit reservoirs and a consistent electron reference. Use the actual nonaqueous HA/A− pair; aqueous pH is not a DMF acidity model. Fixed-charge calculations are not constant-potential calculations, and a static TS is not a sampled solution free-energy barrier.

## 5. Kinetics, uncertainty and prospective status

**[PLANNED CALCULATION]** Begin with reversible *+Q⇌Q*, *+B⇌B*, and evidence-required product, solvent and H states. Impose site balance. Only add a product-forming network after mapping and electron-path verification; retain both Q*+B and B*+Q/outer-sphere alternatives until excluded.

**[MODEL INFERENCE]** For demonstrated fast-equilibrating, single-site competitive binding, θ*=1/(1+KQ aQ+KB aB+KP aP+KL aL+…). This minimal model can fail for coadsorption, multiple sites, slow reconstruction or outer-sphere mediation. Detailed balance requires consistent activity standards and rate units. Barrier uncertainty must preserve forward/reverse thermodynamics; method errors shared by states should not be treated as independent draws.

**[PLANNED CALCULATION]** Predict relative rates, orders, inhibition and state fractions. Do not present reliable absolute current density without active-site density, electrode structure and transport. Calculate transition-state rate control by varying the common transition-state free energy while holding reaction free energy fixed; report intermediate-energy sensitivity separately.

**[PLANNED CALCULATION]** Register a minimum of three held-out condition combinations before fitting: the same Q/B pair at two verified stable potentials and two donor activities, plus a quantitatively defined product-spike condition. Numerical values must come from the actual procedure/window. Holdouts cannot select the functional, reaction network or fit parameters. Freeze predictions and intervals before experiment. If substrate predictions are later justified, split by whole scaffold rather than conformer or repeated job.

**[WORKING HYPOTHESIS]** Current priorities are candidate discriminating predictions, not validated performance rankings: (1) if H2 is necessary, potential/donor changes alter the P advantage alongside state populations; (2) if P stabilizes bound product, product spiking may suppress CuN3P1 more strongly and separate initial-rate from endpoint-yield rankings; (3) if H3 dominates, donor/isotope perturbations distinguish catalysts after medium/occupation controls. Directions and magnitudes are not yet robust. No ML is justified before reliable labels, meaningful candidates and a decision are established; later prefer small-data uncertainty models and simple baselines.

## 6. Five high-value experiments

**[REQUIRES EXPERIMENT]** Collaboration, instruments and measurements are unconfirmed.

| Priority | Experiment and purpose | Controls | Contrary outcome means |
|---|---|---|---|
| 1 | Matched-electrode low-conversion kinetics plus transport scan; establish whether P advantage is chemical | Loading, binder, thickness, area, reference, conductivity, flow/stirring; NC/NP@NC/bare; capacitance is not site count | Ranking disappears with flow/thickness: fix transport/electrode design before more DFT |
| 2 | Potential-step/recovery with product kinetics; optional operando Cu spectra | Capacitive transient and mass transport, dwell time, iR logging, with/without substrate, spectral detection limits | Current transients alone do not establish reconstruction; negative bulk spectra only constrain major populations |
| 3 | Donor and H/D variation; test specific H-transfer/occupation models | Confirm TEOA function; water content, effective acidity, conductivity, viscosity, H/D exchange | No effect weakens a specified rate-limiting H-transfer mechanism; does not exclude H participation |
| 4 | Product spike and separate Q/B reaction orders | Low conversion, product stability, fixed salt/volume/solubility | Missing predicted inhibition falsifies the blocking model rather than inviting unconstrained extra parameters |
| 5 | Solution/electrode cross-transfer, ICP and low-level Cu spikes | Time after potential removal, washed electrode, transfer blanks, detection limits, redeposition; Cu–Cu characterization if needed | Active filtrate may redeposit; inactive filtrate may deactivate; interpretation requires quantitative/time evidence |

## 7. Application and identifiable route comparison

**[LITERATURE FACT]** A more closely matched, identifiable scaffold route appears in original patent [US7307088B2, Scheme 1 / General Synthetic Procedures](https://patents.google.com/patent/US7307088B2/en): prepare a 2-aminobenzamide derivative, then use aldehyde/NaBH4 reductive amination to obtain substituted anthranilic amides. This citation establishes a synthetic route, not the bioactivity, manufacturing use or legal status of the Cu products.

**[WORKING HYPOTHESIS]** If a high-resolution original scheme/SI confirms the Cu product as the visually suggested 2-(arylmethylamino)-N-methylbenzamide, this is the preferred comparator: electrocoupling of Q/B versus the corresponding anthranilamide plus aryl aldehyde. Compare the same final product including all feedstock preparation, natural-alcohol esterification and isolation. The electrochemical route introduces electrode/electrolyte/full-cell-energy requirements; the comparator introduces chemical reductant and residue-removal requirements. No cost, waste or selectivity winner can be declared before matched measurements. Exact N-methyl/natural-product-ester scope is unverified.

**[LITERATURE FACT]** A verified non-electrochemical neighbouring route is hydrated ring opening/formylation of quinazolinones, producing N-arylformyl derivatives. [Original study, 10.1039/D2OB01234K](https://doi.org/10.1039/D2OB01234K).

**[MODEL INFERENCE]** This is an identifiable scaffold-editing comparator, not a verified route to the identical electrocoupling product. A rigorous same-product route ranking remains incomplete until exact product connectivity and an alternative preparation are verified. A route that merely makes quinazolinone starting material must not be misrepresented as a route to the ring-opened product.

**[PLANNED CALCULATION]** A testable same-product process alternative is to prepare the coupled carboxylic-acid intermediate first, then esterify the natural alcohol, instead of functionalizing the alcohol before electrolysis. This sequence is proposed, not literature-established; free-acid electrochemical incompatibility may rule it out. Both routes must include activation-reagent and purification burdens.

| Metric | Evidence class now | Required definition/data |
|---|---|---|
| Selectivity and isolated quality | Requires experiment; published scope is source-reported only | Calibrated Q/B/P and byproduct balance, analytical/isolation purity, repeat batches |
| FE | Unavailable; derived only after z is established | FE=zFnP/Q, distinct from chemical or isolated yield |
| Partial current density | Experimental/derived | jP=FE Itotal/Ageo; 10 mA is not mA cm−2 |
| Concentration | Not verified here | Q/B/salt/TEOA/water concentrations and solubility |
| Space–time yield | Derived under explicit assumptions | Acceptable mass/(working volume × electrolysis time); report external inventory/downtime |
| Lifetime | Requires experiment | Cumulative operation, throughput, drift, regeneration and replacement |
| Cu leaching/product contamination | Requires experiment | Solid/liquid/product metal balance and detection limits; no invented pharmaceutical specification |
| Solvent/electrolyte recovery | Requires experiment | DMF/salt/TEOA recovery and reuse, degradation and make-up streams |
| Separation | Requires experiment | Preactivation, extraction/chromatography/crystallization, residual reagents, water/solvent per mass |
| Continuous flow | Working hypothesis | Adhesion, pressure drop, wetting, residence time, fouling, transport and separation compatibility |
| Electricity/acceptable product | Derived scenario, actual value unavailable | Full-cell voltage, charge, z, FE, recovery and molecular mass; auxiliary/separation energy separate |
| Production cost | Unavailable | Not inferable from DFT without balances, lifetime and separation data |

## 8. Full-cell energy and transparent sensitivity

**[MODEL INFERENCE]** Use `E_cell[kWh/kg acceptable] = integral[V_cell(t) I(t)dt]/(3.6e6 m_acceptable_kg)`.

At approximately constant voltage, with product-generation FE and a separate acceptable-isolation recovery fraction f_rec, `E_cell = z F V_cell/(3600 FE f_rec M_g_per_mol)`, where F=96485.33212 C mol−1. If FE already uses acceptable isolated product, do not count the same recovery loss twice. The same electron current connects both half-reactions and is not doubled. Full-cell voltage cannot be replaced by cathodic potential, one overpotential or a DFT barrier. Add pumping, cooling and separation energy for total process energy.

No absolute kWh/kg is justified without verified z and process inputs. A dimensionless scenario is `E/E0=(V/V0)(FE0/FE)(f0/f_rec)`. Voltage at 0.8 of baseline gives 0.80 energy; FE at 0.8 gives 1.25; recovery at 0.8 gives 1.25; jointly 1.2 voltage, 0.8 FE and 0.8 recovery gives 1.875. These are algebraic sensitivity examples, not chemistry data or a cost forecast.

## 9. Limitations and stop conditions

**[EXECUTED RESULT]** This subtask audited accessible publisher and author-institution material and produced this design. The Cu full-text endpoint redirected to the abstract and SI returned 403; neighbouring-route SI access returned 502/504. No inaccessible structure, electron stoichiometry, barrier or measurement was fabricated.

**[PLANNED CALCULATION]** Stop production for unreproduced baseline, unresolved Q/B/P mapping or balance, method/model uncertainty exceeding the effect, conflict with experimental conditions, indistinguishable observable predictions, or inaccessible catalyst structures. A local minimum proves neither prevalence nor synthetic accessibility. A small fixed-charge molecular pilot is useful engineering/feasibility evidence, not reproduction of the catalyst interface. Replace provisional predictions only with manifest-linked actual evidence.
