# Potential-dependent occupation and selectivity in Cu–N–P organic electroreduction

Execution date: 20 September 2026 (Asia/Shanghai). Stage: evidence audit and local molecular preflight. **The catalytic mechanism and prospective performance have not been established.** Evidence labels distinguish source reports, hypotheses, calculations and experimental needs. This report follows the ten contractual deliverables in order.

## 1. Project charter

**WORKING HYPOTHESIS — Central question.** For the same verified quinazolinone–aryl-nitrile coupling in DMF, does CuN3P1 retain its advantage over CuN4 when potential-dependent occupation by both reactants, donor/solvent and product is included, or can competition erase/reverse the advantage? The principal uncertainty is the identity and population of the working occupied site. Full reconstruction and interfacial transfer are competing explanations to investigate only when needed to distinguish observations.

**LITERATURE FACT.** The source describes reductive *coupling* and ring opening. An author-institution reaction scheme identifies N3-methylquinazolin-4-one, aromatic nitrile, DMF/TEOA/TBAPF6, 10 mA, carbon anode and CuN3P1@NC cathode. These partial conditions do not define a working electrode potential, current density, reaction time, concentration or net electron balance. [Cu article](https://doi.org/10.1002/anie.202505085); [author-institution scheme](https://neunews.neu.edu.cn/info/1341/941741.htm).

**WORKING HYPOTHESIS — Scope.** One reaction family; CuN3P1/CuN4 catalyst pair; competitive occupation as one principal uncertainty; natural-alcohol-derived discovery/fine-chemical diversification as one platform application. Include NC, NP@NC and bare-electrode blanks experimentally. Treat dissolved Cu, clusters, product inhibition and transport as alternatives, not extra catalyst-screening programmes. The natural-alcohol application is reported at laboratory/gram scale; manufacturing use is unverified. No contact, collaboration, private data or instrument access has been established.

**MODEL INFERENCE — Feasibility decision.** Retain this system provisionally. Molecular connectivity permits local preflight, but inaccessible original Cu full text/SI, missing catalyst coordinates/protocol, unresolved product mapping and unverified server resources block published-model reproduction. The narrower question appears testable; whether it remains unresolved after full SI inspection is itself a gate. Do not substitute an invented graphitic fragment for the published model.

**WORKING HYPOTHESIS — Competing explanations.**

| Hypothesis | Existing support / contradictory evidence | Assumptions | Discriminator and falsifier |
|---|---|---|---|
| H1 static coordination | Authors attribute improvement to substrate-N binding; no independently verified target contradiction | Comparable accessible site counts, persistent coordination, kinetic control | Reproduce balanced adsorption and a decision-relevant step; matched-potential kinetics. A robust reversal unexplainable by fixed-site kinetics weakens H1 as a sufficient explanation. |
| H2 working-state occupation | External Cu studies demonstrate possible coordination changes; target populations unmeasured | States are chemically reachable; equilibrated populations distinguished from slow restructuring | Compare Q/B/product/donor states with reservoirs and potential; potential-step/recovery plus population-sensitive measurements. Negligible allowed population change and failed ensemble predictions weaken the specified H2 model. |
| H3 transfer/interface | Redox mediation motivates ET/PCET alternatives; target sampled ET barriers/KIE not verified | Proton source, activity and TEOA function established; appropriate ET regime | Q versus B electron acceptance; donor/isotope perturbations with medium controls. A null response rejects a specified rate-controlling H-transfer model only when the experiment has discriminatory power. |
| A: mobile Cu, support or transport | Plausible alternatives, no verified positive target evidence here | Detection limits and electrode preparation constrain contributions | Quantitative Cu balance/cross-transfer; support controls; flow and film-thickness tests. Activity that follows transferred solution or collapses under transport matching requires redesign. |

Each row is a **WORKING HYPOTHESIS** informed by the cited literature, not a discovered mechanism. The expanded six-field hypothesis matrix is in [the design report](DESIGN_PROCESS_EN.md).

**PLANNED CALCULATION — Success and stop criteria.** Success requires one distinguishable mechanism, a method-robust held-out prediction, or a defensible negative result that changes design. Stop production for inaccessible/unreproducible baseline, chemically implausible models, method spread larger than the interpreted effect, non-identifiable mechanisms, or negligible decision value. The present preflight satisfies none of the full catalytic success claims.

## 2. Evidence matrix and novelty audit

**EXECUTED RESULT.** Five mandatory DOI records and later work were investigated. The combined [evidence matrix](../literature/evidence_matrix.md), [machine-readable matrix](../literature/evidence_matrix.json), [Cu audit](../literature/cu_evidence.json) and search/access records retain source type and gaps. The main Cu bibliographic record is *Phosphorus-Doped Single Atom Copper Catalyst as a Redox Mediator in the Cathodic Reduction of Quinazolinones*, Angew. Chem. Int. Ed. 64, e202505085 (2025), version of record 27 March. Corresponding authors are Wen-Hao Li, Hai-Tao Tang and Dingsheng Wang; Ying-Ming Pan is a coauthor, not a corresponding author in this record.

**LITERATURE FACT.** The other mandatory anchors are CO2-mediated chlorine evolution (Nature), Mn-mediated silane oxidation, Fe-mediated organic electrooxidation, and Mn–Rh relay catalysis. The last concerns **thermal cumene oxidation**, not electroreduction. [Nature](https://doi.org/10.1038/s41586-023-05886-z), [Mn](https://doi.org/10.1002/anie.202315032), [Fe](https://doi.org/10.1002/anie.202404295), [Mn–Rh](https://doi.org/10.1039/D4SC08658A).

**EXECUTED RESULT.** The matrix contains 12 original studies including the core paper: four mandatory conceptual anchors and seven additional catalyst/interface/method neighbors. Only one verified study is the same quinazolinone system; no eight-study same-system body was found. This is a search/access shortfall, not proof of absence. Do not count aqueous CO2 or thermal oxidation as identical-reaction comparators. Under the stricter organic-electroreduction category, four studies are verified (four short of eight); eight are available only with mechanism/method neighbors. Reviews are recorded separately and excluded from the original-study count.

**LITERATURE FACT.** P-doped Cu single atoms were already used for organic furfural electroreduction in 2021. Thus P addition, single-atom Cu, and static adsorption improvement alone are not novel project claims. [Closest earlier Cu/P study](https://doi.org/10.1039/D0GC03999C). Later 2026 Cu-cluster/furfural and Zn/iodide-relay organic work are included as distinct-system developments; bounded searches did not verify a later same-Cu/quinazolinone study.

**MODEL INFERENCE — Novelty boundary.** Rebuilding the published adsorption model is reproduction. Applying known occupation thermodynamics to it is system application. Identifying and experimentally excluding a competing working-state explanation could be mechanistic discovery. An uncertainty-aware prediction protocol is methodology only if it improves decisions against simple baselines. Independent prospective experiments would be new validation. None is claimed as completed or “first”.

**EXECUTED RESULT — Access.** Public Nature, Mn–Rh, Cu-dynamic-site and Pd-interface material/SI were obtained and relevant portions inspected. Wiley Cu/Mn/Fe SI requests failed; accessible abstracts cannot supply their missing protocols. Third-party articles, figures and extracted full texts remain in local audit storage; public deliverables include our summaries, source URLs, access outcomes and hashes, not unlicensed copies. The literature index states the inspection scope rather than implying all-page scientific review.

## 3. Model specification

**EXECUTED RESULT — Molecular preflight model.** The 20-atom molecule C9H8N2O was independently generated from `CN1C=Nc2ccccc2C1=O` with RDKit ETKDGv3, seed 20260920, then MMFF94. Its 2D identity is supported by the author scheme and PubChem CID 17088. Coordinates are computationally generated, not published DFT or crystallographic coordinates. Neutral: charge 0, 84 electrons, zero unpaired electrons. Reduced: charge −1, 85 electrons, one unpaired electron. xTB's UHF input is number of unpaired electrons; it is not an independent validation of spin purity or a Cu oxidation state.

**EXECUTED RESULT — Protocol.** Nonperiodic fixed-charge GFN1/GFN2-xTB 6.7.1; native method parametrization/dispersion, no added D3/D4; no plane-wave cutoff or k mesh. One common optimized neutral gas geometry supplies all eight charge/method/medium single points. Media: gas and equilibrium ALPB(DMF); electronic SCC accuracy 0.1, at most 250 iterations, two threads. A gas GFN2 Hessian checks the neutral reference. No explicit ion, electrode, electrical field, electron reservoir, electrochemical calibration, ZPE/thermal correction or standard-state reaction free energy is applied. ALPB total energies include model solvation terms; they are not pure gas electronic energies.

**MODEL INFERENCE.** `E(q=-1)-E(q=0)` at the same nuclei is a charging diagnostic. Because solvent polarization is equilibrium ALPB, it is not a nonequilibrium vertical electron affinity. It is not a reduction potential, activation barrier or free energy of the actual coupling. The pilot contains no copper and cannot discriminate H1–H3 by itself.

**PLANNED CALCULATION — Published baseline.** First obtain original coordinates/cell, composition, spin/charge, functional/dispersion, basis/pseudopotentials, cutoff/k sampling, solvation, reference potential and thermal conventions. Freeze the same-method targets before execution: exact composition/connectivity, reproduced bond distances within 0.05 Å and relative adsorption values within 0.10 eV where the source precision supports this; convergence changes below 0.03 eV and force target ≤0.03 eV/Å or the paper's stricter criterion. These are reproduction tolerances, not physical error bars. Diagnose failure before new candidates.

**PLANNED CALCULATION — First Cu batch.** B01/B02 are four optimizations: each verified bare CuN4/CuN3P1 model and its published adsorbed reference, followed by source-matched single points. The decision is whether the published sign/scale reproduces. B03 adds Q, aromatic nitrile B, product, one relevant solvent/donor and justified protonated state only after B01/B02 pass; distinguish single occupation from necessary coadsorption. B04 examines two measured potentials bracketing onset, after reference calibration and server checks. Spins are enumerated by electron parity and electronic solutions, not assigned from a site label. No Cu production inputs have been fabricated.

**PLANNED CALCULATION — Ensemble and interface.** Compare adsorption through balanced reactions on each scaffold, not raw energies of different N/P compositions. Use ΔΩ=ΔG−ΣΔniμi with explicit Q/B/product/DMF/TEOA or HA/A− and electron reservoirs. Establish activity conventions and measured nonaqueous reference calibration. If using constant-charge periodic calculations, label them so and apply controlled potential mapping/sensitivity; true constant-potential comparisons require a verified implementation and charge response. CHE is not an organic ET activation model. Marcus theory requires examination of weak electronic coupling, distinguishable diabatic states, solvent response, reorganization and relevant timescales. TS claims require the relevant imaginary mode and downhill connectivity/IRC; static TSs are not sampled solution free energies. Explicit solvent or sampling is triggered only by a concrete uncertainty. Test model boundary, electronic method and selected solvent configurations before interpreting a small P effect.

## 4. Pilot execution package and resource gates

**EXECUTED RESULT — WP0.** The inspected host is Windows 11, **i7-13700H**, 14 cores/20 logical processors, approximately 16 GiB installed RAM; this differs from the supplied i7-13900H description. About 7 GiB memory and 90 GB C-drive storage were free at the audit snapshot. CUDA was unavailable in the imported PyTorch. Existing RDKit/xTB/Psi4 were reused; no package installation or license purchase occurred. Psi4 1.11 imported but no Psi4 energy job was executed. Local scheduler commands were absent, no usable WSL distribution was verified, and no institution SSH host/allocation or license was provided. No institutional login or job submission occurred. [Hardware audit](../environment/hardware_audit.json), [versions](../environment/verified_versions.json).

**EXECUTED RESULT — First batch.** Ten native xTB calls completed: one neutral geometry optimization, eight single points and one Hessian. The optimized geometry had a −50.49 cm−1 mode and was **rejected as a minimum**. Inspection localized the mode mainly to methyl rotation. A documented +60° methyl rotation followed by verytight optimization, Hessian and eight repeated single points produced a second ten-call batch. All raw outputs, failed scientific acceptance and corrected parsing remain available.

**EXECUTED RESULT — Accepted preflight reference.** The repaired structure preserves the expected atom-distance connectivity; its 54 vibrational modes have a lowest frequency of 78.64 cm−1, in addition to six rigid-body zeros. This is only a minimum on the tested gas GFN2 surface. Twenty native calls took about **2.87 s summed subprocess wall time** with two requested threads; launch/analysis time is additional. Peak memory was not measured. This tiny-molecule timing cannot estimate a periodic Cu or interface job.

| Fixed-nuclei diagnostic / eV | GFN1 | GFN2 | Method spread |
|---|---:|---:|---:|
| Gas: E(anion) − E(neutral) | −5.43 | −4.65 | 0.78 |
| ALPB(DMF): same difference | −7.12 | −6.28 | 0.83 |

All entries above are **EXECUTED RESULT** at the stated model, with source IDs in [pilot_summary.json](../results/pilot_summary.json). The output precision is not a physical accuracy statement.

**MODEL INFERENCE — Updated decision.** Both method spreads exceed the preregistered 0.20 eV preflight flag. Retain xTB for bounded geometry preparation; withhold electrochemical/catalyst ranking from these charged-state values. This is a local method-sensitivity warning, not a proof that the paper's mechanism fails. Higher-level work must first pass the source-model gate.

**PLANNED CALCULATION — Three levels.** The following are conditional scheduling scenarios, **not measured HPC forecasts**. Baseline source access and four native baseline timing jobs must precede allocation. GPU need is zero for the planned CPU DFT route; actual software choice depends on entitlement and validated physics.

| Level | Required result and approximate work | Conditional resources and elapsed scenario | Main uncertainty / stop |
|---|---|---|---|
| Minimum viable | Reproduce baseline; constrain fixed-site versus occupation interpretation. Approximately 24–48 optimizations, 12–24 single points, 2–4 selected interconversion/step attempts and 6–12 vibrational/connectivity checks; prune before all are run | Illustrative 44–88 electronic-structure jobs, 32 CPUs, 64–128 GB/job, 50–300 GB retained storage. At an assumed 4–24 h/job: 5,632–67,584 allocated core-hours; two concurrent jobs about 4–44 compute days, plus queue/failures | Unknown actual cell, charge convergence, engine and license; stop if baseline or method gate fails |
| Standard | Minimum level plus 3 frozen condition predictions and required side-reaction/release checks, approximately 120–240 jobs total | Illustrative 32–64 CPUs/job, 64–256 GB, 0.2–1 TB; assumed 4–48 h/job gives 15,360–737,280 core-hours. Four concurrent jobs about 5–120 compute days, plus queue | Network identifiability and sensitivity; no unconstrained fits, no prediction if uncertainty dominates |
| Extended | Only if earlier decisions need solvent sampling or active learning; 4–8 short interface trajectories (10–20 ps each) or 20–50 selected labels, not both by default | At 0.5 fs, 80,000–320,000 total MD steps. Benchmark seconds/step first: with 64 CPUs and illustrative 1–30 s/step, roughly 1,400–171,000 core-hours; 128–512 GB and 0.5–3 TB conditional | Sampling convergence, appropriate potentials and resources. Stop if sampling has little decision value; no ML potential without domain validation |

The arithmetic describes reserved CPU-hours, not useful CPU efficiency. Timing/memory/storage will be replaced by measured first-job medians, tail factors and output sizes. Do not run these campaigns on shared login nodes or submit them without applicable allocation/authorization. The provided Slurm file is a small molecular pilot template only.

## 5. Mechanism comparison

**WORKING HYPOTHESIS.** Retain Q-first electron acceptance, B-first electron acceptance, adsorbed inner-sphere coupling, outer-sphere mediation, stepwise ET/PT, concerted PCET, subsequent cleavage and product release until evidence excludes them. HER, overreduction and partner consumption are selectivity controls. TEOA's donor/base/anodic role must be established rather than assumed from its presence. Product identity and both half-reactions must be balanced before energetics or FE.

**MODEL INFERENCE.** The author adsorption account supports testing H1 but does not by itself establish a kinetic step, working population or whether Q or B is first reduced. External Cu restructuring evidence motivates H2 but does not prove it in DMF. A potential-dependent yield or current alone cannot distinguish coverage, transfer and transport. Orthogonal rate-order, inhibition, state and mass-balance observations are needed. No alternative is presently eliminated by the molecular pilot.

## 6. Predictive results and uncertainty

**EXECUTED RESULT.** There is **no validated catalyst ranking, quantitative rate/selectivity prediction or absolute current density**. [prospective_holdout.json](../config/prospective_holdout.json) reserves three condition groups before model development: doubled B activity, 0.2-equivalent product preaddition and halved calibrated donor activity. Baseline quantities/reference and chemical compatibility are unresolved. All related repeats/calculations remain in the same held-out group. The config bytes are hashed in the delivery manifest.

**WORKING HYPOTHESIS — Ranked by information value, not performance.** (1) If competitive B occupation matters, changing B activity should alter Q order and the P/non-P rate ratio; sign is not yet determined. (2) If P stabilizes product disproportionately, product spikes should inhibit CuN3P1 more strongly; failure weakens that particular blocking model. (3) If proton delivery is controlling after coverage corrections, donor/isotope perturbations should produce distinguishable responses; other medium changes remain alternatives. Uncertainty is currently structural and unquantified; numerical confidence intervals would be fabricated. Freeze numerical predictions before holdout measurement only after sensitivity passes.

**PLANNED CALCULATION.** For a justified physical kinetic network, impose atom/electron and site balances, reversible adsorption/release, explicit activities and justified potential dependencies. Shared transition-state/state energies enforce detailed balance; perturb a shared TS to evaluate rate control and propagate correlated energy uncertainties. The delivered three-state solver and synthetic diagnostics test these principles only; they omit the actual two-reactant mechanism and cannot generate target kinetics. WP4 has not run on physical labels. WP5 remains deferred: no trustworthy candidate dataset or learning decision exists. DFT labels would not constitute experimental validation.

## 7. Experimental discrimination plan

All five entries are **REQUIRES EXPERIMENT**; availability is unconfirmed. [Detailed controls and contrary-outcome interpretation](DESIGN_PROCESS_EN.md).

1. **Initial kinetics with transport/electrode controls:** match catalyst loading, film, electrolyte and measured potential; vary flow/stirring and thickness, include NC/NP@NC/bare electrodes. A disappearing P advantage redirects effort to electrode architecture.
2. **Potential step and recovery:** couple calibrated Q/B/product analysis to time-dependent response; add operando Cu characterization only if available. Hysteresis plus state changes can support a dynamic model; current transients alone cannot.
3. **Donor/H–D response:** establish TEOA function, control water/activity, viscosity, conductivity and isotope exchange. A resolved null effect weakens a specified H-transfer-limited model, not all proton involvement.
4. **Product spike and separate Q/B orders:** use low conversion, verified product and constant medium. Failure of predicted inhibition falsifies the blocking hypothesis instead of licensing extra free parameters.
5. **Cu mass balance and cross-transfer:** timed solution/electrode swaps, washed-film controls, ICP detection limits and dissolved-Cu spikes. Filtrate activity may involve redeposition and a negative result may reflect deactivation; quantify these alternatives.

## 8. Industrial relevance assessment

**LITERATURE FACT.** Natural alcohols are derivatized through 4-cyanobenzoate precursors in the author scheme; gram scale and recycling are reported. These do not provide lifetime, product metal content, productivity or manufacturing economics. An identifiable non-electrochemical reference is the anthranilamide/aldehyde reductive-amination route in [US7307088B2](https://patents.google.com/patent/US7307088B2/en). It is a conditional near-product scaffold comparator until a high-resolution product map establishes identity; no same-product process superiority is claimed. A separate published quinazolinone ring-opening/formylation route gives different products. [Original study](https://doi.org/10.1039/D2OB01234K).

**MODEL INFERENCE.** The application is small-batch fine-chemical diversification, not an established drug-manufacturing route. Compare total acceptable isolated product and separation burden, including DCC/DMAP prefunctionalization, DMF, salt/TEOA recovery and work-up. A reversed order of esterification/electrocoupling is only a proposed alternative and may be incompatible with free acid.

| Metric | Current evidence classification |
|---|---|
| Molecular charging diagnostic | Directly computed; no process performance meaning |
| Reported gram scale, natural-product scope, five recycling runs | Literature fact, source-reported; underlying Cu SI not independently inspected |
| Selectivity, assay/isolated purity, substrate concentration, FE, lifetime, leaching/product metal | Requires experiment or currently unavailable in verified target procedure |
| Partial current density | Derived only from measured FE × current / specified area; 10 mA is not current density |
| Space–time yield | Derived from acceptable isolated mass/(working volume × time); inputs unavailable |
| Solvent/salt recovery, separation burden and flow compatibility | Requires experimental mass balances, electrode adhesion/fouling and recovery data |
| Electricity per acceptable product | Scenario-derived only; full-cell voltage, FE, z, molecular mass and recovery unknown |
| Production cost | Currently unavailable; cannot be inferred from DFT |

**MODEL INFERENCE — Energy.** For approximately constant full-cell voltage,

`E[kWh/kg acceptable] = z F Vcell / (3600 M[g/mol] FE f_recovery)`.

Use the actual integral of Vcell(t)I(t) when voltage varies. FE uses product formed; recovery accounts for acceptable isolation, and must not be applied again if already included in FE. Electron current is not doubled for two half-cells. [18 scenario combinations](../results/process_scenarios.csv) vary Vcell=2/4/6 V, FE=0.3/0.6/0.9, recovery=0.7/0.9 as **assumed scenarios**, reporting only the coefficient multiplying unknown z/M. They are not substitute experimental data. Pumping, cooling, catalyst manufacture and separation are excluded and must be added for process energy.

## 9. Reproducibility package

**EXECUTED RESULT.** The additive project directory preserves the existing repository's pincer work. See [README](../README.md) for commands, [raw-output index](../data/raw/index.json), [configuration](../config/pilot.json), [environment](../environment/verified_versions.json), per-calculation `record.json` and immutable raw inputs/stdout/stderr/native files. Each calculation records identifier, structure provenance/hash, settings, native version, requested resources, status, convergence, parsing and downstream links. Failed scientific acceptance is preserved even when native execution succeeds.

**EXECUTED RESULT.** Scientifically relevant checks cover units, explicit reference conversion, electron/atom balance, complete-output parsing, detailed balance, site conservation, FE/recovery arithmetic, grouped holdout separation and native evidence integrity. The tests use explicitly synthetic numerical fixtures where appropriate; no fixtures are passed off as chemical observations. Independent code review is included. Regenerated summaries are derived from the raw calculations; reruns use fresh immutable output folders.

**PLANNED CALCULATION.** Server probing is read-only and awaits legitimate access. The package provides a scheduler template but makes no remote-access claim. Missing published periodic coordinates prevent a meaningful ready-to-submit Cu deck; the production-gate manifest records exactly which jobs remain unexecuted and why.

## 10. Limitations, negative results and next decision

**EXECUTED RESULT.** Original model reproduction: unexecuted. Cu ensemble/interface/TSs: unexecuted. Physical kinetics: unexecuted. Quantitative prospective validation/ML/industrial trial: unexecuted. Local molecule preflight: 20 completed native calls, one rejected initial minimum, one repaired minimum, unacceptable cross-Hamiltonian charging spread. A preparation-path bug and an unexplained numerical-library test exit are recorded in [negative_results.md](../logs/negative_results.md); neither is hidden as a successful calculation.

**MODEL INFERENCE.** Three major conclusions have distinct evidence and alternatives: (i) the reaction must include aromatic nitrile because the author scheme shows it; exact atom mapping still requires SI; (ii) the fixed-charge preflight is method-sensitive because independently parsed GFN1/GFN2 values disagree beyond the declared gate; calibration or higher-level theory could change this diagnostic; (iii) working-state occupation is worth testing, but static coordination and transfer/transport remain viable because no target state/rate evidence yet separates them. Source-model recovery plus the first matched-potential/occupation experiments would distinguish these alternatives.

**PLANNED CALCULATION — Next minimal action.** Obtain lawful Cu original/SI and model coordinates; verify Q/B/product mapping and reference procedure; audit actual server allocation; then run B01/B02 and measure cost. Only after reproduction and uncertainty gates pass should B03/B04 or kinetic prediction begin. No publication, mechanism confirmation or industrial deployment is promised.


