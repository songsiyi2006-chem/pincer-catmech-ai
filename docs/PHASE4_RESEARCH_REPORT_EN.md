# From a mechanism prototype to falsifiable research: Pincer-CatMech-AI Phase 4 technical report

16 September 2026. This research programme is organized around the originality, importance and evidential expectations of Nature/Science/Cell (CNS), with Nature/Science providing the relevant disciplinary scope. The report does not claim that the project is ready for submission.

## 1. Research assessment and a narrower question

The project currently provides auditable software, several genuine low-cost quantum-chemical calculations and explicit failure records. It does not yet establish an independently validated catalytic mechanism, a complete free-energy network for physical rate predictions, or an AI model that outperforms a simple baseline. The user has no supercomputing access and has supplied no private experimental raw data. This phase therefore delivers a revised research design, local diagnostics, an entry point based on a public experimental structure, and offline HPC bridging tools. It must not be described as a completed experiment-theory validation cycle.

Nature emphasizes originality, importance and relevance beyond a narrow specialist audience. Rigour is necessary but insufficient: more calculations, modules or pages do not establish these contributions. The journal does not publish a checklist that guarantees acceptance. These expectations guide research decisions; the journal name cannot substitute for a scientific conclusion. [Nature journal information](https://www.nature.com/nature/journal-information); [editorial criteria and processes](https://www.nature.com/nature/for-authors/editorial-criteria-and-processes).

Broad novelty claims in the original plan must be withdrawn. K/Na- and base-dependent amine/imine selectivity in Mn borrowing-hydrogen chemistry has direct precedent, as does high-throughput xTB/DFT screening of Mn pincer complexes. Ligand noninnocence and alcohol-assisted proton transfer do not become new discoveries when combined as software modules. Different optimized literature conditions are not single-variable causal experiments; condition tables and dedicated kinetic controls must be distinguished. [ACS Catalysis 2018](https://doi.org/10.1021/acscatal.8b02530); [ZAAC 2021](https://doi.org/10.1002/zaac.202100078).

The proposed question is narrower: in a structurally traceable Mn-PNP borrowing-hydrogen system, can explicit ion pairing, proton transfer and catalyst deactivation jointly explain condition- and time-dependent selectivity, and predict held-out conditions better than simpler models? This is a hypothesis. Neither its correctness nor its broader novelty has been established. The literature audit and falsifiable design are described in NATURE_RESEARCH_CASE.md.

Ligand proton transfer or metal-ligand cooperation does not automatically establish redox noninnocence; the latter requires separate electronic-structure and experimental evidence.

### Concrete requirements for the CNS innovation ambition

The user's further requirement makes flagship-level originality and broader significance project gates. Nature/Science provide the relevant disciplinary scope; Cell primarily addresses important biological questions and cannot be made appropriate by adding a biological label. There is no shared checklist that guarantees CNS acceptance. The current Science author page was inaccessible in this check; the charter explicitly identifies its AAAS guide as the official 2019 edition, without claiming to verify all 2026 policies.

**The single primary innovation candidate is a shared, measurable solution-species mechanism connecting ionic/exchangeable-proton distribution, amine/imine branching and reversible dormancy within one borrowing-hydrogen system, with prospectively frozen intervention and second-scaffold predictions.** K/Na selectivity, base-mediated recovery, hidden states and reaction-progress fitting already have precedents. [Dormancy/base-recovery precedent for the same type of Mn precursor](https://doi.org/10.1039/D5GC05072C); [competing kinetic models and identifiability precedent](https://doi.org/10.1021/jacs.4c10488).

The required evidence chain identifies shared species in the same system, uses quantitative inventories/timescales and independent interventions to reject separate explanations, freezes network rules and structure-to-parameter mappings, and then tests complete trajectories and predictive intervals for unread conditions and a second ligand scaffold. A second crystal molecule of the same precursor is not scaffold transfer. If each new condition requires arbitrary refitting, or separate simpler models predict all results equally well, the coupling or transfer claim is withdrawn.

The claim/prior-art/new-finding/falsification/baseline/stop-rule matrix is in CNS_INNOVATION_CHARTER.md. Targeted review of 12 original studies is not exhaustive novelty clearance. Even perfect model checks and experimental execution do not replace a discovery that changes understanding; that central discovery remains absent.

## 2. Evidence levels and existing results

| Component | Existing evidence | Unsupported conclusion |
| --- | --- | --- |
| Software verification | 556 passing tests in the frozen Phase 3 release | Chemical correctness or publication readiness |
| Spin diagnostics | 48 actual attempts, 8 SCF convergences, 6 states passing the original spin and identity gates | Validated physical spin ordering or spin crossing |
| MECP | 0 eligible original state pairs; 0 accepted MECPs | A crossing channel controls the rate |
| Microsolvation | 9 xTB/ALPB structures; 3 selected clusters passed full-Hessian minimum checks | Dehydration TSs, solution populations or a general proton wire |
| Neutral dehydration pathway | 400 NEB steps across two seeds and one continuation; unchanged force threshold not reached | Accepted activation free energies |
| Kinetics | Numerical verification of 18 species and 18 channels; 40 fixture grids | Measured TOF, lifetime or catalyst rankings |
| EGNN | Test MAE 0.19520 eV; zero-difference baseline 0.19400 eV | Improvement over the baseline or validated barrier prediction |

These are frozen Phase 3 results. New diagnostics are not inserted into the old dataset to alter its historical success rate. Rate parameters used in software fixtures are not fitted literature measurements. Resource-limited termination and scientific acceptance failures remain failures.

Chemical identity is a separate gate. The Fe_bipyridine_pnnoh_iPr prototype has a hydroxypyridine O-H/O− response site; its spin calculations cannot be presented as validation of the user's proposed NH activation. Each scaffold requires an audit of its actual deprotonation site, proton count, total charge, counterions and reaction atom mapping. The amount of tBuOK added is not the activity of free tBuO−.

## 3. Phase 4 calculations and their limits

Values in this section are generated from archived results and checked before delivery. The quantum pilot uses an existing geometry, the same charge and unchanged SCF thresholds. A STO-3G quintet receives a longer time budget; two PBE/def2-SVP jobs diagnose singlet and quintet behaviour. The configuration is 2 threads, 500 MiB of program memory, a 1536 MiB process RSS limit and a total wall budget of 1200 s. These are gas-phase single-point energy and analytic-gradient diagnostics, not geometry optimizations, frequencies or high-accuracy spin gaps in toluene.

| Gas-phase single-point diagnostic | SCF status | S² | Spin-quality gate |
| --- | --- | --- | --- |
| sto3g quintet recovery | converged | 6.193760 | reject |
| def2svp singlet | converged | -0.000000 | pass (diagnostic only) |
| def2svp quintet | timeout | unknown | unknown |

The run attempted 3 jobs: 2 SCF convergences and 1 numerical spin-quality pass. Every record retains accepted_chemical_label=false and production_spin_gap_validated=false. Newly accepted MECPs, TSs and physical rates are all 0. Total energies from different bases must not be mixed, and an isolated passing state does not establish solution spin ordering. Native inputs, outputs, analytic gradients when produced, engine fingerprints and source snapshots are archived.

Thermochemical sensitivity analysis produced 405 rows from 3 previously accepted clusters, 9 temperatures, 3 qRRHO crossover wavenumbers and 5 neutral tBuOH activities. The source results, accepted geometries and spectra were checked by SHA256. All 27 association free energies at the original 100 cm⁻¹ crossover were reconstructed with a maximum difference of 0. This is reanalysis of existing quantum results; the number of new quantum jobs for this table is 0.

| Explicit tBuOH molecules | Standard association free-energy range at 383.15 K for crossover values of 50-150 cm⁻¹ | Range width |
| --- | --- | --- |
| 1 | +3.176 to +3.969 kcal/mol | 0.793 kcal/mol |
| 2 | +8.404 to +10.703 kcal/mol | 2.299 kcal/mol |
| 3 | +10.958 to +14.424 kcal/mol | 3.465 kcal/mol |

The selected clusters retain positive standard association free energies within the scanned conditions. This does not exclude other conformers, ion-assisted pathways or different standard-state behaviour. One selected minimum cannot represent a conformational ensemble. Reusing an ALPB electronic solvation contribution across temperatures is also not a temperature-dependent solvent model.

For a single association equilibrium $C + n A \rightleftharpoons CA_n$, with dimensionless activities referenced to 1 M:

$$\Delta G_{\rm eff}=\Delta G^\circ-nRT\ln a_A.$$

$$\ln(a_{CA_n}/a_C)=-\Delta G^\circ/(RT)+n\ln a_A.$$

A denotes neutral tBuOH, not total tBuOK equivalents. Lowering its activity raises the conditional association free energy shown here. The output is a conditional activity ratio for this equilibrium only. It does not solve bulk mass balances and is not a yield, catalytic rate, transition-state free energy or solution species fraction.

![Figure 1. qRRHO and neutral-alcohol activity sensitivity of accepted xTB/ALPB minima. Lines join deterministic post-processing results, not statistical confidence intervals; n is the explicit tBuOH count.](../examples/plots/phase4_solvation_sensitivity.png)

## 4. Public experimental structure and benchmark entry point

The 2016 Nature Communications Mn-PNP borrowing-hydrogen paper supplies a public experimental reference. Its SI points to the precursor study JACS DOI 10.1021/jacs.6b03709. The original complex 1 CIF hosted by ACS Figshare has been obtained; its molecular formula is C18H37BrMnNO2P2. These crystallographic results belong to the original researchers, not to this project. The original file, download link, hash, parsing rules and checks of the candidate structure are recorded separately. [Borrowing-hydrogen paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5059641/); [JACS precursor paper](https://doi.org/10.1021/jacs.6b03709); [original CIF](https://ndownloader.figshare.com/files/5467289).

Extraction checks expand the 124-atom asymmetric unit through four symmetry operations into a 496-atom unit cell. Eight finite, complete 62-atom molecules match the reported composition and Z=8. No partial occupancy or disorder was found. Periodic molecules were unwrapped, and the distance-derived connectivity matches the published bond table. Two crystallographically independent molecules, Mn1/Mn2, are retained. They represent the same precursor, not two catalytic systems or independent experimental replicates. Each Mn has a P2NC2Br environment, two CO ligands and one NH site. Original NH distances of approximately 0.836/0.846 Angstrom were retained without replacing them by optimized bond lengths.

The source CIF and derived coordinates retain the source's CC BY-NC 4.0 license; the software repository's MIT license does not relicense them. Charge 0 and multiplicities 1/3/5 are computational candidates only. The structure record still requires chemistry review and does not certify a solution-active species or true spin state. Six ORCA single-point candidate inputs were prepared and verified offline. All remain prepared_not_submitted; no ORCA result exists.

Crystal hydrogen positions may be constrained during refinement. A crystal geometry is not necessarily the active species in a hot solution. Even after extraction of a complete molecule, a chemistry researcher must check the NH site, Br coordination, CO count, protonation, charge and candidate spin. Starting from an experimentally corresponding geometry begins benchmark validation; it does not complete mechanism validation.

Reproduce a set of actual reference conditions and species before expanding to designed structures. The user-specified toluene, 383.15 K, baseline 0.05 equiv tBuOK and 0.01-0.20 equiv scan remain a proposed research window. Literature conditions are stored separately and must not be altered to match that proposal. Equivalents can be converted into initial concentrations only when concentration and reference-substrate definitions are available.

## 5. Evidence that distinguishes mechanisms

| Competing explanation | Discriminating observation | Confounding factors to address |
| --- | --- | --- |
| Total base alone controls behaviour | Initial rates, reaction orders and dependence on base amount | Activity, solubility, ion association and induction periods |
| Explicit cation association controls behaviour | K/Na and ion-association perturbations under matched conditions | Water or solvent changes accompanying cation replacement |
| Alcohol/water assists proton transfer | Additive effects at controlled activity; isotope and crossover responses | Equilibrium shifts, exchange and catalyst-solubility changes |
| Deactivation determines endpoint selectivity | Full concentration trajectories, restoration responses and species spectra | Substrate depletion, product inhibition and sampling losses |
| A spin-changing channel is necessary | Validated surfaces, crossing region, SOC and competing kinetics | Numerical contamination, different chemical species and multireference error |

This table is an experimental-design discussion for a research supervisor, not an unsupervised synthesis protocol. Experimental work requires an established laboratory, training and approved operating and waste procedures. Experienced personnel must select equipment, scale and analytical methods.

Shared measurements should include calibrated analyses, complete time series and material balances. Report internal standards and actual calibration curves, raw chromatograms or spectra, integration rules, detection limits, independent reaction repeats, errors, omitted data and failed conditions. Technical replicates do not replace independent experiments. Scale, replicate counts and power analysis should be determined from pilot noise measurements; no precision or sample size is fabricated here.

Register prediction conditions, endpoints, tolerated errors, baselines and stopping rules before data collection. Repeatedly changing parameters after seeing a result is not independent prediction. Baselines should include a minimal single-surface mechanism, a model without explicit ions, a model without deactivation and a non-AI alternative. An additional mechanism can support a necessity claim only if it improves external prediction without arbitrary extra freedom. If parameters are unidentifiable, obtain more discriminating observations before adding hidden states.

## 6. Acceptance rules for professional calculations

Normal program termination, numerical convergence, correct chemical identity, suitable approximations and agreement with experiment are separate checks. Similar values from different programs may reflect shared approximations; cross-program agreement detects only some implementation or input errors.

First, establish species and conformers. Fix atom mapping, charge, electron count, spin and the numbers of explicit ions and solvent molecules. Check coordination, dissociated species and conformational coverage. Total energies of different compositions cannot define a spin gap. NH deprotonation and tBuOH production require a complete stoichiometric ledger.

Second, validate stationary points. Minima require sufficient convergence and no genuine imaginary modes. A candidate TS requires one imaginary mode corresponding to the target chemical change, adequate optimization and bidirectional IRC or equivalent traceable endpoint evidence. A NEB maximum is not an accepted TS. Low-frequency treatment, standard states and solvation terms must be consistent; an electronic single-point energy is not a free energy.

Third, assess methodology. For representative states and decisive barriers, compare at least two reasonable functionals, basis-set improvement, integration grids and SCF stability. Examine sensitivity to dispersion, explicit ions, solvent model, temperature and conformers. For open shells, inspect S², orbital occupations and state identity. Significant multireference character requires suitable specialist methods or narrower conclusions. An offline r2SCAN-3c/CPCM single-point input is a resource and format pilot, not all of these validations.

Fourth, validate kinetics. Connect barriers through a consistent atom-, charge- and mass-balanced network. Specify initial concentrations, solvent/water/base activities and reversibility. Test detailed balance, identifiability, uncertainty propagation and independent time courses. No physical TOF is released without physical barrier evidence. A MECP additionally requires appropriate treatment of spin-orbit coupling and transition probabilities; its energy is not automatically a conventional transition-state rate.

## 7. The human role in bridging to HPC and professional software

The immediate priorities are a supervisor able to guide computational and experimental validation, an institutional account, and a small public-literature replication project. Supply this report, the actual structure source, the minimal pilot bundle and explicit acceptance questions when applying for access. A request for resources should not rest on the claim that large compute alone will produce a Nature paper. While access is pending, work can continue on Linux, SSH, schedulers, quantum-chemistry input, orbital and frequency analysis, literature and structure checks.

Collect and record non-secret platform information: login instructions, Slurm/PBS scheduler type, group account or project allocation, partition limits, cores and memory per node, wall and concurrency limits, scratch paths and cleanup periods, CPU architecture, ORCA/Gaussian versions and legitimate licensing, module-loading commands and MPI requirements. Passwords, private keys and access tokens must not enter the repository or chat. Institutional authorization, login and multifactor authentication are completed by the user.

The current tools implement a Slurm route. A PBS system requires a separately adapted and verified script. Run the offline validator first, check inputs and hashes, and have the administrator confirm software invocation. Quantum calculations belong on compute nodes, not login nodes. Complete a single-point pilot first; inspect original output, exit status, SCF and identity before the user explicitly enables the remaining array.

The pilot measures whether the software runs, its resource use and evidence retention. Estimate core-hours from measured wall time multiplied by allocated CPUs. Core-hours are not the same as charges; pricing, queue rules and remaining allocations are platform-specific. Do not extrapolate one small single point linearly to all TS, frequency or multireference work. Benchmark representative task classes separately.

Return complete input, stdout/stderr, native program output, job script, necessary geometry/orbital/frequency files, software and dependency versions, job identifier, exit status, resource statistics and a hashed manifest. The collector certifies only fields it can actually check. Missing convergence, identity or completeness evidence remains failed or unknown; values must not be invented.

Detailed commands and the platform-information sheet are in HPC_BRIDGE_ZH.md. The offline bundle includes independent Python validation and collection entry points. This phase has not connected to an HPC platform, submitted a paid job or claimed a local ORCA/Gaussian execution.

## 8. Milestones, stopping rules and publication positioning

| Stage | Required output | Action if missing |
| --- | --- | --- |
| Public benchmark | Experimental structure identity, actual conditions, original sources and reproducible inputs | Suspend large designed-structure scans |
| Algorithm and method validation | Reproducible stationary points, sensitivity analysis and representative resource measurements | Repair inputs or narrow conclusions |
| Distinguishable mechanisms | At least one experimental/computational comparison that excludes competing explanations | Continue diagnostics without announcing a discovery |
| Prospective prediction | Frozen held-out conditions, fair baselines and independent validation | Report failures and applicability limits |
| Broader importance | An explanatory or design principle validated across relevant conditions or systems and addressing an important problem | Prefer an appropriate specialist journal |

These are project gates, not Nature acceptance guarantees. Even if completed, the contribution must be judged for its effect on scientific understanding, overlap with prior work and relevance to a wider readership. Honest negative results and rigorous benchmarks can be valuable. Failures must not be hidden to preserve a journal target.

The final full regression run passed 713 tests, with zero failures, errors or skips. These tests verify covered software behavior and evidence constraints, not chemical mechanisms or journal readiness.

## 9. Reproducibility and responsibility

Code, real calculations and exported documents have separate provenance and hashes. Phase 3 stays frozen; Phase 4 uses new paths. Execution records include versions, threads, memory, time limits and failures. The 405-row post-processing table is explicitly not 405 quantum calculations.

Before submission, qualified researchers must audit structures, methods, data permissions, authorship contributions and every scientific conclusion. AI tools are not authors; people retain responsibility. A manuscript should disclose actual assistance, retain reproducible code/data and plan persistent archiving. A Git repository does not replace a permanent data repository or licensing review. [Nature initial-submission guidance](https://www.nature.com/nature/for-authors/initial-submission); [computational-tools reporting guidance](https://www.nature.com/documents/Computational_tools_reporting_guidelines.pdf).

The current conclusion is that the research direction has been narrowed and concrete verification entry points exist. The central chemical discovery remains to be established. Full effort means prioritizing evidence that can change the assessment, not promising absolute truth or journal acceptance.


## 10. Mathematical and physical foundations: specify what can be proved

The additional CNS-level ambition applies to models, equations, methods and experiments. Statistical thermodynamics, elementary kinetics and identifiability provide established foundations, not inventions of this project. Their role is to expose falsifiable mechanisms. Mathematical properties hold under explicit assumptions; chemical species and parameter validity still require computation and experiment.

Let c be species concentrations in mol/L, alpha and beta the reactant and product stoichiometries, and B the elemental and charge ledger. For a homogeneous, isothermal, constant-volume model:

$$N=\beta-\alpha,\qquad BN=0,\qquad \dot{c}=NJ+u.$$

Net flux J has units mol/(L s); u accounts for addition or removal. Bc is constant only in the closed case. Base addition, sampling and gas exchange require explicit external terms. Free ions, ion pairs, alcohol, water, dormant states and relevant salts/byproducts need composition ledgers. An ion/proton inventory means distribution among chemical forms and their availability, not disappearance of potassium or hydrogen atoms. Precipitation or transport limitation requires an extension beyond this homogeneous core.

Standard chemical potentials and dimensionless activities are defined by:

$$a_i=\gamma_i c_i/c^\circ,\qquad \mu_i=\mu_i^\circ+RT\ln a_i.$$

The new implementation is a restricted ideal-mixture foundation, gamma=1. It does not solve ion association or nonideal activities in toluene. An extension must derive activities consistently from an excess free energy and test that model experimentally. Arbitrarily assigned activity coefficients do not preserve the ideal free-energy dissipation proof. Total base equivalents are not free-base activity.

Each channel uses a shared standard transition-state energy G_TS and prefactor for both directions. This prevents independently fitted forward/reverse constants from violating detailed balance. The baseline transition-state-theory assumption is kappa=1:

$$J_\rho^+=c^\circ\frac{k_{\rm B}T}{h}\exp\left[-\frac{G_{{\rm TS},\rho}^\circ-\alpha_\rho^T\mu^\circ}{RT}\right]\prod_i a_i^{\alpha_{i\rho}}.$$

$$J_\rho^-=c^\circ\frac{k_{\rm B}T}{h}\exp\left[-\frac{G_{{\rm TS},\rho}^\circ-\beta_\rho^T\mu^\circ}{RT}\right]\prod_i a_i^{\beta_{i\rho}},\qquad J_\rho=J_\rho^+-J_\rho^-.$$

Dimensionless activities give flux units mol/(L s) for any molecularity; concentration-power rate constants acquire their corresponding order-dependent units. Species and transition-state energies require consistent temperature, solvent, standard state and composition references. Tunnelling or nonadiabatic corrections must preserve the appropriate microscopic reversibility, rather than change just one direction without justification. [IUPAC transition-state-theory definition](https://goldbook.iupac.org/terms/view/T06470).

In the following ratio, $k_\rho^+$ and $k_\rho^-$ are frequency constants with units s⁻¹: $J_\rho^+/c^\circ=k_\rho^+\prod_i a_i^{\alpha_{i\rho}}$ and $J_\rho^-/c^\circ=k_\rho^-\prod_i a_i^{\beta_{i\rho}}$. Their ratio is dimensionless. For an ideal mixture, concentration-power coefficients instead satisfy $k_{\rho,c}^{\pm}=k_\rho^{\pm}(c^\circ)^{1-m_\rho^{\pm}}$, where $m_\rho^+=\sum_i\alpha_{i\rho}$ and $m_\rho^-=\sum_i\beta_{i\rho}$. Their units are $\mathrm{M}^{1-m_\rho^{\pm}}\mathrm{s}^{-1}$. When forward and reverse molecularities differ, these concentration-power coefficients have different units; their bare ratio is not the dimensionless argument of the logarithm below. Setting the numerical standard concentration to 1 M does not remove this distinction.

For positive concentrations in the closed ideal model:

$$\ln(k_\rho^+/k_\rho^-)=-N_\rho^T\mu^\circ/(RT),\qquad \mathcal{A}_\rho=-N_\rho^T\mu.$$

$$f=\sum_i c_i\{\mu_i^\circ+RT[\ln(c_i/c^\circ)-1]\},\qquad \dot{f}=-\sum_\rho J_\rho\mathcal{A}_\rho\leq0.$$

Here f is free-energy density in kcal/L. The flux ratio obeys ln(J+/J−)=A/(RT), so nonnegative channel entropy production follows from (x−y)ln(x/y) being nonnegative. Shared species potentials enforce cycle consistency. For nonzero u, chemical work mu-transpose times u enters the balance; monotonic decay cannot be demanded for a driven system. These mathematical constraints do not establish a novel mechanism. [Rao and Esposito, nonequilibrium reaction-network thermodynamics](https://doi.org/10.1103/PhysRevX.6.041064).

The implementation separates nonnegative-concentration fluxes from logarithmic-potential diagnostics requiring positive concentrations. At a zero-concentration boundary, mass-action consumption vanishes, giving the continuous model an inward-pointing condition. Numerical integrators still need independent step-size, error and nonnegativity checks. The module is not a complete production predictor and does not turn arbitrary numbers into physical catalytic TOF.

## 11. From quantum surfaces to identifiable inference

Conformer or spin states with the same composition and charge may be coarse-grained through a partition function only when interconversion equilibrates sufficiently rapidly compared with chemistry:

$$G_{\rm eff}=-RT\ln\sum_s g_s\exp[-G_s/(RT)].$$

Do not count g_s twice when the corresponding entropy is already in G_s. This equilibrium effective free energy does not by itself define valid effective kinetics. Interconversion must relax rapidly enough relative to the exit reactions under the conditions being modeled, and the reduced model must reproduce the relevant exit fluxes and timescales of the explicit-state network. Multiple computed energies establish neither that separation nor Boltzmann equilibration during a reactive trajectory. Slow spin conversion requires explicit states and transitions; eliminating it can introduce memory or change the effective rate law. An MECP supplies neither spin-orbit coupling nor nonadiabatic probabilities and rates. Those physical inputs are currently absent.

Measurements require observation equations:

$$y_k(t)=h_k(c(t),\eta)+\epsilon_k(t),\qquad m_k=\mathbb{E}[y_k]=h_k(c(t),\eta),\qquad D_{kj}=\frac{\partial m_k}{\partial\theta_j}.$$

The operator h_k includes overlapping signals, response factors, sampling, instrument dead time and calibration parameters eta; the mean expression assumes zero-mean observation error. The current `observational_identifiability` API receives an already computed mean-observation Jacobian D, positive `parameter_scales` s, and positive `observation_noise_scales` sigma. It implements only diagonal row and column scaling followed by SVD:

$$S=\operatorname{diag}(\sigma)^{-1}D\operatorname{diag}(s).$$

Arbitrary positive scales are unit-aware local parameter scales, not automatically log-parameter derivatives. Only the choice $s_j=\theta_j>0$ at the evaluation point gives $D\operatorname{diag}(s)=\partial m/\partial\log\theta$. For correlated errors, diagonal division is not covariance whitening. A caller must independently establish a positive-definite covariance Sigma and a whitening map W satisfying $W\Sigma W^T=I$, transform the observations, predicted means and Jacobian consistently, and verify that transformation externally. The API can then analyze the supplied whitened Jacobian WD with unit row scales; it neither estimates Sigma nor computes or validates W.

For the declared local coordinates $d\theta=\operatorname{diag}(s)d\xi$, the expression $F_\xi=S^TS$ is the actual mean-parameter Fisher information under a Gaussian observation model with known, parameter-independent covariance, using either independent errors with $\Sigma=\operatorname{diag}(\sigma^2)$ or a correctly whitened correlated model. Without those assumptions it is only a local weighted-least-squares/Gauss–Newton information approximation, not a general Fisher-information identity. Parameter-dependent covariance can contribute additional terms, and a non-Gaussian likelihood requires its own information calculation. The present API reports scaled-Jacobian singular values and local numerical rank; it does not fit a likelihood or establish global structural identifiability, statistical precision, or support from unavailable data.

When two unobserved branches enter one product curve only through k1+k2, more observations of that same signal cannot separate k1 and k2. An additional discriminating species signal or intervention is required, together with structural-equivalence checks, profile likelihood/posterior geometry and external validation. Fisher information or expected information gain can guide design, but actual optimization depends on unavailable data and noise; no fictitious optimal condition is reported.

With temperature, prefactor and activities fixed, TST gives d ln k / d deltaG‡=−1/(RT). At 383.15 K, RT=0.7614 kcal/mol, so a hypothetical 1 kcal/mol barrier error corresponds to a rate factor of about 3.72. This is a mathematical uncertainty example, not a measured error estimate. Functional/basis variation is not automatically a statistical 95% confidence interval. Uncertainty propagation should preserve correlations between shared species energies and keep forward/reverse fluxes thermodynamically consistent.

Physical parameterization still needs validated stationary points, composition-consistent comparisons, conformer/ion-pair sampling, method sensitivity, reaction-progress data and species-resolved observations. Implementing these equations does not convert the older 18-channel solver fixture into an experimentally established mechanism.

## 12. Advanced experiments are selected by discriminating power

The detailed supervised experimental design is in ADVANCED_EXPERIMENTS_ZH.md. Advanced measurements matter when they add independent constraints. First establish calibrated GC/qNMR progress curves, mass balance and alcohol/water/elemental inventories; then determine which competing explanations remain indistinguishable.

| Question | Candidate observation combination | Required boundary |
| --- | --- | --- |
| One species process couples selectivity and dormancy? | Full progress curves, operando carbonyl IR, 31P/1H NMR and controlled recovery comparisons | Correlation is not a common cause; calibrate and intervene independently |
| Proton steps control the observed rate? | O-H/O-D and substrate C-H/C-D comparisons, exchange monitoring and component-reaction controls | A single KIE cannot uniquely identify a TS; pre-equilibria and exchange can confound it |
| Total inventory equals active inventory? | Alcohol quantification, water analysis, elemental totals and species observations | ICP-type totals do not measure free-ion activity; sampling may change speciation |
| Metal structure or spin remains decisive? | Reference-supported operando XAS; EPR when applicable | Control beam damage; negative EPR does not prove absence; do not transfer an Fe-specific route to Mn |

Qualified supervisors and facility specialists must determine feasibility, calibration, time resolution, temperature compatibility, sample size and approved operating procedures. Apply for advanced facilities when their observations can change model selection or prediction uncertainty. Equipment access, beamtime and experiments that have not occurred remain proposed work.


## 13. Additional calculations: public precursor and activation bookkeeping

The request for additional discriminating computation triggered a native check of the public starting species. public_precursor_relaxation_001 started zero quantum calls because the available-memory precheck failed; that record is retained. A separate run, public_precursor_relaxation_002, completed six GFN2-xTB/ALPB(toluene) calls: optimization, a fresh gradient and a complete Hessian for each of the two crystal starting molecules. The run took approximately 120.59 seconds with two threads and sequential execution.

| Crystal starting molecule | Fresh maximum force, eV/Angstrom | Lowest internal frequency, cm^-1 | Complete internal modes | Mapping and minimum checks |
| --- | --- | --- | --- | --- |
| Mn1 | 0.00078029 | 33.4689 | 180; no imaginary modes | Full coordination/covalent graph retained; passed |
| Mn2 | 0.00075099 | 33.8110 | 180; no imaginary modes | Full coordination/covalent graph retained; passed |

Both endpoints individually pass the local-minimum checks on this model. Another AI collaboration agent independently checked 74 native files, gradients and complete Hessians. CIF-label mapping followed by translation and proper rotation gives an all-atom RMSD of only 0.00354 Angstrom: the outputs are geometrically almost the same conformation and must not count as two distinct minima. Charge 0 and uhf=0 are candidate settings, not a determination of the true ground spin. Native outputs, gradients, Hessians, source snapshots and hashes are archived. Mn2's log retains a nonconvergence warning from the automatic thermochemistry symmetry search; that thermochemistry was not adopted, while independent gradient and projected-Hessian checks passed. No reaction free energy at 383.15 K was calculated. The 300 K electronic smearing parameter is not a nuclear temperature, and the built-in toluene parameters have not been validated as an exact description at 110 degrees C. There are no new accepted TSs, MECPs, experimental data or physical TOFs.

The activation ledger separately checks element and charge balance. For neutral P0 containing NH and coordinated Br, removal of only the NH proton leaves an anionic complex if Br is retained. A neutral A0 hypothesis must account explicitly for Br release and the KBr/tBuOH products:

$$P0+KOtBu\rightleftharpoons K^+ + [P0-H]^- + tBuOH.$$

$$P0+KOtBu\rightleftharpoons A0+KBr+tBuOH.$$

These are balanced overall hypotheses, not established elementary steps. Contact ion pairs, aggregation, salt dissolution/precipitation and nonideal activities in toluene remain unresolved. A gas-phase KBr molecule cannot substitute for a salt chemical potential in solution or a solid. activation_hypotheses.json records conservation without invented energies.

A separate synthetic audit evaluated 200 positive concentration states. The largest free-energy dissipation identity defect was 1.42e-14 kcal/(L s); inventory conservation, equilibrium flux, cycle closure and energy-reference shifts were also checked. Abstract A/B/C labels and arbitrary test free energies were used, with no physical network fitting or ODE integration. This checks implementation under stated assumptions and does not establish catalytic chemistry.

COMPUTATION_DECISION_PLAN.md defines G0-G7 prerequisites, discriminating outcomes and stopping rules. Activation identity and key competing pathways precede a small set of higher-level calculations with opposing observable predictions. JACS 2022 already examined base-dependent benzyl-alcohol binding and amido-alkoxide equilibria for Mn-PNP; species effects plus a second scaffold cannot be repackaged as a first discovery. The remaining proposed contribution still requires a shared branching-dormancy mechanism in this borrowing-hydrogen system and added predictive value. [Original study and Figure 7](https://pure.tudelft.nl/ws/portalfiles/portal/120918272/jacs.2c00548.pdf).
