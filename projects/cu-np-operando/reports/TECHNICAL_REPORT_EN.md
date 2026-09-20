# Cu–N–P electrosynthesis under operating conditions: local execution report

Date: 20 September 2026. The user explicitly authorized completing and publishing locally executable work without HPC access or the original supporting information (SI). This report corresponds to the [Chinese version](TECHNICAL_REPORT_ZH.md). The [advanced charter](RESEARCH_CHARTER_ZH.md) describes the longer research programme, not completed results.

## 1. Delivered outcome

This delivery contains a literature and novelty audit, bounded reproducible molecular DFT runs, a neutral reaction/transport numerical kernel, evidence gates and independent code review. **It is an executable local research foundation. It does not establish the Cu electrocatalytic mechanism, a ranking reversal, originality or industrial superiority.**

| Evidence class | Executed result | Supported interpretation |
|---|---|---|
| LITERATURE | 14 source records, 13 novelty-matrix rows and 35 actual search queries; SI unavailable | Research anchor and prior-art constraints; partial G0 |
| CALCULATION | Four launched Psi4 jobs: two converged neutrals and two timed-out anions; two further preflight memory rejections | Local resource measurements and fixed-nuclei neutral electronic energies |
| SOFTWARE VERIFICATION | Reversible neutral surface cycle, film transport, conservative plug flow, separate electrical accounting and archive validation | Implementation and numerical consistency in the tested domain |
| HYPOTHESIS | CuN4/CuN3P1 potential–concentration–transport ranking boundary and product occupation | Questions requiring physical parameters and prospective tests |
| NOT RUN | Target periodic constant-potential paths, physical microkinetics, target AI training, experiments, full TEA/LCA | No corresponding achievement may be inferred |

No complete neutral/anion pair exists. DFT attachment energies, an electron-affinity comparison and a cross-functional charging-energy spread are therefore unavailable. Subtracting neutral absolute energies from different functionals would not repair that missing comparison.

## 2. Research connection and industrial motivation

The anchor is Wang et al., [Phosphorus-Doped Single Atom Copper Catalyst as a Redox Mediator in the Cathodic Reduction of Quinazolinones](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202505085). Its official abstract supports a CuN3P1 microcoordination environment in electroreductive coupling/ring opening and already discusses adsorption/desorption modulation. Attribution should preserve the collaboration involving corresponding authors Hai-Tao Tang, Wenhao Li and Dingsheng Wang.

The proposed question is whether an intrinsic coordination advantage survives product occupation, transport and stability effects under matched structures, measured site densities and calibrated nonaqueous electrode potentials. A validated operating boundary could guide flow-reactor conditions and electrical consumption per isolated product, while revealing failures of a single adsorption descriptor.

This is industrial motivation, not a measured benefit. No real flow-reactor dataset, electrode lifetime, feedstock/electrode price, separation duty or process-safety assessment was supplied. Cost savings, scale-up yield and LCA advantages cannot be reported. The implemented energy function covers full-cell electrical input only, excluding pumps, separation and auxiliaries.

## 3. Literature identity, conflicts and novelty constraints

The official SI endpoints returned HTTP 403. Exact structures and atom/electron mapping for 1a, 2a and 3a remain unverified. A CuN3P1 coordination label does not uniquely specify the support defect, periodic cell or coordinates. Neither the reaction nor catalyst geometry was reconstructed by analogy; coupling/ring opening must not be silently replaced with ordinary hydrogenation.

The Wiley-supplied [article preview](https://www.researchgate.net/publication/390023854_Phosphorus-Doped_Single_Atom_Copper_Catalyst_as_a_Redox_Mediator_in_the_Cathodic_Reduction_of_Quinazolinones) contains caption-level conditions including DMF, 10 mA, 5 h, an undivided cell and the additive string **TDEA**. Earlier records used **TEOA**. The discrepancy remains explicit and neither acronym is used to infer an additive structure or proton/electron balance. See the [source audit](../literature/review_zh.md) and [conflict record](../literature/evidence_conflicts.json). Isolated yield is not faradaic efficiency.

Prior art narrows the possible contribution:

- [BEAST DB](https://doi.org/10.1021/acs.jpcc.4c06826), with its [author preprint](https://arxiv.org/abs/2405.20239), already presents bias-dependent energetic-span rankings. This is not a kinetic ranking for the target organic reaction, but it rules out a broad first-discovery claim about potential changing rankings.
- [Hao et al., JCTC 2025](https://pubs.acs.org/doi/10.1021/acs.jctc.5c00914) already couples constant-potential DFT, microkinetics and Nernst–Planck–Poisson transport. Combining these method classes is not itself original.
- [Oliver et al., ACS Central Science 2025](https://doi.org/10.1021/acscentsci.4c01733) studies interactions between mechanism and transport in organic electrosynthesis. A flow-dependent yield alone does not establish a new scientific principle.

A defensible future contribution would require a specific, reproducible and uncertainty-qualified boundary for source-matched CuN4/CuN3P1 chemistry, then validation under prospectively held-out conditions. If no robust crossing exists, that negative result must be retained. The directed literature and public-patent search is not systematic novelty clearance or a freedom-to-operate analysis. Originality remains unestablished.

## 4. Resources and inherited evidence

The inspected host has an i7-13700H, 14 cores/20 threads, approximately 15.73 GiB RAM and Intel Iris Xe graphics; no CUDA device was detected. Available RAM was much lower than total RAM. Existing Psi4 1.11, xTB 6.7.1 and Python environments were reused without installation. No usable periodic constant-potential engine or HPC allocation was verified.

Final numerical checks used Conda base: Python 3.14.6, NumPy 2.5.2 and SciPy 1.18.1. Direct invocation of another unactivated Conda environment exited inside native linear algebra; configuring its DLL search path restored execution. Failure and recovery logs remain in [review/](review/); an environment crash is not counted as a passing test.

The inherited baseline is commit `6197929f41074cb4b21fb2792c7dcd55772a86c0`. Its 31 tests and hashes for 20 native xTB records/executed sources were rechecked. The earlier GFN1/GFN2 fixed-nuclei charging-energy disagreements remain 0.781648 eV in gas and 0.831980 eV with ALPB-DMF, exceeding the original 0.20 eV gate. Neither the acceptance threshold nor negative results was revised. A molecular preflight does not validate a Cu active site.

## 5. New native DFT calculations

The object is an independently generated 3-methylquinazolin-4(3H)-one model, C9H8N2O, SMILES `CN1C=Nc2ccccc2C1=O`. It inherits the earlier GFN2-xTB neutral minimum. This is not SI confirmation of substrate 1a, a Cu-containing structure or a DFT-optimized geometry. Coordinate SHA256: `29351eee2bde9d27f5d9857c9d24d5ed74faa327a9e0a69ed52117ae34e688b9`.

Protocol: Psi4 1.11; PBE0/B3LYP; def2-SVPD; fixed nuclei in gas; DF-SCF with def2-universal-jkfit; 75/302 grid; energy/density convergence 1e-8. Neutral: 84 electrons, singlet. Anion: 85 electrons, doublet. Each job used two threads, a 512 MiB Psi4 setting and a 300 s wall-time limit. Jobs ran sequentially. RSS thresholds were reduced from 1024 to 896/768 MiB as available memory changed; individual records contain exact limits. Sampling is not an OS-enforced allocation quota.

| Job | Outcome | Wall time / s | Peak RSS / MiB | Electronic energy / Eh |
|---|---|---:|---:|---:|
| DFT001 PBE0 neutral | SCF converged | 108.659 | 612.523 | −531.5663974136239 |
| DFT002 PBE0 anion | Timeout; excluded as a label | 300.061 | 622.555 | unavailable |
| DFT003 B3LYP neutral | Resource preflight rejection; not launched | — | — | unavailable |
| DFT003b B3LYP neutral | SCF converged | 117.232 | 612.547 | −532.1735781682701 |
| DFT004 B3LYP anion | Timeout; excluded as a label | 300.135 | 622.176 | unavailable |

One earlier preflight rejection is separately archived. Thus there were four actual quantum launches and two preflight rejections. A timeout proves neither the absence of an anion nor a defect in the electronic-structure method. Partial SCF energies are not accepted scientific labels. Inputs, native outputs, executed-source snapshots, resource records, successful Molden exports and failure logs are in the [DFT archive](../data/dft_pilot/README.md). Its delivery list covers 58 files, about 6.01 MB, plus the list itself. Eight reconstructible scratch files, about 849.08 MB, remain local with an exclusion/hash inventory.

![Actual native job outcomes](../results/figures/native_job_status.png)

No DFT optimization, frequencies, SCF stability analysis or basis/grid convergence was performed. A single converged solution with low spin contamination does not establish physical accuracy. There is no solvent, electrode, Cu site, reference potential, barrier or thermal free energy. Even a future complete gas-phase fixed-nuclei charge pair would not directly provide a solution reduction potential. No xTB energy is mixed into a DFT difference.

## 6. Numerical model and independent checks

The kernel is strictly the neutral single-site cycle `A + * ⇌ A* ⇌ P* ⇌ P + *`. Forward and reverse Eyring rates share state/transition-state energies; adsorption uses dimensionless activities. Matrix-tree stationary coverages are checked against an independent linear solve. Imported rate constants must also satisfy positive equilibrium constants and logarithmic cycle closure.

Film transport jointly solves `J = km_A(cA_bulk−cA_surface) = km_P(cP_surface−cP_bulk) = Γ·TOF`. The result feeds a conservative first-order implicit axial plug-flow discretization. Independent Radau integration and a closed-form equal-rate solution supply verification. Every parameter is an arbitrary mathematical fixture, not a Cu barrier, training label or reactor prediction.

| Check | Observed result | Interpretation |
|---|---:|---|
| Outlet absolute error with 20/40/80/160 cells | 0.219978 / 0.110789 / 0.0555975 / 0.0278498 mol m⁻³ | First-order spatial convergence |
| Successive refinement error ratios | 1.9855 / 1.9927 / 1.9963 | Approach the expected factor of two |
| 160-cell error / inlet concentration | 2.78498×10⁻⁴ | Fixture discretization error, not chemical prediction error |
| Radau versus analytic outlet | 5.40×10⁻¹³ mol m⁻³ | Independent solver comparison |
| Stationary surface residual | 6.66×10⁻¹⁶ s⁻¹ | Algebraic residual for the tested fixture |
| Film-flux balance residual | 3.39×10⁻²¹ mol m⁻² s⁻¹ | Fixture conservation check |

![Numerical grid verification](../results/figures/neutral_grid_convergence.png)

Integrated surface sources close against net outlet formation to approximately 2.38×10⁻²² and −1.11×10⁻²¹ mol/s for forward and reverse cases. Product already present at the inlet is excluded from net formation. Assumptions include isothermal single-phase operation, fixed volumetric flow, uniform accessible sites and no axial dispersion. There are no potential, current, ionic migration or deactivation equations.

Electrical accounting is independent of the neutral network: `Ecell = ∫Ucell I dt`, normalized by isolated product mass. FE uses formed moles before isolation and externally validated electron stoichiometry. The neutral network cannot supply that stoichiometry. Its flux is not combined with the separate electrical fixture to fabricate electrocatalytic performance.

## 7. Review, acceptance and reproducibility

Three working agents handled literature/novelty, quantum/resources and kinetics/transport, followed by cross-review. Resolved issues included metadata incorrectly unlocking scientific rankings, successful exit status for a rejected gate, malformed data, imported rates bypassing thermodynamic closure and finite-input overflow. DFT archive checks additionally reject nonfinite energies, settings/version conflicts and ambiguous duplicate successes. Process termination and failure recording were repaired. Historical executed scripts remain unchanged; runner repairs apply to future launches.

The final integrated suite contains **42 checks: 41 passed and one explicitly skipped**. The skipped integration check requires an actual complete physical charge pair; no pair was fabricated to satisfy it. The suite comprises 20 kinetics/transport, eight evidence-gate and 14 native-record/process-control checks. See the [validation summary](../results/validation_final/summary.json) and [full log](../results/validation_final/tests.txt). This is neither experimental reproduction nor a proof over the entire floating-point domain.

`physical_ranking_ready` remains false because no validated target-reaction model exists. Complete files, matching hashes or a CALCULATION label can satisfy metadata checks but cannot establish physical validity. All nine scientific inputs remain unaccepted. CI now includes the new suite and frozen-file verification; passing locally does not establish remote CI success.

The [RUNBOOK](../workflows/RUNBOOK.md) gives reproduction and fresh-job instructions. Checks do not launch Psi4. Archived native records are frozen; independent reruns require new output directories. `results/file_manifest.json` binds bytes within this project, excluding caches and itself. It is not an authorship signature or scientific authenticity certificate.

## 8. Charter status and continuation

| Gate | Current status | Missing evidence |
|---|---|---|
| G0 literature/identity/resources | Partial | SI, exact structures/mapping, additive resolution and systematic novelty review |
| G1 benchmark | Not passed | Even the molecular DFT charge pair is incomplete; target Cu/solvent/potential benchmarks are unexecuted |
| G2 mechanism | Not run | Target paths, transition states, connectivity and competing chemistry |
| G3 multiscale model | Engineering components verified only | Neutral fixtures cannot replace calibrated electrochemical chemistry and transport |
| G4 new predictions | Not run | Ranking boundaries under method/parameter uncertainty and held-out validation |
| G5 design benefit | Not run | Matched-constraint productivity, energy, lifetime and process data |
| G6 AI benefit | Not run | Reliable physical labels, strong baselines, fixed budgets and blind testing |
| G7 experiment | Not run | Real independent experiments under preregistered new conditions |

The smallest defensible continuation is: obtain lawful SI and source-matched coordinates; resolve reagent identities and atom/electron balance; configure real compute and a constant-potential engine; converge one Cu/adsorption benchmark; establish a small competing-path set; measure site density and transport; lock uncertainty-qualified prospective predictions; then test them independently. Complex target surrogates should wait for reliable labels; TEA claims should wait for actual cost and lifetime inputs.

The [experimental interface](EXPERIMENT_INTERFACE_ZH_EN.md) specifies conditional product-addition, matched-flow and site-retention measurements. None was performed, and no numerical experimental outcome is claimed. This completes the user-defined local delivery while leaving the scientific programme at the engineering-prototype and documented-negative-result stage.
