# From Designed Pincer Structures to Identifiable Borrowing-Hydrogen Kinetics

## 1. Abstract and frozen campaign status

**Frozen bounded-campaign record.** Real local quantum computation began at 2026-09-12T17:01:17.890189+00:00 and the main search ended at 2026-09-12T20:01:18.364829+00:00, spanning 3.0001 hours including initial gas-phase preparation and the subsequent ALPB-toluene campaign. This is a report of the actual calculations, including unsuccessful searches. The requested complete transition-state, kinetic, and learned-prediction outcomes were not achieved within this run.

Borrowing-hydrogen N-alkylation of aniline with benzyl alcohol was investigated through a factorial library of 24 scaffold-informed designs: Ru(II), Mn(I), Fe(II), and Co(I), three proton-responsive ligand backbones, and Ph/iPr phosphine substitution. The solution-model conformer controller completed 812 native optimization jobs. Its full chemical-identity and donor-coordination screens retained best-found geometries for 70 of 72 intended catalyst/state combinations, including 22 active designs. Both Co bipyridine active designs rearranged through ligand oxygen–carbonyl coupling and are excluded from the intended active-state set. These designed analogues are not a library of experimentally established catalysts, and sampled electronic minima do not establish global minima or electronic ground spins.

Independent full-Hessian auditing selected 70 catalyst-state minimum certificates and 9 molecular reference species. Their consistent 13-temperature records contain 910 catalyst rows and 117 reference rows; balanced dehydrogenation, imine hydrogenation and separated-ion activation reactions yield 858 reaction-temperature rows. A separately certified K–OtBu contact pair and an N-deprotonated hemiaminal cluster provide additional explicitly scoped molecular evidence. Reaction free energies are thermodynamic differences, not activation barriers.

The user-selected conditions are ALPB toluene, 383.15 K baseline, potassium tert-butoxide at 0.05 equivalent, and a requested total-base scan of 0.01–0.20 equivalent. Actual paths were attempted for two Ru, two Mn, and two Co representatives and for explicit condensation clusters. At freeze, 0 dehydrogenation transition-state records passed the complete configured acceptance gates. No complete six-step source-verified reaction network was available, so no real TOF grid, steady-state Campbell degree of rate control, trained surrogate, or Pareto catalyst ranking is claimed. The deliverable includes executable evidence gates, real failed and successful numerical records, endpoint thermodynamics, and an explicit list of the scientific targets still unmet.

## 2. Chemical question and literature context

The net transformation is

$$\mathrm{PhCH_2OH+PhNH_2\longrightarrow PhCH_2NHPh+H_2O}.$$

Its convenient overall stoichiometry hides several mechanistically distinct operations. The alcohol must lose two hydrogen equivalents to become benzaldehyde. Carbon–nitrogen bond formation then connects that carbonyl compound with aniline, commonly through a hemiaminal and an imine. Reduction of the imine returns hydrogen and produces the secondary amine. The catalyst can participate in one or both hydrogen-transfer stages, while condensation may be influenced by water, acid–base chemistry, or the coordination environment. A mechanistic model must account for the coupled sequence rather than identifying a favorable isolated dehydrogenation barrier and treating it as the complete catalytic rate.

Defined manganese PNP complexes have been experimentally demonstrated for N-alkylation of amines with alcohols, establishing a primary-literature basis for considering this reaction family. That evidence supports the research question; it does not establish that every cross-metal analogue assembled here is experimentally accessible or active. [Elangovan et al., Nature Communications 2016, DOI 10.1038/ncomms12641](https://www.nature.com/articles/ncomms12641).

The distinction between inner-sphere and outer-sphere hydrogen transfer must remain a question. In an inner-sphere mechanism, substrate coordination and subsequent hydride migration can be essential. In an outer-sphere mechanism, hydrogen transfer can occur without a persistent substrate–metal bond, often through a ligand-assisted arrangement. Neither a PNP donor set nor a ligand proton alone proves that the outer-sphere pathway dominates. Earlier Co-PNP calculations by Hou and collaborators provide a relevant example in which a non-bifunctional inner-sphere mechanism was favored; that study was affiliated with Sun Yat-sen University. [Hou et al., Dalton Transactions 2015, DOI 10.1039/C5DT02163D](https://pubs.rsc.org/en/content/articlelanding/2015/dt/c5dt02163d/unauth).

Guangxi Normal University work with Cheng Hou as corresponding author has examined the differing roles of proton-responsive sites in borrowing-hydrogen chemistry and has combined DFT with machine learning for related catalyst questions. These publications motivate careful energetic and descriptor definitions, not automatic transfer of their mechanistic assignments or accuracy statistics to the present library. [Mei et al., 2022](https://pubs.rsc.org/en/content/articlelanding/2022/dt/d2dt02597c); [Mo and Hou, 2025](https://pubs.rsc.org/en/content/articlelanding/2025/qo/d5qo01139f). The present project is independent and does not claim endorsement or participation by that research group.

## 3. Definition of the 24-member design space

The design is the Cartesian product $4\times3\times2$. The metal dimension contains Ru(II), Mn(I), Fe(II), and Co(I). The cobalt branch is specifically Co(I); Co(III) is not silently included under the same identifier. The ligand dimension contains an aliphatic MACHO-type PNP scaffold, a pyridine-derived PNN scaffold, and a hydroxylated bipyridine-derived PNN scaffold. Each is constructed with phenyl or isopropyl phosphine substituents. A catalyst identifier encodes these three choices, making it possible to connect every numerical result to a unique proposed molecular graph.

The following table specifies all designs independently of whether an optimization has succeeded. Each metal row summarizes six backbone/substituent combinations. A design count establishes neither literature identity nor successful optimization, and it does not imply exploration of every possible proton-responsive site.

| Metal branch | Backbone identifiers | Substituent identifiers | Number of designs |
| --- | --- | --- | --- |
| Ru(II) | macho_pnp; pyridine_pnn; bipyridine_pnnoh | iPr; Ph | 6 |
| Mn(I) | macho_pnp; pyridine_pnn; bipyridine_pnnoh | iPr; Ph | 6 |
| Fe(II) | macho_pnp; pyridine_pnn; bipyridine_pnnoh | iPr; Ph | 6 |
| Co(I) | macho_pnp; pyridine_pnn; bipyridine_pnnoh | iPr; Ph | 6 |

The bipyridine motif is informed by published ruthenium PNN(O) chemistry featuring a hydroxypyridine unit and a bipyridyl–methylphosphine framework. The exact metal, substituent, charge, and ancillary-ligand combinations generated here remain designed analogues. [de Boer et al., Organometallics 2017, DOI 10.1021/acs.organomet.7b00111](https://pubs.acs.org/doi/10.1021/acs.organomet.7b00111).

Ancillary ligands are explicit model choices. The Ru and Fe active designs contain a carbonyl and a hydride, the Mn designs contain two carbonyls, and the Co designs contain a carbonyl. Together with the chosen deprotonated ligand, these give nominal 16-electron active sets in the generator's counting convention. Formal counts organize hypotheses; they are not measurements of orbital occupancy. Substituting a metal can alter preferred geometry, accessible oxidation states, spin, and ligand lability, so the factorial design is an exploration of model structures rather than a controlled experimental series whose members differ only in one physical property.

## 4. Proton-responsive sites, hydrogen inventory, and charge

The three backbones do not share an identical N–H activation chemistry. In `macho_pnp`, the responsive position is the central amine/amido nitrogen. In `pyridine_pnn`, the represented responsive position is a phosphinomethyl carbon; the deprotonated graph is an aromatic carbanion resonance representation related to the familiar aromatization/dearomatization concept. In `bipyridine_pnnoh`, the modeled position is the hydroxypyridine oxygen. A benzylic carbon provides another possible response site in that scaffold, but the existence of that site does not demonstrate that its competing pathway was searched.

Three states are constructed per design. The active state contains the selected deprotonated ligand. The hydrogenated state equals the active composition plus two hydrogen atoms: one proton is attached to the specified ligand site and one additional hydride is placed at the metal. The protonated reference equals the active composition plus one proton and carries a charge greater by one. Thus

$$\mathrm{Cat_{hydrogenated}=Cat_{active}+H_2},\qquad
\mathrm{Cat_{protonated}^{+}=Cat_{active}+H^{+}}.$$

The equations are compositional bookkeeping, not calculated equilibria. Comparing the absolute electronic energies of these three different compositions cannot establish their populations. Hydrogen uptake requires a hydrogen reference chemical potential. Protonation requires a proton-donor/base equilibrium or another explicitly defined reference. Electron transfer would require an electron chemical potential and a separate charge-state treatment. Omitting these references would convert a formally large electronic-energy difference into an apparently decisive but meaningless stability ranking.

Under the newly specified activation condition, the appropriate model reaction is a proton transfer involving tert-butoxide and the protonated ligand, with tert-butanol as the conjugate-acid product. Potassium may be a spectator, a contact-ion-pair partner, or part of an organized transition structure. A calculation that omits K+ and contains tert-butoxide remains an anionic cluster calculation, and its charge must be preserved in every endpoint and image. Calling that cluster “tBuOK activation” without the counterion qualification would hide an important model approximation in a low-polarity solvent.

A separate potassium tert-butoxide contact-pair reference has now been obtained with ALPB toluene. Its optimized K–O distance is 2.151712546 Å, its tert-butoxide covalent graph is retained, and its full internal Hessian certifies a minimum. The convergence rescue first used a 1000 K electronic-occupation calculation to obtain restart information, then performed a fresh 300 K electronic calculation followed by actual unconstrained optimization and a Hessian at that occupation setting. No high-temperature electronic energy was substituted into the reported reference. Earlier SCC failures remain recorded. This establishes one contact-pair model, not salt solubility or bulk speciation; without a separately defined free-K+ chemical potential, it cannot determine a salt dissociation fraction. Sources: `data/ionic_reference_retry/summary.json` and `data/ionic_reference_retry/distance2.2_warm1000/result.json`.

## 5. Electronic-structure level and the revised solvent conditions

GFN2-xTB supplies real semiempirical quantum energies and analytical gradients. Its parametrized Hamiltonian includes electrostatic and dispersion contributions designed for efficient molecular calculations. This makes it appropriate for a bounded search and for generating hypotheses that can later be examined at higher levels; it does not make all organometallic barriers quantitatively reliable. The original method paper defines the model, not a universal error bar for the present reaction. [Bannwarth, Ehlert, and Grimme, JCTC 2019, DOI 10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176).

Standard GFN1-xTB and GFN2-xTB use spin-independent energy expressions for their restricted open-shell treatment. Consequently, a nominal number of unpaired electrons does not establish that competing spin states have been energetically resolved by an appropriate spin-dependent model. The generator records a nominal low-spin occupation assumption rather than a ground-spin conclusion. This limitation is especially relevant when interpreting Fe and Co structures. [Official xTB spin-polarization documentation](https://xtb-docs.readthedocs.io/en/latest/spgfn.html).

The early observed optimization set used gas-phase GFN2-xTB and a 300 K electronic occupation temperature. Electronic smearing temperature is distinct from the nuclear thermochemical temperature; the two must not be identified. The revised study condition is ALPB toluene, with the user-specified dielectric context approximately 2.38, and a baseline nuclear thermochemical temperature of 383.15 K. The kinetic grid remains 340–440 K in eleven points. Solution-model calculations must be recorded separately from the initial gas-phase preparation set.

For the revised workflow, ALPB `gsolv` output is selected so that no solution standard-state conversion is silently combined with the explicit RRHO conversion. The xTB `reference` and `bar1mol` conventions contain their own reference-state adjustments and are not interchangeable with this choice. Full RRHO thermochemistry referenced to 1 atm is converted to 1 M once, at the actual temperature. [Official xTB solvation reference states](https://xtb-docs.readthedocs.io/en/latest/gbsa.html#reference-states). Solvent model, reference-state mode, charge, and electronic method belong in every energy comparison.

The 383.15 K baseline is not one of the eleven temperatures separated by 10 K between 340 and 440 K. It therefore requires its own thermochemical entry. Re-evaluating partition functions at several temperatures from one accepted geometry is a defined fixed-geometry approximation; it does not imply that geometry, solvent response, or solution speciation were reoptimized at every temperature. If solvent parameters are held fixed across the grid, that assumption should accompany the results rather than being mistaken for a measured temperature-dependent solvent model.

The implemented thermochemical table contains thirteen temperatures: 298.15 K, the 383.15 K baseline, and 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, and 440 K. The accepted geometry, Hessian, and ALPB excess contribution are reused while statistical partition functions and the concentration correction are evaluated at each temperature. Reported enthalpy and entropy are model RRHO/qRRHO quantities assembled on this fixed ALPB potential. No independent solvation enthalpy or derivative of solvent free energy with respect to temperature was calculated. These fields therefore do not constitute a complete enthalpy–entropy decomposition of the real solution.

The bookkeeping is $G_{\rm model}^{\circ}(T)=H_{\rm model}^{\rm RRHO}(T)-TS_{\rm full}^{q}(T)+\Delta G_{\rm conc}(T)$, with the ALPB contribution already included once in the model potential entering $H$. Complete molecular entropy includes translation, rotation, vibration, and the explicitly assumed electronic contribution. Legacy field names `H_298` and `G_298_qRRHO_sol` are retained for compatibility; the row's `temperature` determines its actual evaluation temperature. A 383.15 K row is not a room-temperature calculation because a field name contains “298.” Likewise, the legacy `G_qRRHO_gas` field denotes the value before the concentration shift and may already contain ALPB stabilization. Scientific comparisons use the recorded method, solvent, and temperature rather than inferring physical meaning from those names.

## 6. Execution design, resource limits, and telemetry

The campaign controller started at 2026-09-12 17:01:17.890189 UTC and recorded a three-hour computation window ending at 20:01:17.890189 UTC. Its early configuration used two worker processes with two threads per worker. The time window is an execution budget, not a certificate that a specified number of minima or transition states can be found. Difficult saddle searches can consume much more time than ordinary geometry relaxation, and a resource limit must produce a recorded incomplete outcome rather than a forced acceptance.

An individual job has a defined catalyst, state, starting seed, backend, charge, occupation setting, and output directory. Completion, convergence, and scientific acceptance are separate fields. A backend can terminate without satisfying the force threshold; a force-converged structure can lose its intended coordination; a coordinated stationary point can still have an unwanted negative curvature. These distinctions allow later auditing to determine where progress stopped.

Telemetry should report elapsed time, completed and active jobs, convergence counts, best-energy changes within identical compositions, and worker CPU activity. Energy drift is meaningful only when the compared calculations share their method, solvent, charge, and composition. Reporting a lower energy from a larger molecule as an improvement would be invalid. The ten-minute telemetry interval requested for the campaign is operational monitoring and does not replace the individual native calculation records.

The execution design favors real bounded computation over simulated progress. Retries retain earlier errors and vary specified numerical or geometrical conditions. They do not overwrite an unconverged result with an assumed successful value. When the remaining time is inadequate for a complete Hessian or a pair of downhill connectivity checks, the corresponding status remains incomplete. This prevents a publication table from acquiring a barrier solely because a plot or regression expects one.

## 7. Molecular assembly and conformational exploration

The assembly procedure begins from explicit organic connectivity and identified donor/proton-response atoms. RDKit distance geometry generates candidate ligand coordinates, and constraints help place donor atoms in an intended starting arrangement. Organic force-field preparation is a geometry-building step; its energy is not substituted for GFN2-xTB energy. Metal and ancillary fragments are added with recorded atom mappings before electronic optimization. The original graph, atom order, and hydrogen inventory are retained for subsequent path construction. [RDKit distance-geometry interface](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html).

A successful embedded conformation is not necessarily sterically acceptable. Atom pairs that become unrealistically close during assembly can make the initial electronic problem ill-conditioned or correspond to the wrong bonding pattern. The generator therefore rejects severe overlaps and permits new seeds or torsion choices. Such rejection is an assembly failure, not a finding that a particular catalyst is chemically impossible. Likewise, a donor-distance check after optimization is a geometric retention test, not a complete bonding analysis.

The search can compare multiple seeded or torsion-perturbed starting structures under the same electronic model. The lowest accepted energy among visited structures is the best found conformer. Global-minimum language would require substantially stronger evidence of exhaustive coverage. Repeatedly finding one geometry is useful convergence evidence but does not exclude remote ligand conformations or dissociated states.

CREST provides a literature framework for automated exploration of low-energy conformational space using fast quantum methods. The present report does not call a custom seeded search a completed CREST run merely because xTB is its energy backend. A CREST claim requires actual CREST execution records. [Pracht, Bohle, and Grimme, PCCP 2020, DOI 10.1039/C9CP06869D](https://pubs.rsc.org/en/content/articlehtml/2020/cp/c9cp06869d).

For a validated set of conformers of one species, ensemble populations would satisfy $p_k\propto g_k\exp[-G_k/(RT)]$, and an ensemble free energy can be formed from the corresponding partition sum. Such weighting requires consistent free energies and degeneracies. Using electronic minima alone with arbitrary counting of duplicate geometries can bias populations, particularly when soft conformations differ mainly in entropy. The current best-structure table is therefore a structural-search record, not a complete thermal ensemble.

## 8. Data lineage and levels of evidence

The main early sources are `data/campaign/control.json`, `data/campaign/records.jsonl`, `data/campaign/assembly_failures.jsonl`, and `data/datasets/best_conformers.json`. The control document records the controller settings and declared time window; append-only records describe individual attempts; the best-conformer collection is a derived selection that can change when an improved structure is accepted. A result quoted from that collection must carry its observation time. The table in the next section is a fixed report snapshot and should not be interpreted as a live view.

Native output and optimized coordinates have separate SHA256 values. The first identifies the numerical calculation record; the second identifies the geometry used for descriptors and later calculations. Hashes establish file identity, not chemical correctness. They become useful when a geometry is rerun, a parser changes, or an electronic method is upgraded, because the analyst can determine whether a difference reflects source data or postprocessing.

Five evidence levels are distinguished. An assembled graph establishes a proposed composition. A converged optimization establishes a numerical stationary-search outcome subject to force and structural checks. A full acceptable Hessian establishes local curvature within the chosen model. A transition-state acceptance additionally requires the intended unstable displacement and connection to the specified endpoints. A validated reaction network then requires all relevant state free energies and barriers, consistent references, and an appropriate concentration model. These levels are cumulative but not interchangeable.

The thermochemistry bridge records statuses such as `verified_minimum`, `numerical_uncertainty`, `not_a_minimum`, `failed`, and `deadline_exhausted`. An accepted free energy is emitted only after its stated structural and vibrational criteria are met. A Gibbs-energy column cannot be filled from a force-converged structure alone. Similarly, a NEB maximum without a certified saddle is a path-energy observation, not an accepted activation free energy. This lineage makes absence informative: a missing number identifies a missing scientific prerequisite instead of becoming an unexplained blank.

Identity screening extends to the full ligand graph and ancillary ligands. It checks original nonmetal covalent connectivity, unintended new ligand–ancillary bonds, retained metal–hydride and metal–carbonyl coordination, and the declared proton/hydrogen inventory. The current geometric screen uses 1.28 times the ASE covalent-radius sum for nonmetal bonds, alongside explicitly recorded coordination-distance ranges. These are screening conventions, not universal bond-order thresholds. For intended-catalyst thermochemistry, the accepted minimum, exact accepted-geometry hash, exact NPZ source hash, and full chemical-identity screen must all pass. The explicit downstream selector is `eligible_for_intended_catalyst_thermochemistry` in `data/campaign/thermochemistry/identity_audit.json`. A subsequently improved best conformer does not inherit an earlier geometry's Hessian because their catalyst identifiers coincide.

The retained original Ru-MACHO active-state certificates illustrate why these checks are separate. Their full Hessians can describe mathematical minima even though ligand C–H connectivity has changed. In the retained iPr and Ph source geometries, the affected C–H distances are approximately 2.350 and 2.370 Å, while the corresponding hydrogen lies approximately 1.563 and 1.572 Å from Ru. These are computed hydrogen-migration geometries under the chosen model, not valid realizations of the originally intended active states. They are withheld from that state's thermochemical comparison and preserved as side-reaction evidence. The geometry-specific audit binds the chemical interpretation and Hessian certificate to the same source; ligand-local indices must be mapped to full molecular indices before distances are compared. The separate conformer audit clarification, `data/campaign/neb/ru_conformer_identity_index_clarification.json`, records this mapping and verifies that all fourteen corresponding current geometry hashes match their original campaign records.

A separately diagnosed `Co_bipyridine_pnnoh_Ph` candidate forms a ligand-O to ancillary-carbonyl-C contact. Zero-based atoms 25 and 46 have a separation of 1.496300406 Å and an xTB Wiberg bond order of 0.833712668; the native bond-order file uses one-based indices 26 and 47. Both copied geometry and bond-order files match the SHA256 values in `data/cobalt_carbonyl_rescue/side_reaction_diagnostic.json`. Together these observations support an unintended covalent connection in that calculated structure. They do not demonstrate an experimental side product, an oxidation-state assignment, or a productive catalytic pathway. This candidate is excluded from the intended active-state ranking even if its optimization terminates normally.

The kinetic adapter additionally requires explicit solvent and solvation-reference metadata, reviewed full chemical identity, elemental composition, and molecular charge. It checks atom and charge balance for each elementary step and verifies that each saddle contains the complete declared reactant and shuttle composition. Accepted saddle metadata must also attest reaction-specific motion and endpoint connectivity. These checks consume reviewed calculation evidence; the adapter does not perform a new quantum calculation or infer chemical validity merely from a matching artifact hash. A computed model must retain its complete source snapshot, so declaring an evidence-level string alone cannot replace missing steps.

## 9. Frozen structural results and acceptance boundaries

The following frozen table uses only ALPB-toluene best-found active geometries from `data/datasets/catalyst_descriptors.csv`. There are 22 identity-preserving entries among 24 designs. Angles denote the actual terminal donors: P–M–P for PNP, P–M–N for both PNN backbones. Buried-volume values use 50,000 Monte Carlo samples; per-geometry Wilson intervals and full-precision coordinates are retained in the CSV and `data/structures/`. A dash means unavailable, never zero. Every activation-free-energy cell remains unavailable because no corresponding source-verified free-energy barrier dataset has been produced.

| Metal | Backbone | R | Angle / deg | Vbur / % | ΔG‡ / kcal mol−1 |
| --- | --- | --- | ---: | ---: | --- |
| Ru | PNP | Ph | 162.618 | 61.250 | — |
| Ru | PNP | iPr | 163.554 | 63.568 | — |
| Ru | pyridine PNN | Ph | 160.073 | 64.980 | — |
| Ru | pyridine PNN | iPr | 159.636 | 67.382 | — |
| Ru | bipyridine PNN(O) | Ph | 164.168 | 60.472 | — |
| Ru | bipyridine PNN(O) | iPr | 164.046 | 61.834 | — |
| Mn | PNP | Ph | 163.861 | 63.338 | — |
| Mn | PNP | iPr | 155.406 | 65.342 | — |
| Mn | pyridine PNN | Ph | 162.547 | 65.420 | — |
| Mn | pyridine PNN | iPr | 162.476 | 66.586 | — |
| Mn | bipyridine PNN(O) | Ph | 160.534 | 61.546 | — |
| Mn | bipyridine PNN(O) | iPr | 160.445 | 62.596 | — |
| Fe | PNP | Ph | 167.484 | 65.054 | — |
| Fe | PNP | iPr | 166.351 | 67.998 | — |
| Fe | pyridine PNN | Ph | 167.626 | 67.452 | — |
| Fe | pyridine PNN | iPr | 163.651 | 69.628 | — |
| Fe | bipyridine PNN(O) | Ph | 165.897 | 62.876 | — |
| Fe | bipyridine PNN(O) | iPr | 162.565 | 65.022 | — |
| Co | PNP | Ph | 165.093 | 65.774 | — |
| Co | PNP | iPr | 168.909 | 69.936 | — |
| Co | pyridine PNN | Ph | 169.543 | 66.290 | — |
| Co | pyridine PNN | iPr | 169.238 | 71.884 | — |
| Co | bipyridine PNN(O) | Ph | — | — | — |
| Co | bipyridine PNN(O) | iPr | — | — | — |

The two missing Co bipyridine active entries have a specific computational explanation: an unintended ligand O–carbonyl C connection persists across alternative starting geometries. The Ph case includes an O–C distance of 1.4963 Å and native xTB Wiberg bond order 0.8337 in the saved diagnostic structure. This is evidence of rearrangement on this approximate potential, not proof of a solution reaction yield or experimental decomposition rate. Those structures may possess positive Hessians as mathematical minima while still failing the intended catalyst identity.

Electronic-energy selection and thermochemical certification serve different purposes. The descriptor table reports the lowest sampled eligible electronic geometry, while the thermochemistry CSV selects the lowest certified G at 383.15 K and retains that same certified geometry at all thirteen temperatures. These choices are explicitly linked by geometry hashes and are not merged as if they were necessarily the same conformer. Unexpected Co PNN angles and several iPr/Ph volume reversals therefore remain geometrical observations whose chemical explanation requires additional calculations.

The append-only search records distinguish assembly failures, electronic nonconvergence, donor loss, covalent rearrangement, Hessian rejection, and finite-budget path failure. Current-round controller counters differ from cumulative log counts after a resumed run; the frozen status uses the full historical record. Native evidence archives retain the actual optimization and gradient output, including failed calculations. `data/CAMPAIGN_STATUS.json` and the source-verified thermal export record the final counts without erasing the early gas-phase preparation history.

![Active-state buried volume; hatched cells failed intended identity](../examples/plots/steric_buried_volume_map.svg)

## 10. Steric observables and what they can establish

For donor coordinates $\mathbf r_1,\mathbf r_2$ and metal coordinate $\mathbf r_M$, define $\mathbf a=\mathbf r_1-\mathbf r_M$ and $\mathbf b=\mathbf r_2-\mathbf r_M$. The angle is

$$\theta=\arccos\left[\operatorname{clip}\left(\frac{\mathbf a\cdot\mathbf b}{|\mathbf a||\mathbf b|},-1,1\right)\right].$$

The clipping operation addresses floating-point roundoff, not undefined geometry. Coincident donor and metal coordinates remain invalid. The descriptor is invariant under rigid translation and rotation, but choosing different donor pairs changes its chemical meaning. PNN structures have one phosphorus donor; inventing a P–M–P angle for them would require an atom that is absent. Their terminal P–M–N angle must retain that label in tables, regression features, and any visualizations.

Buried volume measures the fraction of a metal-centered probe sphere occupied by the union of selected ligand atomic spheres:

$$\%V_{\rm bur}=100\frac{\operatorname{Vol}[B_R(M)\cap\cup_j B_{s r_j}(j)]}{4\pi R^3/3}.$$

The default probe radius is 3.5 Å, and Bondi radii are scaled by 1.17, with hydrogens excluded unless explicitly requested. These settings follow the documented SambVca convention, while this implementation uses stochastic integration rather than the server's spatial mesh. [SambVca parameter documentation](https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html).

Uniform sphere sampling uses independent uniform variates, a uniformly distributed direction, and a radial coordinate proportional to $u^{1/3}$. With hit fraction $\hat p$, the estimated standard error in percentage points is $100\sqrt{\hat p(1-\hat p)/N}$. A Wilson interval is recorded because a naive symmetric interval behaves poorly near zero or complete occupancy. This is integration uncertainty conditional on one geometry and one radius convention. It excludes conformational, electronic-structure, radius, and ligand-selection uncertainty.

Descriptors should therefore be compared at matched model settings and interpreted alongside geometry. A small numerical error bar does not imply a precise prediction of accessibility in solution. A conformer ensemble may expose different quadrants even when its scalar buried volumes are similar, and an association transition state can depend on directional pockets invisible to a single aggregate percentage.

Angle type also creates a statistical interpretation issue. A terminal P–M–P angle and a terminal P–M–N angle are both angular measurements, but their ranges arise from different donor identities and backbone constraints. A regression that treats them as one anonymous numerical column may learn the backbone distinction instead of a transferable geometric effect. Retaining donor labels, testing within-backbone comparisons, and examining model performance on held-out ligand groups are necessary to determine what information the feature actually carries.

## 11. Statistical thermodynamics and qRRHO free energies

For an internal mode reported as wavenumber $\tilde\nu_i$ in cm⁻¹, use $\nu_i=100c\tilde\nu_i$ in hertz and define $x_i=h\nu_i/(k_BT)$. The harmonic energy ladder is $\epsilon_n=h\nu_i(n+1/2)$, which gives

$$q_i=\sum_{n=0}^{\infty}e^{-\beta\epsilon_n}=\frac{e^{-x_i/2}}{1-e^{-x_i}},\qquad
U_i=N_Ah\nu_i\left[\frac12+\frac{1}{e^{x_i}-1}\right].$$

With $S_i=R[\ln q_i+T\partial_T\ln q_i]$, the zero-point contribution cancels from the entropy expression:

$$S_i^{\rm HO}=R\left[\frac{x_i}{e^{x_i}-1}-\ln(1-e^{-x_i})\right].$$

As frequency tends to zero, $S_i^{\rm HO}\sim R(1-\ln x_i)$ diverges. Flexible ligand motions can therefore dominate a harmonic entropy even though they are better understood as hindered internal rearrangements. Grimme's interpolation introduces a free-rotor limit to reduce this sensitivity. It is a controlled modeling convention, not an exact solution of coupled torsional quantum dynamics. [Grimme, Chemistry—A European Journal 2012, DOI 10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497).

For the rotor surrogate,

$$\mu_i=\frac{h}{8\pi^2\nu_i},\quad
\mu_i'=\frac{\mu_iB_{\rm av}}{\mu_i+B_{\rm av}},\quad B_{\rm av}=10^{-44}\ \mathrm{kg\,m^2},$$
$$q_i^{\rm FR}=\left(\frac{8\pi^3\mu_i'k_BT}{h^2}\right)^{1/2},\qquad
S_i^{\rm FR}=R\left(\frac12+\ln q_i^{\rm FR}\right).$$

The reduced inertia is essential to the finite low-frequency limit. The quantity $h/(8\pi^2\nu_i)$ alone is unreduced and would not reproduce that limit. The default weighting function and corrected vibrational entropy are

$$w_i=[1+(\tilde\nu_0/\tilde\nu_i)^4]^{-1},\quad\tilde\nu_0=100\ \mathrm{cm^{-1}},\qquad
S_{\rm vib}^{q}=\sum_i[w_iS_i^{\rm HO}+(1-w_i)S_i^{\rm FR}].$$

Only the vibrational entropy is replaced. The complete molecular entropy satisfies

$$S_{\rm total}^{q}=S_{\rm total}^{\rm RRHO}-S_{\rm vib}^{\rm HO}+S_{\rm vib}^{q},\qquad
G^q=H_{\rm RRHO}-TS_{\rm total}^{q}.$$

The supplied total enthalpy already contains electronic energy, zero-point energy, and thermal contributions. Adding electronic energy or ZPVE again would double count them. If energies are in kcal mol⁻¹ and entropy in cal mol⁻¹ K⁻¹, the product $TS$ is divided by 1000. An accepted minimum contributes all positive internal modes; an accepted first-order saddle excludes its single unstable reaction coordinate. Negative modes are never made positive by taking absolute values. Parser-aware treatment also prevents applying qRRHO twice to an output that already includes it.

## 12. Standard states, activation references, and free-energy consistency

The pressure and concentration reference states differ even for an ideal gas. Starting from $\mu(p)=\mu(p^\circ)+RT\ln(p/p^\circ)$ and $C_g=p/(RT)$ gives

$$\Delta G_{\rm conc}=RT\ln\left(\frac{C^\circ RT}{p^\circ}\right).$$

The logarithm must have a dimensionless argument. In SI units, 1 M is 1000 mol m⁻³, while 1 atm is 101325 Pa. At 298.15 K the correction is 1.8943284455 kcal mol⁻¹. At 383.15 K and at each grid temperature, it is evaluated again; a rounded room-temperature correction is not reused. The revised ALPB `gsolv` convention and this conversion together define the chosen bookkeeping, with the concentration shift added once.

For a reaction, the shift is multiplied by the change in molecularity, $\Delta n=\sum_j\zeta_j$, where products have positive stoichiometric coefficients. Association of two molecules into one complex therefore receives $-\Delta G_{\rm conc}$ as a reaction correction. An activation free energy referenced to two separated reactants has the same counting issue. A transition state referenced to a preassociated complex has a different standard activation free energy and must not be substituted into a bimolecular Eyring expression without redefining the step.

Chemical references also matter for activation. Comparing a deprotonated catalyst with its protonated precursor requires the tert-butoxide/tert-butanol pair or another complete proton-transfer reaction. Comparing an active catalyst with its hydrogenated state requires a complete hydrogen inventory. Atoms, charge, and reference conditions must balance before a numerical difference becomes an interpretable reaction energy.

ALPB represents a continuum contribution; it does not determine the concentration of free tert-butoxide, the structure of potassium ion pairs, or the number of bound solvent molecules. A mixed table of gas-phase minima and ALPB transition states would confound the model difference with the barrier. The acceptance gate therefore requires a consistent calculation protocol for all quantities in a reaction-network snapshot. Numerical convergence and unit consistency are necessary for this comparison, but solvent and speciation adequacy remain separate scientific questions.

Independent molecular references provide a useful thermodynamic result before a complete catalyst network is available. Define $A$ as benzyl alcohol, $D$ as benzaldehyde, $N$ as aniline, $Q$ as the neutral hemiaminal, $I$ as the imine, $W$ as water, and $P$ as N-benzylaniline. These symbols retain the same meanings in the network below. Standard reaction free energies are calculated as $\Delta_rG^{\circ}=\sum_i\nu_iG_i^{\circ}$, with positive product coefficients and negative reactant coefficients. The three evaluated reference reactions are $D+N\rightarrow Q$, $Q\rightarrow I+W$, and $A+N\rightarrow P+W$.

| Reaction of independently optimized references | $\Delta n$ | $\Delta_rG^{\circ}$ at 383.15 K / kcal mol⁻¹ |
| --- | ---: | ---: |
| Hemiaminal formation, $D+N\rightarrow Q$ | −1 | −7.51178867 |
| Hemiaminal dehydration, $Q\rightarrow I+W$ | +1 | +5.52597660 |
| Net borrowing-hydrogen substitution, $A+N\rightarrow P+W$ | 0 | −4.23531646 |

All thirteen temperatures are exported in `data/datasets/reference_reaction_thermodynamics.csv`. The accompanying `reference_reaction_thermodynamics.manifest.json` records seven source species, stoichiometric coefficients, atom and charge balances, result and geometry hashes, NPZ hashes, accepted spectra, metadata-correction histories, and the reproducible export script. Each source passed a fresh comparison of its molecular graph, accepted geometry, complete all-positive internal spectrum, and every temperature row against its saved Hessian. The source-hash correction associates each reported hash with the NPZ actually named by that source path; the independent original raw-Hessian checksum remains available. This is verified postprocessing of actual GFN2-xTB calculations, not a new quantum calculation or a synthetic dataset.

At the chosen standard states, the model favors hemiaminal formation and gives positive dehydration free energy, while the net substitution remains exergonic. The sum for $D+N\rightarrow I+W$ is −1.98581207 kcal mol⁻¹ at 383.15 K. These are differences between stable molecular-reference free energies, not activation barriers. Their signs cannot identify a rate-determining step, generate a rate constant, or establish catalyst activity. Actual reaction free energy also includes the activity correction:

$$\Delta_r G=\Delta_r G^\circ+RT\ln\prod_i a_i^{\nu_i}.$$

Water accumulation, substrate depletion, association, and ionic speciation therefore matter when translating reference thermodynamics into a reaction mixture. Catalyst and shuttle species cancel from a net balance only after a chemically closed mechanism has been established, and this cancellation does not determine kinetic order in base.

### Certified catalyst reaction thermodynamics

All table entries are kcal/mol at 383.15 K, for separate 1 M molecular species and the fixed ALPB potential. Alcohol dehydrogenation denotes Cat + benzyl alcohol → Cat–H₂ + benzaldehyde; activation denotes protonated Cat(+) + tBuO(−) → Cat + tBuOH. Each sum is atom- and charge-balanced and uses the complete certificate attached to its own geometry.

Both Ru PNP alcohol-dehydrogenation reactions are exergonic in this model, whereas both Mn PNP alcohol-dehydrogenation reactions are endergonic. The imine-hydrogenation endpoint reactions show the opposite sign pattern for these same four designs. This supports a specific thermodynamic distinction between these chosen molecular states. It does not establish a lower Ru barrier or faster Ru turnover: the association basins, true saddles, imine reduction, condensation and catalyst populations have not been fully determined. Several pyridine-PNN reactions are more favorable than their PNP counterparts; replacing the responsive atom and ligand topology also changes the molecular state, so this cannot be assigned solely to metal or substituent electronics.

The strongly negative separated-ion proton-transfer values describe neutralization of deliberately isolated charged references in a low-dielectric model. They must not be interpreted as measured pKa values, a dissolved tBuOK activation equilibrium, or evidence that 0.05 equivalent total salt generates a specified free-anion concentration. Contact pairing and the counterion chemical potential are absent from that reaction definition. The complete 13-temperature series is available for thermal-model sensitivity, while the omitted solvent temperature derivatives remain a distinct limitation.

Imine hydrogenation uses Cat–H₂ + imine → Cat + product amine. Adding this reaction to alcohol dehydrogenation and the two certified condensation reference reactions cancels the catalyst states and intermediates exactly. The independently exported cycle residual checks this thermodynamic bookkeeping at every temperature. Favorable hydrogen uptake can be accompanied by less favorable hydrogen delivery; neither isolated endpoint difference establishes a rate or selects a catalyst.

| Design | Alcohol dehydrogenation ΔG | Imine hydrogenation ΔG | Separated-ion deprotonation ΔG |
| --- | ---: | ---: | ---: |
| Ru / PNP / Ph | -10.775 | +8.526 | -62.348 |
| Ru / PNP / iPr | -11.398 | +9.148 | -70.884 |
| Ru / pyridine PNN / Ph | -18.043 | +15.793 | -44.000 |
| Ru / pyridine PNN / iPr | -20.092 | +17.842 | -46.760 |
| Ru / bipyridine PNN(O) / Ph | -21.305 | +19.055 | -55.810 |
| Ru / bipyridine PNN(O) / iPr | -21.393 | +19.144 | -55.049 |
| Mn / PNP / Ph | +2.589 | -4.839 | -87.579 |
| Mn / PNP / iPr | +6.007 | -8.257 | -85.690 |
| Mn / pyridine PNN / Ph | -29.891 | +27.641 | -42.664 |
| Mn / pyridine PNN / iPr | -34.350 | +32.100 | -41.239 |
| Mn / bipyridine PNN(O) / Ph | -15.129 | +12.879 | -71.003 |
| Mn / bipyridine PNN(O) / iPr | -13.874 | +11.625 | -69.603 |
| Fe / PNP / Ph | -2.628 | +0.378 | -66.410 |
| Fe / PNP / iPr | -3.655 | +1.405 | -60.589 |
| Fe / pyridine PNN / Ph | -15.663 | +13.414 | -31.458 |
| Fe / pyridine PNN / iPr | -18.785 | +16.536 | -35.108 |
| Fe / bipyridine PNN(O) / Ph | -7.145 | +4.896 | -62.274 |
| Fe / bipyridine PNN(O) / iPr | -6.790 | +4.541 | -59.328 |
| Co / PNP / Ph | +1.705 | -3.955 | -48.902 |
| Co / PNP / iPr | +10.376 | -12.625 | -46.660 |
| Co / pyridine PNN / Ph | -1.227 | -1.023 | -43.555 |
| Co / pyridine PNN / iPr | +4.185 | -6.434 | -41.237 |


![Certified endpoint free energies; no transition-state barrier](../examples/plots/pes_comparison_ru_vs_mn.svg)

## 13. Reaction-path construction and climbing-image NEB

A path search requires chemically compatible endpoints with the same atom ordering, total composition, and charge. For alcohol dehydrogenation, the hydrogen originally bonded to oxygen must map to the selected ligand response site, while the appropriate carbon-bound hydrogen maps toward the metal. The carbonyl product must retain the benzyl carbon and oxygen identities. Adding or deleting a hydrogen between endpoints would change the problem rather than solve the intended reaction.

The implemented path protocol uses seven internal images between two endpoints, yielding nine images in total. This distinction avoids interpreting “seven intermediate images” as seven total structures. Initial interpolation can be improved to reduce atomic overlap, but its coordinates remain guesses. Each image must obtain energy and forces from the actual chosen quantum backend. A curve drawn through interpolated coordinates is not a potential-energy profile unless those image energies were evaluated.

For image $i$ with tangent $\hat\tau_i$, the ordinary NEB force combines a perpendicular physical force with a parallel spring force:

$$\mathbf F_i^{\rm NEB}=-\nabla E(\mathbf R_i)+[\nabla E(\mathbf R_i)\cdot\hat\tau_i]\hat\tau_i
+k\left(|\mathbf R_{i+1}-\mathbf R_i|-|\mathbf R_i-\mathbf R_{i-1}|\right)\hat\tau_i.$$

The projection discourages artificial shortening of the path across curved valleys, while springs distribute images. After suitable path relaxation, the selected climbing image removes its spring force and reverses the physical force along the tangent:

$$\mathbf F_i^{\rm CI}=-\nabla E(\mathbf R_i)+2[\nabla E(\mathbf R_i)\cdot\hat\tau_i]\hat\tau_i.$$

The method targets a saddle along the represented path, but it neither proves that the globally lowest pathway has been found nor establishes the saddle's vibrational order. [Henkelman, Uberuaga, and Jónsson, JCP 2000, DOI 10.1063/1.1329672](https://doi.org/10.1063/1.1329672).

Endpoint optimization, chemical identity, band convergence, and the candidate's unmodified force norm are recorded separately. A geometry from an unconverged band may be used as a dimer initial guess, with that provenance retained explicitly. This is an initialization strategy, not a waiver of transition-state certification: the refined candidate must independently meet the true-force criterion, full-Hessian mode count, reaction-specific displacement criterion, and correct two-sided downhill connectivity before an accepted barrier can be formed. A certified saddle obtained by refinement does not retroactively make the original band converged. The highest energy image of an unconverged band remains a path-search observation and is never reported directly as a barrier. Step limits, electronic failures, and exhausted deadlines remain explicit outcomes.

An actual Mn-MACHO-PNP-Ph candidate provides a geometry-specific diagnostic of these acceptance boundaries. The record completed at 2026-09-12 19:13:22 UTC used frame 110 of a preceding dimer trajectory and retained the screened catalyst identity. Its complete Cartesian Hessian is 240 × 240 for 80 atoms; projection of six external directions leaves 234 internal modes. Four signed curvature frequencies are negative: −378.9759, −27.7247, −18.9928, and −7.7317 cm⁻¹. The independent true maximum force is 0.2800261 eV Å⁻¹. This is therefore a nonstationary curvature diagnostic, not an accepted first-order saddle, and the number of negative directions must not be presented as the index of a stationary saddle. The exact candidate and NPZ hashes, full spectrum, force, and source trajectory are retained in `data/campaign/neb/Mn_macho_pnp_Ph_diagnostic_1789240138_30924/diagnostic.json`.

The leading negative mode has proton-transfer, hydride-transfer, and combined overlaps of 0.1060991, 0.8215832, and 0.6559704, respectively. Although its direction is marked concerted and its frequency lies inside the dehydrogenation window, the proton overlap is below the separate 0.15 threshold; its large hydride component cannot replace that requirement. Together with the substantial force and additional negative modes, this prevents TS certification and Gibbs-barrier construction. By the final freeze, dimer continuation initialized from this diagnosed mode had ended without TS certification. Its final full-Hessian diagnostic appears in the table below, alongside the original candidate; each retains its own force, mode and identity records. These unsuccessful bounded searches do not establish that the Mn pathway is impossible.

Two explicit numerical sources can supply a complete Hessian. The ASE-driven `full_hessian` in `src/pincer_catmech/kinetics/neb_ts_search.py` forms every Cartesian column from central differences of actual quantum forces, requiring 6N force evaluations. The native provider in `scripts/run_neb_campaign.py` instead calls the full unbiased xTB Hessian implementation in `src/pincer_catmech/quantum/thermochemistry.py`; it requires normal termination, zero frozen atoms, scale 1.0, and a finite unprojected 3N × 3N matrix. Both are full numerical Hessians, with their displacement conventions recorded, rather than selected-mode approximations. For this candidate, native and independent-gradient electronic energies differ by −3.2871 × 10⁻⁸ eV, within the 0.001 eV consistency gate. `native_hessian_provider.json` in the same diagnostic directory preserves that comparison and hashes of the native input, command, output, and raw matrix. The saved geometry, NPZ, and five native artifacts were independently hash checked during report preparation. This validates provenance and surface consistency; it does not assert element-by-element agreement between two separately evaluated Hessian schemes or override stationary-point acceptance.

### Actual nonstationary curvature diagnostics

These fixed-candidate full-Hessian calculations are diagnostics at nonstationary geometries. Columns report unmodified maximum force (eV/Å), negative-mode count, and the most negative mode's frequency (cm⁻¹) and separate proton/hydride overlaps. They are not stationary-point certificates or activation free energies. Later mode-guided searches retain separate results in the final NEB inventory.

| Candidate | Force | Negative modes | Lowest frequency | Proton overlap | Hydride overlap |
| --- | ---: | ---: | ---: | ---: | ---: |
| Mn PNP Ph / 19:13:22 UTC | 0.2800 | 4 | -378.976 | 0.106 | 0.822 |
| Mn PNP Ph / 19:56:03 UTC | 0.2220 | 2 | -338.863 | 0.117 | 0.822 |
| Ru PNP iPr / 19:44:30 UTC | 0.3869 | 5 | -778.913 | 0.823 | 0.040 |
| Ru PNP Ph / 19:34:30 UTC | 0.5983 | 4 | -373.320 | 0.046 | 0.084 |

## 14. Saddle refinement, Hessians, and corrective searches

Dimer refinement uses two nearby configurations to estimate the direction of negative curvature through force differences. Rotation seeks a low-curvature direction, and translation reverses the force component along that direction while relaxing orthogonal components. The procedure can refine a candidate without building a full Hessian at every step, but final curvature certification still requires additional analysis. [Henkelman and Jónsson, JCP 1999, DOI 10.1063/1.480097](https://henkelmanlab.org/pubs/henkelman99_7010.pdf).

For Cartesian component indices $a,b$, a central force difference gives

$$H_{ab}=\frac{\partial^2E}{\partial x_a\partial x_b}≈-\frac{F_a(\mathbf x+\delta\mathbf e_b)-F_a(\mathbf x-\delta\mathbf e_b)}{2\delta}.$$

The symmetric Hessian is mass weighted as $\widetilde H=M^{-1/2}HM^{-1/2}$. Translational and rotational directions are projected before diagonalization, so genuine internal negative modes are not confused with numerical rigid-body artifacts. With consistent SI conversion, the signed wavenumber associated with eigenvalue $\lambda_k$ is $\operatorname{sgn}(\lambda_k)\sqrt{|\lambda_k|}/(2\pi c)$. The antisymmetric component before symmetrization provides a numerical-consistency diagnostic; symmetrization alone must not hide unreliable force differences.

The requested dehydrogenation criterion requires exactly one internal imaginary frequency in the interval −1800 to −300 cm⁻¹. In addition, its displacement must overlap both the proton-transfer and hydride-transfer coordinates in the intended concerted direction. A single imaginary torsion would satisfy the mode count but fail reaction identity. Conversely, a physically meaningful saddle outside the requested frequency interval is not accepted under this campaign's declared rule, although it may deserve separate scientific examination.

Corrective searches displace structures along unwanted internal modes and repeat refinement within explicit attempt and time bounds. The histories retain displacement direction, amplitude, forces, curvature, and outcome. Following the candidate in both signs of the reaction mode and relaxing downhill tests connection to the specified endpoint basins. A successful mode count without appropriate descent connectivity remains insufficient. The actual free-energy barrier is then constructed from thermochemistry of the certified saddle and its correctly chosen reactant reference, not from the numerical height of an unverified NEB image.

Condensation requires its own reaction-mode definition. The thermochemistry bridge accepts an explicit `ts_imaginary_window_cm1=None`, which disables only the dehydrogenation-specific frequency window. It still requires exactly one negative internal frequency and no zero internal modes. Depending on the explicitly modeled condensation step, the unstable displacement must describe the relevant proton transfer or transfers together with C–O cleavage and C–N contraction, followed by correct downhill connection. No hydride motion is imposed on a mechanism that contains none. The complete signed spectrum and chosen window convention remain in the saved evidence.

## 15. Elementary network and the eight-equation formulation

The modeled species are alcohol $A$, free catalyst $C$, alcohol complex $X$, aldehyde $D$, hydrogenated catalyst $H$, aniline $N$, hemiaminal $Q$, imine $I$, water $W$, imine–hydrogenated-catalyst complex $Y$, and product $P$. Six reversible steps are defined:

$$A+C\rightleftharpoons X,\quad X\rightleftharpoons D+H,\quad
D+N\rightleftharpoons Q,\quad Q\rightleftharpoons I+W,\quad
I+H\rightleftharpoons Y,\quad Y\rightleftharpoons P+C.$$

Writing every step reversibly ensures that reverse rates are derived from the same state free energies rather than implicitly set to zero. The stoichiometric matrix $S$ has species as rows and the six steps as columns. With species ordered $(A,C,X,D,H,N,Q,I,W,Y,P)$,

$$S=\begin{pmatrix}
-1&0&0&0&0&0\\-1&0&0&0&0&1\\1&-1&0&0&0&0\\
0&1&-1&0&0&0\\0&1&0&0&-1&0\\0&0&-1&0&0&0\\
0&0&1&-1&0&0\\0&0&0&1&-1&0\\0&0&0&1&0&0\\
0&0&0&0&1&-1\\0&0&0&0&0&1
\end{pmatrix},\qquad\dot{\mathbf c}=S\mathbf r.$$

Three useful conserved inventories are

$$C_T=C+X+H+Y,$$
$$B_T=A+X+D+Q+I+Y+P,\qquad N_T=N+Q+I+Y+P.$$

The second is the benzyl-derived organic fragment, and the third is the aniline-derived fragment. Total carbon associated with the reacting substrates is $7B_T+6N_T$, while substrate nitrogen is $N_T$. Ligand carbon and nitrogen remain conserved through $C_T$. This makes elemental interpretation explicit instead of calling every species concentration a carbon balance.

Eliminating $C$, $A$, and $N$ through these inventories leaves the requested eight dynamic concentrations $(X,D,H,Q,I,W,Y,P)$. Their derivatives are

$$\dot X=r_1-r_2,\quad\dot D=r_2-r_3,\quad\dot H=r_2-r_5,\quad\dot Q=r_3-r_4,$$
$$\dot I=r_4-r_5,\quad\dot W=r_4,\quad\dot Y=r_5-r_6,\quad\dot P=r_6.$$

Other conservation relations can exist in a closed network, so eight coordinates need not be the mathematically minimal dimension. This representation satisfies the requested interface while reconstructing all eleven concentrations and exactly preserving the three explicit inventories. It does not replace missing elementary chemistry; all six step barriers are still required before a research simulation is allowed.

The bounded condensation investigation demonstrates a limitation of the declared six-step network. Four actual ALPB cluster searches were completed. The neutral 43-atom hemiaminal/tBuOH model retained both endpoint identities, but its nine-image NEB did not converge in 180 steps; its best recorded post-warmup projected force was 0.130327 eV Å⁻¹, above the 0.07 eV Å⁻¹ threshold. No saddle or barrier was accepted from that band. In the 57-atom anionic model, unconstrained optimization transferred N–H to the originally mapped tert-butoxide oxygen, so the intended neutral-hemiaminal reactant basin was rejected. This is an observed optimization outcome, not proof of a barrierless solution reaction.

The resulting N-deprotonated hemiaminal cluster was then studied as a distinct chemical model, $Q_N^-+2BH$, where $B^-\equiv\mathrm{tBuO^-}$ and $BH\equiv\mathrm{tBuOH}$. Two geometry-specific full-Hessian records for this one model contain 165 positive internal modes, with lowest frequencies 23.6346 and 24.6279 cm⁻¹. They are repeated numerical characterizations of one chemical identity, not two independent catalysts or a population estimate. Standalone accepted-minimum records and failed parent searches are distinguished in `data/condensation_search/summary.json` and the `endpoint_0_verified_minimum.json` files. In two product-preparation attempts, including one assembled from separately accepted imine, water, alcohol, and alkoxide references, unconstrained optimization reformed C–O connectivity and returned the proton to the original alcohol oxygen. Both product identities were rejected, and no condensation activation free energy or TOF was generated.

The compositional extension suggested by these calculations is $Q+B^-\rightleftharpoons Q_N^-+BH$, followed by a separately examined elimination balance $Q_N^-+BH\rightleftharpoons I+W+B^-$. The actual cluster implementation includes an additional mapped $BH$ on both sides of the second balance. These equations define an explicit acid–base hypothesis and conserve atoms and charge; they do not assert that either balance has been certified as one elementary step. Adding $Q_N^-$, explicit base/conjugate-acid species, or their association complexes requires a revised stoichiometric matrix, concentration interpretation, and conservation audit. The new intermediate cannot be silently inserted into the original eight-equation model while retaining its six original barriers. Potassium remains absent from these anionic path models, even though an independent contact-pair reference is now available.

## 16. Eyring molecularity, detailed balance, and the base channel

For a step with reactant molecularity $m$, standard concentration $C^\circ$, and standard activation free energy $\Delta G^{\ddagger\circ}$, transition-state theory gives

$$k_f=\frac{k_BT}{h}\exp\left[-\frac{\Delta G_f^{\ddagger\circ}}{RT}\right](C^\circ)^{1-m}.$$

Thus a unimolecular rate constant has units s⁻¹, while a bimolecular rate constant has units M⁻¹ s⁻¹. Equivalently, the forward flux is

$$r_f=C^\circ(k_BT/h)e^{-\Delta G_f^{\ddagger\circ}/RT}\prod_i(c_i/C^\circ)^{\alpha_i}.$$

The dimensionless activity form exposes the standard-state factors that disappear numerically when $C^\circ=1$ M but remain essential to the units.

One absolute transition-state free energy defines both barriers: $\Delta G_f^\ddagger=G_{TS}-\sum_i\alpha_iG_i$ and $\Delta G_r^\ddagger=G_{TS}-\sum_i\beta_iG_i$. Their difference is the standard reaction free energy. Consequently, $k_f/k_r$ contains the corresponding equilibrium factor and molecularity units. Independently assigning forward and reverse barriers would break this relationship. Negative apparent activation free energies require scientific review or a separate capture treatment; the Eyring gate does not silently reinterpret them as validated barrierless association.

The revised base chemistry involves tert-butoxide and tert-butanol as a proton shuttle during hemiaminal dehydration. If both are present as independent spectator species on both sides of one elementary event, their standard free energies enter both reference sums, and their activities multiply both forward and reverse fluxes. This creates concentration dependence only after an actual transition structure for that molecularity is supplied. A label saying “base assisted” is not enough to justify an arbitrary first-order multiplier.

The total-base grid is 0.01, 0.05, 0.10, and 0.20 equivalents. Total added tBuOK is not automatically free tert-butoxide. The grid API therefore requires an explicit free-base fraction or speciation assumption, and a two-species shuttle requires a separately stated tert-butanol concentration. Potassium pairing and aggregation are unresolved model questions until supported by calculation or experiment. Without the relevant path free energy and concentration interpretation, the base-dependent kinetic grid remains unavailable.

The present shuttle approximation holds the specified free-base and conjugate-acid activities fixed during a batch trajectory. Its exported condition record preserves total base, free-base fraction, actual free-base concentration, and conjugate-acid concentration separately. This defines a conditional kinetic model with fixed shuttle reservoirs; it does not solve potassium pairing, acid–base redistribution, activation, or aggregation equilibria. A changing free-base fraction requires an extended speciation model and corresponding balances. A nonzero base input without a computed base-assisted path is rejected rather than silently ignored.

## 17. Stiff integration, observables, and numerical verification

The difference between rapid association and slower chemical conversion can create widely separated time scales. The implementation therefore uses the implicit Radau or BDF methods in `scipy.integrate.solve_ivp`. Relative and absolute tolerances are explicit, and failure of the integrator is distinct from a scientifically invalid input dataset. [SciPy `solve_ivp` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).

The solver evolves the eight reduced variables and reconstructs the eliminated concentrations at every evaluation. Small negative trial values can arise during implicit iterations; nonnegative concentrations are used in elementary flux evaluations, while the completed trajectory is checked for materially negative or nonfinite values. This numerical accommodation is bounded by explicit diagnostics and is not permission to repair a failed integration by clipping a large negative concentration after the fact.

For a batch duration $t_f$ and total catalyst $C_T$, the average net turnover frequency is

$$\mathrm{TOF}_{avg}=\frac{P(t_f)-P(0)}{C_Tt_f},\qquad
\mathrm{TOF}_{inst}(t_f)=\frac{\dot P(t_f)}{C_T}.$$

Both have units s⁻¹, but they answer different questions. A high initial rate can coexist with a lower average rate after substrate depletion or water accumulation. A final product concentration divided by catalyst concentration is a turnover number, not a turnover frequency unless divided by a stated time. Conditions must include initial substrate concentrations, loading, base interpretation, temperature, and duration before catalysts are compared.

Analytical tests verify the eight-equation machinery using a declared first-order release problem whose exact concentration is proportional to $1-e^{-kt}$. Additional tests check conservation, rate units, reversible ratios, zero-flux sensitivity, evidence rejection, and both stiff solvers. Those tests are software validation; they do not supply missing chemical barriers. The observation that an ODE implementation reproduces an analytical curve cannot be represented as a computed borrowing-hydrogen profile for one of the 24 designed catalysts.

![Kinetic grid availability; no numerical TOF result](../examples/plots/microkinetic_tof_heatmap.svg)

## 18. Campbell sensitivity and the distinction between batch and steady state

The energy-based degree of rate control measures the logarithmic rate response to an individual transition-state free energy with all other transition-state and intermediate energies fixed:

$$X_j=-RT\left(\frac{\partial\ln r}{\partial G_{TS,j}}\right)_{G_{k\ne j},G_{intermediates}}.$$

This definition is valuable because it concerns sensitivity of the complete mechanism, not simply the largest barrier measured from a conveniently chosen zero. The primary degree-of-rate-control framework also explains why several steps can share control and why intermediate stabilization can inhibit turnover. [Campbell, ACS Catalysis 2017, DOI 10.1021/acscatal.7b00115](https://www.osti.gov/pages/servlets/purl/1534882).

In the implementation, a central finite difference perturbs one transition state by $\pm\delta G$. Both its forward and reverse rate constants are multiplied by the same appropriate exponential factor, preserving that step's equilibrium constant. Perturbing only the forward constant would alter thermodynamics and would not evaluate the stated derivative. Intermediate free energies and initial concentrations remain fixed throughout the comparison.

A closed batch system changes composition with time and ultimately approaches zero net rate at equilibrium. Its finite-time average-TOF sensitivity is therefore explicitly called `finite_time_apparent_DRC_of_average_TOF`. It is not labeled a verified steady-state Campbell DRC. A nonpositive net TOF makes the logarithm undefined and produces an unavailable sensitivity rather than an arbitrary large number. No sum-to-one identity is imposed on the finite-time quantity.

This distinction matters for the requested conclusion about condensation versus dehydrogenation. A slow early condensation event may build an intermediate reservoir; later control can shift as imine, water, or hydrogenated catalyst accumulates. A rigorous claim about a rate-controlling step requires a complete valid network and a specified observable. At the present frozen evidence level, the relative roles of condensation and dehydrogenation are identifiable research questions, not established outcomes.

The finite-difference derivative introduces its own numerical question. A useful sensitivity result should be stable when the perturbation is reduced while solver tolerances remain tight enough to resolve the rate change. An excessively large energy shift measures a nonlinear response over an interval, whereas an excessively small shift can amplify integration noise. The analytical release test evaluates the derivative against a known expression, but a real multistep network still needs its own sensitivity-convergence assessment before small differences between coefficients are interpreted.

## 19. Surrogate learning, grouped validation, and ordinal Pareto analysis

The surrogate target is an accepted computed dehydrogenation activation free energy at one specified temperature and electronic/solvent protocol. Its inputs may include geometrical angles, buried volume, a documented formal d-electron count, and specified topological descriptors. The target is not inferred from a ground-state energy or an optimization success flag. Each label requires an accepted first-order saddle, reaction connectivity, complete reactant references, source identities, and consistent standard states. A training set also uses one activation-reference convention and reactant molecularity. A barrier relative to one preassociated complex cannot be mixed with a barrier relative to separated catalyst and alcohol, even when their total elemental compositions match; the missing association free energy would change the prediction target. This convention and the complete stoichiometric references remain in the saved training manifest.

The training gate requires at least ten accepted labels from distinct catalysts and at least five independent ligand groups. The campaign's initially requested six representative transition-state searches would not, even if all succeeded, satisfy this minimum. That is a data limitation rather than a reason to duplicate conformers or invent labels. `train_surrogate` returns `insufficient_computed_data` before training if those requirements are unmet.

Groups are defined by backbone and phosphine substituent and kept together across metals. Five-fold grouped cross-validation therefore tests prediction for held-out ligand combinations rather than allowing closely related cross-metal structures to populate both training and test sets. Duplicate source labels are rejected, and preprocessing or model fitting occurs only on the training portion of each fold. This addresses one important leakage route; it does not guarantee extrapolation to entirely new ligand chemistry. [Scikit-learn grouped cross-validation documentation](https://scikit-learn.org/1.1/modules/cross_validation.html#cross-validation-iterators-for-grouped-data).

A group-bootstrap gradient-boosting ensemble produces prediction means and a spread across members. That spread is reported as ensemble disagreement, not as a calibrated confidence interval or a direct estimate of experimental barrier error. Out-of-fold MAE and RMSE characterize the supplied calculated labels. An out-of-training-range flag identifies descriptor extrapolation, while the absence of that flag cannot prove chemical similarity.

Pareto analysis requires valid kinetic predictions at matched temperature, loading, base conditions, and duration. The comparator also checks the recorded initial concentration vector, actual free-base and conjugate-acid concentrations, and computational/solvent protocol. The declared temperature and loading must agree with the calculation's temperature and conserved initial catalyst/benzyl inventory. It maximizes average TOF while minimizing an explicitly ordinal metal-preference score: Mn and Fe are assigned the most preferred category, Co the next, and Ru the least preferred category in accordance with the user's ranking. These scores are not market prices. Without valid TOFs, there is no defensible Pareto-optimal catalyst to identify; geometry alone cannot populate that frontier.

## 20. Interpretation, uncertainty, and the next scientific decision

The completed bounded calculations demonstrate that the workflow can construct and electronically relax multiple designed pincer analogues and extract reproducible geometrical descriptors. They also reveal concrete limitations in starting-structure generation and differences in optimized donor arrangements that deserve chemical inspection. These are useful computational observations even when a complete catalytic ranking is unavailable. Their value depends on preserving the distinction between what was calculated and what the final research objective would require.

The strongest immediate comparison is within a fixed composition, charge, method, solvent convention, and structural identity. Best-found conformers can be compared and revisited, donor distances can be checked, and significant geometry rearrangements can be identified. Cross-composition absolute-energy ranking is not meaningful. Cross-metal catalytic ranking further requires valid transition states, compatible thermochemistry, and an adequate treatment of spin and speciation. The present GFN2 model cannot settle all of those questions by increasing the number of optimizer steps.

The toluene/tBuOK condition required a separate set of ALPB energy and force calculations with consistent reference states, which were carried out in this campaign. The resulting isolated-ion and potassium-contact-pair models retain distinct chemical interpretations; their numerical characterization does not determine the solution populations. The original gas-phase calculations remain legitimate preparatory evidence, but they are not retroactively relabeled as solution calculations. The same principle applies to unsuccessful paths: an unconverged or wrong-mode candidate remains part of the audit trail and does not become an accepted barrier when the report is formatted.

Uncertainty has several distinct sources: numerical electronic convergence, conformer coverage, local curvature, approximate Hamiltonian, solvent/speciation treatment, network completeness, integration error, and surrogate generalization. A small error in one category does not bound another. The Wilson interval on buried volume, for example, says nothing about an omitted potassium ion, and a low cross-validation error against computed labels says nothing about a missing competing mechanism.

A complete scientific conclusion would require accepted solution-model state free energies, certified pathways for every included elementary step, a physically interpretable base model, and enough independent labels for prediction. The report therefore preserves unavailable barriers, TOFs, rate-control coefficients, and surrogate rankings as unavailable at the frozen evidence level. It provides a substantive computational record and the formal conditions for advancing the study, while withholding the mechanistic and performance claims that the current evidence cannot support.
