# Pincer Catalysis Across Spin, Solvation, Kinetics and Geometry

## 1. Scope, evidence classes and the baseline chemical question

**Phase 3 scientific and mathematical whitepaper. English edition.** This report accompanies the executable `pincer-catmech-ai` repository. It develops a constrained crossing optimizer, an explicit alcohol cluster generator, an eighteen-state kinetic model and a native PyTorch equivariant neural network. Its central purpose is to make the mathematical objects, the calculations actually performed, and the remaining chemical uncertainties independently inspectable.

The campaign concerns alcohol–amine coupling through alcohol dehydrogenation, aldehyde–amine condensation and imine hydrogenation, with competing catalyst sequestration or destruction. The baseline thermochemical temperature is 383.15 K and the semiempirical solvent protocol is ALPB(toluene). The requested potassium tert-butoxide loading is 0.05 equivalents, with a 0.01–0.20 equivalent scan. Total added base, free tert-butoxide concentration and neutral tert-butanol activity are different quantities. No measured equilibrium model connecting them has been supplied.

Four evidence classes govern every subsequent interpretation. **Computed model results** come from retained native electronic-structure calculations with their stated approximations. **Analytical results** follow from explicitly stated mathematical models. **Software fixtures** use deliberately assigned coefficients to verify algorithms. **Chemical hypotheses** describe mechanisms that require additional calculation or experiment. These labels are not interchangeable. An electronic minimum is not a transition state; a solver trajectory is not a measured kinetic profile; symmetry compliance is not predictive accuracy.

The existing campaign ledger contains zero accepted transition states. Phase 3 therefore has no complete, accepted activation-free-energy network from which physical catalyst turnover frequencies or catalyst rankings can be calculated. It does contain nine new optimized microsolvation clusters, three selected Hessian-characterized minima, twenty-seven temperature-dependent association free energies, and a trained auxiliary conformer-energy benchmark. The latter does not beat its simple equal-reference baseline. These negative and incomplete outcomes are part of the result.

The four thematic parts occupy Sections 2–9, 10–18, 19–27 and 28–34. Sections 35–36 describe reproducibility and the conditions required for stronger conclusions. Equations use named species rather than ambiguous numerical indices wherever practical; machine-readable tables retain the full numerical precision.

## 2. Part I — Why compare Mn, Fe and Co with Ru?

The intended comparison concerns differences in accessible electronic configurations, ligand participation and reaction pathways across 3d and 4d metal complexes. The metal symbol alone cannot determine the spin state. Charge, donor field, coordination geometry, ligand protonation, covalent bonding and electronic method all enter the energy difference between competing states. A low-spin Ru reference can be a useful starting hypothesis, but this report does not treat every Ru complex as necessarily closed shell or every Mn, Fe or Co complex as necessarily undergoing spin crossover.

Let an electronic state be described by its nuclear geometry, electron count and wavefunction or density. A formal oxidation-state assignment is an electron-bookkeeping convention; it does not directly reveal the distribution of spin density. A ligand that exchanges protons during hydrogen transfer is mechanistically participating. A ligand that accepts substantial redox density is electronically non-innocent. These are related possibilities, but neither follows merely from the presence of a pincer ligand.

For a nominal spin quantum number S, multiplicity is

$$
M=2S+1.
$$

The pure-spin expectation is

$$
\langle S^2\rangle=S(S+1)=\frac{M^2-1}{4}.
$$

Thus singlet, triplet and quintet determinants target different electronic sectors when the electron-count parity allows them. Comparing their energies requires the same nuclear geometry and Hamiltonian. Comparing separately optimized structures answers a different question: the relative energies of different minima. Following a reaction across two surfaces requires still more information about their intersection and coupling.

The implemented crossing workflow is motivated by the constrained-surface problem discussed by [Harvey, Aschi, Schwarz and Koch](https://doi.org/10.1007/s002140050309). That methodological precedent does not establish spin crossover for the current catalyst matrix. Current calculations are reported as diagnostic calculations at their actual level of theory, with no transfer of small-model results into pincer labels.

## 3. Matched electronic surfaces and unit discipline

For N atoms collect the Cartesian coordinates into a vector R with 3N components. The crossing optimizer stores R in bohr, electronic energies in Hartree, gradients in Hartree per bohr and Hessians in Hartree per bohr squared. Molecular files and the graph model use angstrom instead. Conversion occurs at explicit interfaces: differentiating with respect to a coordinate in a different unit without converting the derivative changes the optimization problem numerically.

Two matched state evaluations provide

$$
E_1(\mathbf R),\quad E_2(\mathbf R),\quad \mathbf g_1=\nabla E_1,\quad \mathbf g_2=\nabla E_2.
$$

Define the average energy A, gap c, average gradient g and gap gradient d:

$$
A=\frac{E_1+E_2}{2},\qquad c=E_1-E_2.
$$

$$
\mathbf g=\frac{\mathbf g_1+\mathbf g_2}{2},\qquad \mathbf d=\mathbf g_1-\mathbf g_2=\nabla c.
$$

Matched means identical atom order, charge, basis, solvent model, functional and nuclear coordinates. Only the electronic-state specification changes. An apparent crossing obtained by subtracting energies from different compositions or solvents has no defined meaning in this optimization. A failed SCF evaluation is also not a very high energy state: it is missing evidence and causes the evaluation to fail.

The regular crossing seam is the set of geometries where c vanishes while d is nonzero. Locally it removes one degree of freedom. At such geometries the two state energies and their average coincide. The optimization problem is therefore

$$
\min_{\mathbf R} A(\mathbf R)\quad \mathrm{subject\ to}\quad c(\mathbf R)=0.
$$

The symmetry of A under state exchange is convenient. Minimizing E1 on the same seam produces the same stationary geometries, because A differs from E1 by a multiple of the constraint. This changes a Lagrange multiplier but not the constrained physical question. A near-zero d is a singular case and cannot be repaired by dividing it by a small arbitrary number and declaring success.

## 4. Deriving the MECP first-order conditions

Introduce a scalar multiplier lambda and define the Lagrangian

$$
L(\mathbf R,\lambda)=A(\mathbf R)+\lambda c(\mathbf R).
$$

Stationarity with respect to the multiplier enforces the crossing. Stationarity with respect to the geometry requires that the average force has no component along any allowed tangent displacement:

$$
\partial_\lambda L=c=0.
$$

$$
\nabla_{\mathbf R}L=\mathbf g+\lambda\mathbf d=\mathbf 0.
$$

These are the first-order Karush–Kuhn–Tucker conditions for a regular equality constraint. The normal component of the average gradient need not vanish: it is balanced by the multiplier. Requiring both unconstrained state gradients to vanish would demand a substantially different and generally unnecessary condition.

At a regular trial geometry, choose the multiplier that minimizes the squared stationarity residual. Expanding the norm gives

$$
q(\lambda)=\mathbf g^T\mathbf g+2\lambda\mathbf d^T\mathbf g+\lambda^2\mathbf d^T\mathbf d.
$$

Taking the derivative and setting it to zero gives

$$
q'(\lambda)=2\mathbf d^T\mathbf g+2\lambda\mathbf d^T\mathbf d=0.
$$

$$
\lambda_*=-\frac{\mathbf d^T\mathbf g}{\mathbf d^T\mathbf d}.
$$

The denominator is positive whenever the seam is regular. The minimizer is unique because the second derivative is twice that positive denominator. Substituting lambda-star back into the residual will produce the tangent projection derived next.

This derivation distinguishes three objects that are sometimes conflated: an energy-gap residual, a first-order constrained stationary point and a local constrained minimum. The first two are necessary for an MECP candidate. The third requires curvature information. Even a locally certified MECP does not determine the transition probability between surfaces; the spin–orbit coupling and nuclear dynamics remain separate quantities.

A historical constrained crossing-seam method is described by [Farazdel and Dupuis](https://doi.org/10.1002/jcc.540120219). The present implementation exposes its KKT/SQP equations explicitly rather than attaching an unverified historical name to a particular block solver.

## 5. Tangent projection and the corrected convergence criterion

Normalize the gap gradient and define a tangent projector:

$$
\mathbf n=\frac{\mathbf d}{\|\mathbf d\|},\qquad P=I-\mathbf n\mathbf n^T.
$$

The scalar identity n-transpose-n equals one immediately gives

$$
P^T=P,\qquad P^2=P,\qquad P\mathbf d=0.
$$

Substituting the optimal multiplier yields

$$
\mathbf g+\lambda_*\mathbf d=\mathbf g-\mathbf d\frac{\mathbf d^T\mathbf g}{\mathbf d^T\mathbf d}=P\mathbf g.
$$

Because P annihilates the gradient difference, the two projected state gradients are identical. Thus the code checks the actual absolute energy gap and the tangent-gradient maximum independently:

$$
|c|\le10^{-4}\ \mathrm{Hartree}.
$$

$$
\|P\mathbf g\|_\infty\le3\times10^{-4}\ \mathrm{Hartree/bohr}.
$$

The gradient tolerance is the recorded default, not a universal precision guarantee. The absolute value on the gap is essential: a large negative gap must fail. Both conditions must pass at the same accepted geometry.

The directive's expression involving an absolute gradient difference dotted with the position vector is not the MECP stationarity condition. A vector norm cannot be dotted with a vector; componentwise absolute values introduce origin dependence. Even a signed dot product need not vanish. Consider the analytic bond-length surfaces E1=(r−1)²/2 and E2=(r−3)²/2. Their crossing lies at r=2. The gap gradient in two-atom Cartesian coordinates is (−2u,2u), giving d·R=2r=4 although degeneracy and tangent stationarity hold exactly. This is a mathematical counterexample, not a computed molecule.

Under translation, energies and internal-coordinate gradients do not change. Under rotation, g and d rotate together, and the projection transforms consistently. Testing the full optimization trajectory under these transformations is more informative than checking a coordinate-dependent scalar. The implementation includes such covariance tests, including rejection of negative-gap and equal-gradient false convergence cases.

## 6. Scalar penalties, projected restoring vectors and globalization

A complete rational scalar corresponding to the specified penalty term is

$$
F_{\mathrm{sat}}=A+\alpha\frac{c^2}{c^2+\epsilon}.
$$

Here alpha has energy units and epsilon has squared-energy units. Applying the quotient rule explicitly,

$$
\frac{d}{dc}\frac{c^2}{c^2+\epsilon}=\frac{2c(c^2+\epsilon)-2c^3}{(c^2+\epsilon)^2}=\frac{2c\epsilon}{(c^2+\epsilon)^2}.
$$

Therefore

$$
\nabla F_{\mathrm{sat}}=\mathbf g+\frac{2\alpha\epsilon c}{(c^2+\epsilon)^2}\mathbf d.
$$

At large absolute gaps the restoring coefficient decays as the inverse cube of the gap. A saturated energy penalty can consequently have a small gradient while the crossing constraint remains badly violated. It is implemented and finite-difference tested, but cannot alone certify a crossing.

A Harvey-type effective vector separates tangent descent from normal restoration:

$$
\mathbf G_H=P\mathbf g+\kappa c\mathbf d.
$$

The two components are orthogonal, so

$$
\|\mathbf G_H\|^2=\|P\mathbf g\|^2+\kappa^2c^2\|\mathbf d\|^2.
$$

For a regular seam and positive kappa, vanishing of this vector implies both first-order conditions. The vector is not automatically the derivative of a global scalar function. Its formula must not be treated as permission to discard Hessian terms arising from a position-dependent unit gradient in another proposed scalar.

The implemented line search instead uses the explicitly defined merit function

$$
\Phi=A+\rho|c|+\frac{\alpha_q}{2}c^2.
$$

For nonzero gap its directional derivative along p is

$$
D\Phi[\mathbf p]=\mathbf g^T\mathbf p+\rho\operatorname{sgn}(c)\mathbf d^T\mathbf p+\alpha_q c\mathbf d^T\mathbf p.
$$

At zero gap, the absolute-value term becomes rho times the absolute directional constraint derivative. A full linearized constraint step has d·p=−c, making the last two terms negative. Multiplier-aware scaling, a capped increasing quadratic penalty and Armijo backtracking globalize the local step. Their atomic-unit numeric cap is 1000. No penalty setting relaxes the separate convergence tests.

## 7. SQP elimination, damped BFGS and true seam curvature

At each iterate approximate the Lagrangian Hessian by a positive-definite matrix B. The local quadratic subproblem is

$$
\min_{\mathbf p}\ \mathbf g^T\mathbf p+\frac12\mathbf p^TB\mathbf p,
\qquad c+\mathbf d^T\mathbf p=0.
$$

Its two block equations are Bp+d lambda-new=−g and d-transpose-p=−c. Solve Bu=g and Bv=d, then substitute p=−u−lambda-new v into the second equation:

$$
\lambda_{\mathrm{new}}=\frac{c-\mathbf d^T\mathbf u}{\mathbf d^T\mathbf v},\qquad \mathbf p=-\mathbf u-\lambda_{\mathrm{new}}\mathbf v.
$$

Since B is positive definite and d is nonzero, d-transpose-v is positive. The code solves linear systems rather than forming a matrix inverse. A trust-radius cap limits Cartesian displacement before backtracking; after rescaling, the linearized gap correction is fractional.

For accepted displacement s, form the secant vector using the same multiplier at both endpoints:

$$
\mathbf y=(\mathbf g_{k+1}+\lambda_{k+1}\mathbf d_{k+1})-(\mathbf g_k+\lambda_{k+1}\mathbf d_k).
$$

Let b=s-transpose-Bs and q=s-transpose-y. Set theta to one when q is at least 0.2b; otherwise use theta=0.8b/(b−q). Define

$$
\widetilde{\mathbf y}=\theta\mathbf y+(1-\theta)B\mathbf s.
$$

This makes s-transpose-y-tilde at least 0.2b and supports the damped update

$$
B_{k+1}=B_k-\frac{B_k\mathbf s\mathbf s^TB_k}{b}+\frac{\widetilde{\mathbf y}\widetilde{\mathbf y}^T}{\mathbf s^T\widetilde{\mathbf y}}.
$$

A positive search matrix does not prove positive physical curvature. With Q spanning directions perpendicular to the seam normal and rigid motions, the relevant physical matrix is

$$
H_{\mathrm{seam}}=Q^T(\nabla^2A+\lambda_*\nabla^2c)Q.
$$

The optional curvature routine finite-differences matched analytic gradients along tangent directions, reports asymmetry and eigenvalues, and requires first-order convergence before labeling a local minimum. Negative tangent curvature identifies a constrained saddle. Positive nontrivial eigenvalues support local minimality only, with no guarantee that another region of the seam is lower.

To construct Q, the code normalizes each nonzero crossing-normal, translation and rotation column before its singular-value decomposition. Without this normalization, a large molecular radius can give rotation columns a much greater scale than a small but regular gap gradient; the numerical rank decision can then lose the crossing constraint. Normalization preserves the excluded span while avoiding that scale imbalance. The finite differences keep lambda-star fixed, so they estimate the full Lagrangian curvature, including the constraint contribution. An analytic check uses A=−2x−y²/2+z²/2 and c=x+y²/2 at the origin: lambda-star is two, and both constrained eigenvalues are positive one even though the Cartesian second derivative of A along y is negative one.

Let H-red denote the raw reduced finite-difference matrix. Before interpreting its symmetrized eigenvalues, the numerical minimum check also requires

$$
\max_{i,j}|H_{\mathrm{red},ij}-H_{\mathrm{red},ji}|\leq10^{-5}\ \mathrm{Hartree}/\mathrm{bohr}^{2}.
$$

This default absolute threshold is configurable as `antisymmetry_tolerance`; the result records `hessian_symmetry_pass` and the threshold in Hartree/bohr-squared. Positive nontrivial eigenvalues must exceed the default curvature threshold of 1e−5 Hartree/bohr-squared. `positive_curvature` describes the symmetrized eigenvalues, while `minimum_verified` additionally requires the symmetry check and both first-order conditions. Symmetrization alone must not conceal a nonconservative or noisy gradient field: an adversarial test with symmetric eigenvalues [1,1] and maximum antisymmetry 40 correctly fails the minimum check.

This is a local numerical check at one finite-difference displacement, not a displacement-convergence certificate or a Hessian noise bound. A symmetric gradient bias can still pass the symmetry check. Step-size sensitivity, electronic-state continuity and uncertainty in small curvatures remain separate physical validations; none of these analytic examples supplies catalyst energies or a nonadiabatic rate.

## 8. Physical backend capability and retained diagnostic outcomes

Changing nominal occupation in native GFN2-xTB is not sufficient to produce a validated competing-spin Hamiltonian. The official implementation distinguishes its spin-independent energy expression from the spin-polarized extension invoked through `--spinpol --tblite`. This distinction determines which outputs may become spin labels. [Official xTB spin-polarization documentation](https://xtb-docs.readthedocs.io/en/latest/spgfn.html).

The same 1.21 angstrom O2 geometry produced native GFN2 energies of −7.906752352880 and −7.904118275554 Hartree for requested multiplicities one and three. The largest gradient difference was approximately 6×10−15 Hartree/bohr. Both spin-polarized tblite invocations terminated with SIGSEGV and exit code 3. These four calls and their original output are retained in `data/phase3/spin/xtb_capability_001/`. They are capability diagnostics and do not describe the pincer catalysts.

The preserved Psi4 environment supplies analytic state-specific RKS/UKS gradients, or RHF/UHF for an HF request. Each state receives a fresh guess with fixed orientation and center of mass. The implementation records geometry identity, energy, gradient, SCF convergence and the unrestricted spin diagnostic. Its orbital-overlap expression is

$$
\langle S^2\rangle=S_z(S_z+1)+N_\beta-\sum_{ij}|(C_\alpha^TS_{AO}C_\beta)_{ij}|^2.
$$

A deviation from the nominal pure-spin value above 0.1 is flagged for review. This is a screening threshold rather than a proof of state quality. Unrestricted DFT determinants, broken symmetry and near-degeneracy require additional interpretation. The provider follows the documented [Psi4 SCF](https://psi4.github.io/psi4docs/master/scf.html) and [analytic-gradient interfaces](https://psi4.github.io/psi4docs/master/opt.html).

The FeH2 diagnostic illustrates the distinction between execution and convergence. An initial singlet PBE/def2-SVP SCF failed within 120 iterations. A triplet/quintet optimization made 34 matched-pair calls in 300.204 seconds but exhausted its budget. Its last accepted gap was −1.578660260293×10−4 Hartree and tangent-gradient maximum 0.0140413212782 Hartree/bohr. No MECP was certified. A small intermediate gap had insufficient tangent stationarity. FeH2 is an artificial starting model, with no authority to supply pincer energies.

### Completed target matrix and crossing-search gate

The completed Fe/Co/Mn scan contains 54 multiplicity slots across 18 designs: 48 states had actual target-geometry attempts and 6 slots lacked the intended active geometry. SCF converged for 8 states; 6 passed the declared S² and identity diagnostic. There is 1 raw matched vertical gap and 0 diagnostically eligible pairs. Missing or failed values remain blank in the full-precision CSVs; quality-rejected raw values are retained with flags.

All states in this final matrix use gas-phase PBE/STO-3G, a 35-radial by 110-spherical integration grid, density fitting, SAD guesses, and late SOSCF activation at 1e−4. Energy and density convergence targets are 1e−8 and 1e−6 respectively. Each state has at most 120 SCF iterations, 120 seconds, 500 MiB configured memory, a 600 MiB observed worker-RSS ceiling and two threads. These measured resource bounds do not establish method accuracy. Co bipyridine active geometries rejected in the preceding campaign are not replaced by rearranged structures. Earlier resource/protocol pilots and an interrupted quintet attempt remain separate from the final 54-slot matrix.

Q means SCF converged and the stated S²/identity screens passed; U means SCF converged but quality was rejected; T is a time limit; M is missing source geometry; F is another retained failure. Q is a vertical diagnostic, not a validated chemical spin label. Multiplicity M is unrelated to molar concentration units here.

| Catalyst design | M=1 | M=3 | M=5 |
|---|---|---|---|
| Co_bipyridine_pnnoh_Ph | M | M | M |
| Co_bipyridine_pnnoh_iPr | M | M | M |
| Co_macho_pnp_Ph | T | T | T |
| Co_macho_pnp_iPr | T | T | T |
| Co_pyridine_pnn_Ph | T | T | T |
| Co_pyridine_pnn_iPr | T | T | T |
| Fe_bipyridine_pnnoh_Ph | Q | T | T |
| Fe_bipyridine_pnnoh_iPr | Q | U | T |
| Fe_macho_pnp_Ph | T | T | T |
| Fe_macho_pnp_iPr | T | T | U |
| Fe_pyridine_pnn_Ph | T | T | T |
| Fe_pyridine_pnn_iPr | T | T | T |
| Mn_bipyridine_pnnoh_Ph | Q | T | T |
| Mn_bipyridine_pnnoh_iPr | Q | T | T |
| Mn_macho_pnp_Ph | T | T | T |
| Mn_macho_pnp_iPr | T | T | T |
| Mn_pyridine_pnn_Ph | Q | T | T |
| Mn_pyridine_pnn_iPr | Q | T | T |

| Catalyst | M | Energy / Hartree | S² | Quality flag |
|---|---:|---:|---:|---|
| Fe_bipyridine_pnnoh_iPr | 1 | -2533.7122436149 | -0.00000000 | Q |
| Fe_bipyridine_pnnoh_iPr | 3 | -2533.7437027437 | 2.66144576 | U |
| Fe_bipyridine_pnnoh_Ph | 1 | -2756.8422312402 | -0.00000000 | Q |
| Fe_macho_pnp_iPr | 5 | -2713.7056752371 | 6.16690202 | U |
| Mn_bipyridine_pnnoh_Ph | 1 | -2756.7258829205 | 0.00000000 | Q |
| Mn_bipyridine_pnnoh_iPr | 1 | -2533.5901894093 | 0.00000000 | Q |
| Mn_pyridine_pnn_Ph | 1 | -2687.2566169593 | 0.00000000 | Q |
| Mn_pyridine_pnn_iPr | 1 | -2464.1132928431 | 0.00000000 | Q |

Source: `data/phase3/spin/pincer_vertical_diagnostic_003/spin_states.csv` and `vertical_gaps.csv`. The target MECP orchestrator then rechecked the completed matrix, native hashes, current structural identity, engine fingerprints and spin quality. Its actual receipt is `no_eligible_pair`, with 0 new pair evaluations and 0 new state attempts; first-order crossing = False, local minimum verified = False. A missing admissible pair cannot establish that a crossing is absent. The prospective search shares an eight-pair/1920-second budget with curvature checks; a full catalyst curvature calculation generally requires substantially more evaluations and cannot be certified by exhausting that budget.

![Final physical spin-state attempts, with raw rejected gaps distinguished](../examples/plots/phase3_spin_diagnostics.png)

The figure retains spin-contaminated raw gaps with rejection labels. Its energies are gas-phase small-basis diagnostics, not ALPB(toluene) spin free energies, activation barriers or intersystem-crossing rates.

## 9. What a target spin result would establish

Target pincer results are kept separate from O2 and FeH2 diagnostics. The matrix spans the available active Fe, Co and Mn geometries at nominal singlet, triplet and quintet specifications. Missing geometry slots must appear as missing records. A geometry that fails the prior coordination and identity gate cannot be silently reconstructed or promoted to an accepted molecular structure to complete a rectangular table.

The resource workflow first measures a PBE/def2-SVP target calculation. A larger-basis timeout motivates a separately labeled PBE/STO-3G gas-phase diagnostic protocol with a coarse integration grid. Coarse-grid or minimal-basis values are not ALPB(toluene) free energies. Distinct protocols remain in distinct native directories and summaries; a later success cannot overwrite an earlier failure. The final matrix counts and per-target data are reported by `data/phase3/spin/summary.json` and its linked `spin_states.csv` and `vertical_gaps.csv`. The completed target matrix and its explicit crossing-search gate are tabulated in Section 8.

For a matched vertical pair,

$$
\Delta E_{12}^{\mathrm{vertical}}(\mathbf R_0)=E_1(\mathbf R_0)-E_2(\mathbf R_0).
$$

This compares states at one geometry. It is neither a reaction barrier nor the height of a minimum-energy crossing. To calculate the latter, optimize the constrained geometry, certify the local seam minimum, identify its relationship to reactant and product valleys and use a consistent reference energy. To infer a physical intersystem-crossing rate, additionally characterize coupling and dynamics. An energy gap alone does not supply a probability per second.

No named “Chaban–Gordon–Dyall” nuclear MECP algorithm was verified in the primary sources checked for this campaign. The located [Chaban–Schmidt–Gordon article](https://doi.org/10.1007/s002140050241) concerns orbital optimization. The implementation therefore uses the descriptive KKT/SQP name, alongside explicitly attributed crossing-seam literature.

Before chemical ranking, the unresolved checks include basis and functional sensitivity, SCF stability, solvation, spin contamination, ligand identity and actual reaction connectivity. These are evidence requirements for the next scientific claim; they do not negate the completed derivation or the usefulness of preserved failed attempts for improving the search.

The acceptance interface requires a spin-state-sensitive electronic-structure provider. Its legacy Boolean name, `spin_dependent`, denotes the capability to resolve distinct electronic spin surfaces; it should not be read as requiring an explicitly spin-dependent operator in every nonrelativistic calculation. A spin-free electronic Hamiltonian can yield different-spin energies through the spatial symmetry and antisymmetry of its allowed electronic states. Capability must therefore be assessed from the method and state treatment, not from a nominal occupation setting alone.

`SurfacePair.validate` rejects explicitly failed supplied SCF, molecular-identity or diagnostic-quality checks, flagged spin contamination, and declared nonanalytic gradients. If S-squared metadata is supplied, it must be finite and have a valid reference. A reference supplied with multiplicity M must agree with (M²−1)/4, and eligibility requires

$$
\left|\langle S^2\rangle-\frac{M^2-1}{4}\right|\leq0.1.
$$

This is the project's conservative screening policy, not a universal spin-purity criterion or a proof of orbital stability, root continuity or quantitative accuracy. Failing supplied metadata cannot be bypassed by satisfying the mathematical crossing residuals. Missing quality metadata also proves nothing: metadata-free analytic toy surfaces remain available specifically for software verification.

Supplied state identities, multiplicities, energies and gradients must agree with the returned pair; the stated geometry and atomic units must agree with the actual evaluation. Supplied pair/state method, basis, charge, solvation and related provenance must be matched. Optimization and curvature evaluation reject changes in tracked state identities or the supplied electronic-structure protocol between geometries. These checks establish consistency of recorded metadata, not independent orbital-overlap tracking. The targeted spin suite passed 85 tests, including 57 added adversarial cases for quality rejection, protocol consistency, scaled projections and constrained curvature. Optimizer convergence still returns a first-order crossing with `minimum_verified=false`; the optional second-order result has the numerical limits described in section 7.

## 10. Part II — Explicit solvent as a local structural variable

A continuum solvent model supplies an effective response of the surrounding medium. It does not enumerate which individual tert-butanol molecule donates a proton, accepts a hydrogen bond or bridges two substrate sites. Adding explicit molecules therefore asks a different structural question while retaining the continuum outside the cluster. The present model combines GFN2-xTB with ALPB(toluene); it does not assert that continuum methods universally fail or that an explicit cluster automatically supplies a correct solution mechanism. Method sources are [GFN2-xTB](https://doi.org/10.1021/acs.jctc.8b01176) and [ALPB](https://doi.org/10.1021/acs.jctc.1c00471).

The substrate is the neutral hemiaminal with validated SMILES **`OC(Nc1ccccc1)c1ccccc1`**, formula C13H13NO and 28 atoms. The explicit solvent is `CC(C)(C)O`, formula C4H10O and 15 atoms. The one-, two- and three-alcohol clusters therefore have 43, 58 and 73 atoms. All have total charge zero and zero nominal unpaired electrons.

The verified hemiaminal reference has an NH orientation directed away from its OH oxygen. A simple rigid one-alcohol bridge could not meet the requested directional contacts. The generator consequently builds declared open OH/NH-anchored arrangements. With three alcohols, one extends an existing solvent OH contact to form an open relay. It records closure probes and failed distances instead of labelling those arrangements as cyclic dehydration transition states.

There is no metal catalyst, potassium ion or free tert-butoxide anion in these calculations. This limitation matters chemically: the base can alter protonation and electrostatic stabilization, while a catalyst can alter the electronic structure of bound reactants. The current clusters isolate one tractable local association question. They show which sampled neutral motifs are minima under the specified approximate Hamiltonian and how their standard association thermochemistry changes with temperature.

The distinction between an initial docking topology and an optimized hydrogen-bond graph is preserved. A contact that rearranges during relaxation remains a physical outcome of that calculation rather than being restored with an unreported geometric restraint.

## 11. Deriving the directional docking construction

Let D be a mapped donor atom and H its covalently attached hydrogen. The outward donor direction is

$$
\mathbf u=\frac{\mathbf R_H-\mathbf R_D}{\|\mathbf R_H-\mathbf R_D\|}.
$$

For a trial random vector z, remove its component along u and normalize the remainder:

$$
\mathbf v=\frac{\mathbf z-(\mathbf z\cdot\mathbf u)\mathbf u}{\|\mathbf z-(\mathbf z\cdot\mathbf u)\mathbf u\|}.
$$

Degenerate near-parallel draws are rejected. With an angular deviation theta between zero and fifteen degrees, define an acceptor direction and position:

$$
\mathbf d=\cos\theta\,\mathbf u+\sin\theta\,\mathbf v.
$$

$$
\mathbf R_{O_A}=\mathbf R_H+\ell_{HB}\mathbf d,\qquad \ell_{HB}=1.85\ \mathrm{angstrom}.
$$

Since u and v are orthonormal, d has unit length. The H-to-acceptor separation is therefore exactly the prescribed docking distance. The donor–H–acceptor angle is close to linear, by construction. A proper rotation Q orients the incoming solvent fragment around its reference oxygen:

$$
\mathbf R'_a=Q(\mathbf R_a-\mathbf R_O)+\mathbf R_{O_A},\qquad Q^TQ=I.
$$

All intrafragment distances are preserved because an orthogonal matrix preserves vector norms. The generator can therefore explore orientation without distorting an alcohol's covalent bonds. Unselected interfragment pairs must satisfy a steric lower bound of 0.65 times the sum of their van der Waals radii. Intended hydrogen bonds are screened independently, and the complete mapped covalent graph excludes accidental covalent attachments.

A bounded placement search allows twelve whole-placement attempts; each of the nine actual seeds succeeded on its first attempt. Initial contacts are checked against the requested 1.6–2.1 angstrom interval. Post-optimization retention uses the separately declared 1.4–2.6 angstrom H-to-acceptor interval and a donor–H–acceptor angle of at least 120 degrees, while also discovering all final interfragment N/O–H···N/O contacts.

These numerical cutoffs define the construction and diagnostic screens. They do not constitute a universal hydrogen-bond definition and do not prove proton transfer. Their purpose is to make success, failure and rearrangement reproducible from the same atom mapping and geometry.

## 12. Balanced association energies and consistent references

Let S denote the hemiaminal and A the neutral alcohol. The balanced association reaction is

$$
S+nA\rightleftharpoons SA_n.
$$

The corresponding electronic-model association energy is

$$
\Delta E_{\mathrm{assoc}}(n)=E(SA_n;\mathbf R_n^*)-E(S;\mathbf R_S^*)-nE(A;\mathbf R_A^*).
$$

Stars mark the sampled relaxed geometries under the same method and solvent convention. An energy comparison that omits n alcohol reference energies would compare different atom counts. An energy comparison that uses a gas-phase alcohol reference against a solvated cluster would mix protocols. Neither is the quantity calculated here.

Fresh same-protocol reference single points give

$$
E(S)=-1127.825737362139\ \mathrm{eV}.
$$

$$
E(A)=-482.70190329033903\ \mathrm{eV}.
$$

The selected one-alcohol cluster has energy −1611.004389644162 eV. Subtracting E(S)+E(A) gives an association energy of approximately −0.47675 eV, reported from full-precision arithmetic as −10.9940929259 kcal/mol. The table conversion is 1 eV=23.060547830619 kcal/mol. Large absolute energies are retained so the balanced subtraction can be reconstructed, while comparisons use the smaller differences appropriate to the scientific question.

Native xTB TOTAL ENERGY under this ALPB protocol includes the excess solvation contribution. Accordingly, the expression here is an ALPB-corrected total model association energy. It is often described as an electronic association energy, but it is not a pure gas-phase interaction energy and it is not a Gibbs free energy. The solvent term must not be added again. [Official xTB solvation documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html).

The directive called this balanced difference a cooperative stabilization energy. Its exact interpretation is association: more alcohols introduce more interactions, so more negative values do not by themselves demonstrate many-body cooperation. A separate frozen-fragment decomposition is required to test additivity, as derived next. This terminology prevents the attractive total interaction from becoming an unsupported claim about a special proton-wire mechanism.

## 13. Deriving the many-body remainder and deformation cost

Extract m=n+1 fragments from one selected cluster geometry without relaxing them. For each frozen monomer compute Ei and for every unique pair compute Eij with the same electronic and continuum protocol. Define each pair increment as

$$
\varepsilon_{ij}^{(2)}=E_{ij}-E_i-E_j.
$$

The complete frozen interaction is the cluster energy minus all monomer energies:

$$
E_{\mathrm{int}}^{\mathrm{frozen}}=E_{1\ldots m}-\sum_iE_i.
$$

Subtracting all pair contributions defines the beyond-pair remainder. Expanding that subtraction and counting each monomer in m−1 pairs gives

$$
E_{\mathrm{beyond\ pair}}=E_{1\ldots m}-\sum_{i<j}E_{ij}+(m-2)\sum_iE_i.
$$

For two fragments the expression is identically zero. For three fragments it becomes the usual complete three-body increment:

$$
\varepsilon_{123}^{(3)}=E_{123}-E_{12}-E_{13}-E_{23}+E_1+E_2+E_3.
$$

For four fragments it includes all three-body increments and the four-body increment together. Triple-subset calculations were not used to separate them, so no uniquely isolated four-body value is claimed.

The monomers may be strained relative to their isolated minima. Define their deformation cost and add it to the frozen interaction:

$$
E_{\mathrm{def}}=\sum_i E_i(\mathbf R_i^{\mathrm{cluster}})-E(S;\mathbf R_S^*)-nE(A;\mathbf R_A^*).
$$

$$
\Delta E_{\mathrm{assoc}}=E_{\mathrm{def}}+\sum_{i<j}\varepsilon_{ij}^{(2)}+E_{\mathrm{beyond\ pair}}.
$$

This identity is checked before unit conversion. For n=1, 2 and 3, deformation costs are +0.51122, +0.81026 and +0.78791 kcal/mol. The pair sums are −11.50531, −23.52359 and −36.06279 kcal/mol. The beyond-pair remainders are 0, +1.54707 and +3.55836 kcal/mol. Positive remainders reduce the attraction predicted by summing pairs; they do not demonstrate positive cooperative stabilization.

Each ALPB subset has its own cavity and reaction-field response. The remainder therefore includes continuum nonadditivity along with electronic nonadditivity. It cannot be attributed solely to hydrogen-bond polarization. There were 4, 7 and 11 frozen-fragment/full-cluster evaluations for the three sizes, including the redundant two-fragment full-pair check.

## 14. The complete nine-cluster association dataset

All nine new geometry optimizations converged with the original atom-mapped covalent identity retained. Each final structure is connected through its observed hydrogen-bond graph. One initially selected relay in the three-alcohol seed 02 changed; the final graph is stored separately from the initial chosen contacts. No geometry in this table is a transition state.

| Alcohols | Seed | Total energy / eV | Association / kcal mol−1 | Hessian selection |
|---:|---:|---:|---:|---|
| 1 | 00 | −1610.988973487 | −10.638588 | No |
| 1 | 01 | −1610.914578195 | −8.922992 | No |
| 1 | 02 | −1611.004389644 | −10.994093 | Minimum |
| 2 | 00 | −2094.093995332 | −19.934723 | No |
| 2 | 01 | −2094.147399853 | −21.166260 | Minimum |
| 2 | 02 | −2094.132699191 | −20.827255 | No |
| 3 | 00 | −2577.306805782 | −31.716522 | Minimum after repair |
| 3 | 01 | −2577.282715014 | −31.160975 | No |
| 3 | 02 | −2577.275422945 | −30.992816 | No |

The table is the full sampled set, rounded only for print. Source: `data/phase3/solvation/association_energies.csv`. Forces after optimization range from approximately 0.0004249 to 0.0008675 eV/angstrom. Those forces establish numerical relaxation under the configured stopping conditions, while the selected full Hessians supply the additional stationary-point classification.

![Actual GFN2-xTB/ALPB cluster association energies and thermochemistry](../examples/plots/phase3_micro_solvation.png)

**Figure 1. Computed neutral-cluster results.** The plotted association energies, association free energies and nonadditivity values come from the retained xTB calculations. They describe microsolvated minima. None of the plotted vertical differences is a dehydration activation barrier.

Three seeds per size constitute a bounded orientation sample, not an exhaustive search of all conformers. Selecting the lowest sampled electronic energy is a reproducible rule; it does not establish the global free-energy minimum. The seed-to-seed range itself illustrates the relevance of local arrangement. A larger ensemble could change the preferred motif or the association free energy through statistical weighting. Those possibilities are not represented by silently increasing the precision of the current three selected values.

The three-alcohol energy in this electronic table belongs to the initially selected unrestrained minimum candidate. Its later Hessian repair uses a slightly lower geometry for thermochemistry; the distinction is documented in the next section.

## 15. Full Hessians, external motions and stationary-point classification

An optimized geometry can have a small gradient while still being a saddle. The selected cluster from each size therefore receives a full numerical Hessian from two-sided differences of analytic gradients. For Cartesian components a and b,

$$
H_{ab}\simeq\frac{g_a(\mathbf R+\delta\mathbf e_b)-g_a(\mathbf R-\delta\mathbf e_b)}{2\delta}.
$$

The calculation used native `--hess --acc 0.05`, a 0.005 bohr displacement, scale one and no frozen atoms. Full unprojected Hessians are retained. This follows the [official xTB Hessian workflow](https://xtb-docs.readthedocs.io/en/latest/hessian.html); approximate or biased Hessians are not accepted in place of the specified full matrix.

Symmetrize the numerical matrix, mass-weight it and project out external motions. With diagonal atomic mass matrix M,

$$
F=M^{-1/2}\frac{H+H^T}{2}M^{-1/2}.
$$

An SVD of mass-weighted translations and rotations identifies six external directions for these nonlinear clusters. Let U span their orthogonal complement. The internal normal-mode problem is

$$
U^TFU\mathbf v_k=\lambda_k\mathbf v_k.
$$

With consistent atomic-unit conversion, the signed frequency is proportional to the signed square root of lambda. Negative eigenvalues must remain visible rather than being removed from the thermal calculation as inconvenient data.

The 43- and 58-atom clusters passed their first full-Hessian checks, with lowest positive frequencies 15.08035 and 17.76954 cm−1. The initial 73-atom cluster contained one −2.56861 cm−1 mode and was classified as numerical uncertainty. A documented displacement of at most 0.10 angstrom along that mode followed by unconstrained extreme optimization produced a second full Hessian with zero imaginary modes and a lowest frequency +1.42092 cm−1. Both attempts remain archived.

The repair lowered the Hessian-level electronic energy by 0.00030413 eV, about 0.0070 kcal/mol. The repaired geometry supplies the three-alcohol Gibbs values; the original selected geometry supplies its electronic association and many-body tables. The very soft final positive mode warns that low-frequency treatment and conformer sampling matter, even though the declared minimum gate passes. It is a numerical characteristic to investigate, not a dehydration mode.

## 16. Step-by-step qRRHO and one-molar standard states

For a positive vibrational frequency, define x=hc times the wavenumber divided by kBT. The harmonic oscillator partition function includes zero-point energy:

$$
q_k^{HO}=\frac{e^{-x_k/2}}{1-e^{-x_k}}.
$$

Using S=R[ln q+T times the temperature derivative of ln q] gives the harmonic vibrational entropy

$$
S_k^{HO}=R\left[\frac{x_k}{e^{x_k}-1}-\ln(1-e^{-x_k})\right].
$$

Low frequencies make the harmonic entropy sensitive to soft intermolecular motions. The implemented reduced-inertia qRRHO interpolation, motivated by [Grimme's thermochemical treatment](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201200497), smoothly mixes oscillator and free-rotor entropy. Define

$$
\mu_k=\frac{h}{8\pi^2c\widetilde\nu_k},\qquad \mu_k^{eff}=\frac{\mu_kB_{av}}{\mu_k+B_{av}}.
$$

$$
w_k=\frac{1}{1+(\widetilde\nu_0/\widetilde\nu_k)^4},\qquad \widetilde\nu_0=100\ \mathrm{cm}^{-1}.
$$

The free-rotor factor and entropy are

$$
q_k^{FR}=\sqrt{\frac{8\pi^3\mu_k^{eff}k_BT}{h^2}},\qquad S_k^{FR}=R[\ln q_k^{FR}+1/2].
$$

$$
S_{vib}^{qRRHO}=\sum_k[w_kS_k^{HO}+(1-w_k)S_k^{FR}].
$$

The repository includes translation, rotation and electronic entropy, states symmetry number one, retains harmonic zero-point and thermal enthalpy, and replaces vibrational entropy by this interpolation. The reference-state conversion follows directly from the ideal-gas concentration p°/(RT). Moving to concentration c° changes chemical potential by RT ln[c°RT/p°]. Thus

$$
G^{1M}=E_{model}+H_{thermal}-TS^{qRRHO}+RT\ln(c^\circ RT/p^\circ).
$$

Here c° is one molar and p° is 101325 Pa, with coherent volume units inside the dimensionless logarithm. H-thermal includes zero-point energy. ALPB excess solvation is already present in E-model and is not added again.

For association S+nA to SAn, the change in molecule count is −n. Consequently its net reference-state correction is −nRT ln(c°RT/p°), not one positive correction per reaction. The reported association free energy subtracts all n+1 reactant standard free energies from the cluster standard free energy. These thermochemical manipulations use the computed spectrum but do not model actual bulk concentration or ion pairing.

## 17. Complete temperature-dependent association thermochemistry

The three selected cluster minima and two verified isolated references provide the twenty-seven values below. All entries are qRRHO association free energies in kcal/mol at the one-molar standard state, using GFN2-xTB/ALPB(toluene). The temperature extension reuses the same cached full Hessians. It changes partition-function evaluation, not the underlying electronic structure or solvent parametrization.

| Temperature / K | One alcohol | Two alcohols | Three alcohols |
|---:|---:|---:|---:|
| 298.15 | +0.658254 | +3.586112 | +4.105509 |
| 360.00 | +2.943523 | +8.415315 | +11.079699 |
| 370.00 | +3.312597 | +9.195240 | +12.205561 |
| 380.00 | +3.681538 | +9.974882 | +13.330908 |
| 383.15 | +3.797726 | +10.220408 | +13.685284 |
| 390.00 | +4.050341 | +10.754229 | +14.455731 |
| 400.00 | +4.419004 | +11.533271 | +15.580024 |
| 410.00 | +4.787521 | +12.312001 | +16.703782 |
| 420.00 | +5.155889 | +13.090409 | +17.826998 |

Source: `data/phase3/solvation/association_thermochemistry.csv`. The positive association free energies at the baseline contrast with the negative total-model association energies. This is a thermodynamic result within the stated approximation: favorable interactions do not automatically overcome the free-energy cost of assembling several independently translating species into one cluster.

The equilibrium relation uses dimensionless activities:

$$
K_n^\circ=\exp[-\Delta G_{\mathrm{assoc},n}^\circ/(RT)].
$$

$$
a_{SA_n}=K_n^\circ a_Sa_A^n.
$$

The exponent n makes cluster populations sensitive to alcohol activity. The table alone cannot determine those populations because actual activities, competing clusters and chemical speciation are unresolved. A standard one-molar association free energy should not be inserted directly into a statement that a particular cluster dominates under the requested 0.05-equivalent total-base condition.

A conformer ensemble would require a consistent partition sum over multiple distinct minima, including their degeneracies and free energies. Three electronic-energy seeds per size with one characterized minimum do not supply that ensemble. The retained table is valuable precisely because its scope is narrower and explicit: it is a reproducible comparison of selected local motifs, rather than a bulk equilibrium measurement.

## 18. From open minima to a testable proton-shuttle mechanism

The completed calculations establish attractive local association, positive beyond-pair remainders for the larger clusters, and positive standard association free energies over the reported temperatures. They establish neither a six-membered nor an eight-membered cyclic dehydration transition state. A set of hydrogen bonds drawn around an intermediate cannot be promoted to a transition state by naming it a proton wire.

For a proposed shuttle path, the required activation free energy is a balanced difference between a characterized saddle and its correctly assembled reactant reference:

$$
\Delta G_n^\ddagger=G^\circ(TS_n)-G^\circ(S)-nG^\circ(A).
$$

If the path proceeds from a preassociated cluster, the intracluster barrier instead subtracts G(SAn). The two reference choices are connected by association thermochemistry:

$$
\Delta G_n^\ddagger=\Delta G_{\mathrm{assoc},n}^\circ+\Delta G_{\mathrm{cluster},n}^\ddagger.
$$

This identity shows why electronic stabilization of a cluster cannot be used as a surrogate barrier lowering. The transition state may gain a different amount of stabilization, and assembling the cluster carries its own free-energy cost. Concentration-dependent rate expressions must also use the same molecularity and standard state.

A convincing dehydration saddle would require a small stationary gradient, one appropriate imaginary mode, its explicit displacement pattern and connectivity to the intended hemiaminal and imine-plus-water basins. The proton motion may be concerted or asynchronous; that outcome should emerge from the path rather than being imposed by the prose. Potassium, tert-butoxide, catalyst participation and water can change this landscape, so their eventual inclusion should be composition-balanced and separately recorded.

The bounded dehydration search is now complete. Two independent 43-atom hemiaminal/tBuOH seeds were evaluated with actual GFN2-xTB/ALPB(toluene), followed by one immutable continuation of the first seed's final nine-image band. The continuation preserved all mapped image coordinates and the source SHA256, recomputed every force and repeated both endpoint Hessians without resetting interpolation. The campaign used one xTB worker and one thread; physical calculations finished before the original 1200-second deadline. All three path runs remained unconverged, so the final counts are **0 accepted transition states and 0 activation-free-energy barriers**.

| Path run | Actual FIRE/NEB steps | Final NEB maximum force, eV/angstrom | Required maximum force, eV/angstrom | Final status |
| --- | ---: | ---: | ---: | --- |
| Original seed 02 | 100 | 0.375882 | 0.07 | Unconverged |
| Original seed 01 | 100 | 1.347282 | 0.07 | Unconverged |
| Immutable seed-02 continuation | 200 additional | 1.206630 | 0.07 | Unconverged |

These are 400 actual FIRE/NEB steps across two seeds and one continuation, not three independent starting conformers. Each original path certified its reactant and product endpoint through free optimization, the complete mapped covalent graph, the C-N distance interval and a full native Hessian with 123 positive internal frequencies. This gives four original endpoint certifications plus two repeated certifications during continuation; it does not establish six independently discovered minima. Correct C=N product geometries were fitted from accepted imine, water and tBuOH component references while preserving the identities of the transferred N-H and alcohol O-H protons. The seed-02 preparation's H-H close-contact flag remains in the record: free optimization removed the contact and certified the correct product basin. Preparation coordinates themselves are not stationary-point evidence.

The tested channel contains neutral tBuOH and no explicit potassium, tert-butoxide or metal catalyst. It therefore does not resolve the 0.05-equivalent total-tBuOK boundary, a free-base concentration or a base kinetic order. Endpoint thermochemistry includes the 383.15 K baseline, with 300 K electronic smearing kept distinct. Because the NEB force gate failed, saddle refinement, a saddle Hessian, imaginary-mode acceptance and opposite-endpoint descents were not reached. A maximum on any of these unconverged bands remains a sampled potential-energy diagnostic; no numerical activation barrier was assigned to the kinetic model or the EGNN training labels.

The final `data/phase3/proton_wire/summary.json` and `verification.json` retain every outcome. `sampled_band_profiles.csv` contains 27 actual image-profile rows, and `endpoint_thermochemistry.csv` contains 78 endpoint-temperature rows, including the repeated endpoint checks. Two complete native ZIPs, each below 95 MiB, preserve stdout, launch records, gradients, XYZ structures, full Hessians and electronic restart/cache files. The combined manifest verifies 51,322 archive-qualified member versions covering 51,321 current evidence-tree files; the extra version preserves an earlier helper source. All archive/member SHA256s and the exact nine-image restart coordinates passed verification. Checksums establish file identity, not a catalytic mechanism or a six-/eight-membered transition state.

Current data justify building and testing such hypotheses. They do not justify a numerical catalytic rate. This separation also governs machine learning: these minimum energies can support carefully matched conformer-energy tasks, whereas barrier training requires actual transition-state free-energy labels. The next part therefore develops the kinetic equations and verifies their software without using these association values to fill the missing activation barriers.

## 19. Part III — Eighteen states and the distinction between inventory and activity

The master kinetic model has eighteen dynamic concentrations. Each is measured in mol/L, time is in seconds and the standard concentration is one molar. The state vector is

$$
\mathbf c=(C,CA,CHB,CH,CHI,CP,D,CO,X,A,B,N,H,I,P,W,Z,U)^T.
$$

The symbols are bookkeeping names for defined model states. In particular, H means hemiaminal concentration; it is not a hydrogen atom or hydride concentration. The symbol CO denotes the poisoned catalyst concentration; it does not mean dissolved free carbon monoxide.

| Symbol | State meaning | Metals per entity |
|---|---|---:|
| C | Active catalyst | 1 |
| CA | Catalyst–alcohol | 1 |
| CHB | Hydrogenated catalyst–aldehyde | 1 |
| CH | Hydrogenated catalyst | 1 |
| CHI | Hydrogenated catalyst–imine | 1 |
| CP | Catalyst–desired product | 1 |
| D | Hydrogenated catalyst dimer | 2 |
| CO | CO-poisoned catalyst | 1 |
| X | Degraded catalyst | 1 |
| A | Benzyl alcohol | 0 |
| B | Benzaldehyde | 0 |
| N | Aniline | 0 |
| H | Hemiaminal | 0 |
| I | Imine | 0 |
| P | Desired secondary amine | 0 |
| W | Water | 0 |
| Z | Benzene | 0 |
| U | Free tert-butoxide | 0 |

Neutral tert-butanol enters the reaction compositions as a reservoir with prescribed dimensionless activity a. It appears on both sides of each proposed shuttle event and has zero net stoichiometric change. It is not a consumed species omitted from a material balance. Free base U is dynamic because cleavage consumes it, while base-assisted condensation retains it on both sides.

The model does not solve ion pairing, potassium coordination or the activation equilibrium that produces free base. Consequently, a numerical initial U corresponding to 0.05 equivalents is a **free-base fixture condition** and cannot be represented as an experimentally established consequence of adding 0.05 equivalents tBuOK. The scope of the kinetic state space is part of the model definition, not an experimentally verified speciation map.

## 20. All eighteen net rates and their molecularity

Let fj and bj denote the forward and reverse coefficients of path j, with one-based numbering in this report. The first six net rates describe the proposed borrowing-hydrogen cycle:

$$
r_1=f_1 C\cdot A-b_1(CA).
$$

To avoid any ambiguity from adjacent species symbols, the same first rate is written explicitly as r1=f1·C·A−b1·(CA), where (CA) is the bound state. The remaining algebraic expressions use parentheses around multi-letter bound states:

$$
r_2=f_2(CA)-b_2(CHB),\qquad r_3=f_3(CHB)-b_3(CH)\cdot B.
$$

$$
r_4=f_4(CH)\cdot I-b_4(CHI),\qquad r_5=f_5(CHI)-b_5(CP).
$$

$$
r_6=f_6(CP)-b_6 C\cdot P,\qquad r_7=f_7 B\cdot N-b_7H.
$$

$$
r_8=f_8H-b_8 I\cdot W,\qquad r_9=U(f_9 B\cdot N-b_9H).
$$

$$
r_{10}=U(f_{10}H-b_{10}I\cdot W),\qquad r_{11}=a(f_{11}H-b_{11}I\cdot W).
$$

$$
r_{12}=a^2(f_{12}H-b_{12}I\cdot W),\qquad r_{13}=a^3(f_{13}H-b_{13}I\cdot W).
$$

$$
r_{14}=f_{14}(CH)^2-b_{14}D,\qquad r_{15}=f_{15}C\cdot B.
$$

$$
r_{16}=f_{16}(CA)\cdot B,\qquad r_{17}=f_{17}C\cdot U,\qquad r_{18}=f_{18}(CA)\cdot U.
$$

The complete pathway mapping is listed below using multiplication dots for separate concentrations.

| Path | Net rate | Interpretation |
|---:|---|---|
| 1 | f1·C·A − b1·CA | Alcohol coordination |
| 2 | f2·CA − b2·CHB | Alcohol hydrogen transfer |
| 3 | f3·CHB − b3·CH·B | Aldehyde release |
| 4 | f4·CH·I − b4·CHI | Imine coordination |
| 5 | f5·CHI − b5·CP | Imine hydrogenation |
| 6 | f6·CP − b6·C·P | Product release |
| 7 | f7·B·N − b7·H | Hemiaminal addition |
| 8 | f8·H − b8·I·W | Direct dehydration |
| 9 | U·(f9·B·N − b9·H) | Base-assisted addition |
| 10 | U·(f10·H − b10·I·W) | Base-assisted dehydration |
| 11 | a·(f11·H − b11·I·W) | One-alcohol shuttle |
| 12 | a²·(f12·H − b12·I·W) | Two-alcohol shuttle |
| 13 | a³·(f13·H − b13·I·W) | Three-alcohol shuttle |
| 14 | f14·CH² − b14·D | Hydrogenated dimerization |
| 15 | f15·C·B | Decarbonylation |
| 16 | f16·CA·B | Bound-alcohol decarbonylation |
| 17 | f17·C·U | Phosphine-arm cleavage |
| 18 | f18·CA·U | Bound-alcohol arm cleavage |

Paths 15–18 use an explicitly declared irreversible-sink approximation, with reverse coefficients exactly zero. This approximation is not a deduction from a reaction energy. The dimerization event consumes two CH entities; its macroscopic f14 convention already defines the event rate, so no extra one-half factor is inserted.

For a reaction side with n dynamic molecules and m reservoir molecules, a full conventional coefficient has units M^(1−n−m)/s. Substituting reservoir concentration by c°a to the power m produces an effective coefficient with units M^(1−n)/s. Because c° is numerically one molar, this unit transformation may leave the stored coefficient's number unchanged. The units nevertheless matter in deriving rates or comparing protocols.

Powers a, a² and a³ are the molecularity hypotheses for the specified shuttle compositions. They are not measured kinetic orders. A catalyst-mediated proton relay could follow a different pre-equilibrium dependence, which would require a different state model and new evidence.

## 21. Expanding the complete eighteen-equation ODE

For each reaction, subtract the reactant stoichiometric coefficients from the product coefficients to construct a column of S. Rows follow the state order in Section 19. Multiplying this eighteen-by-eighteen matrix by the net rate vector yields the concentration derivative:

$$
\dot{\mathbf c}=S\mathbf r(\mathbf c).
$$

The implementation also writes the full derivative explicitly so the signs and stoichiometric multiplicities can be reviewed against the reaction list. The nine metal-containing states obey

$$
\dot C=-r_1+r_6-r_{15}-r_{17}.
$$

$$
\dot{CA}=r_1-r_2-r_{16}-r_{18}.
$$

$$
\dot{CHB}=r_2-r_3.
$$

$$
\dot{CH}=r_3-r_4-2r_{14}.
$$

$$
\dot{CHI}=r_4-r_5,\qquad \dot{CP}=r_5-r_6.
$$

$$
\dot D=r_{14},\qquad \dot{CO}=r_{15}+r_{16},\qquad \dot X=r_{17}+r_{18}.
$$

The nine remaining dynamic species obey

$$
\dot A=-r_1+r_{16}+r_{18}.
$$

$$
\dot B=r_3-r_7-r_9-r_{15}-r_{16}.
$$

$$
\dot N=-r_7-r_9.
$$

$$
\dot H=r_7+r_9-r_8-r_{10}-r_{11}-r_{12}-r_{13}.
$$

$$
\dot I=-r_4+r_8+r_{10}+r_{11}+r_{12}+r_{13}.
$$

$$
\dot P=r_6,\qquad \dot W=r_8+r_{10}+r_{11}+r_{12}+r_{13}.
$$

$$
\dot Z=r_{15}+r_{16},\qquad \dot U=-r_{17}-r_{18}.
$$

These are eighteen equations despite compactly sharing some display lines. The factor of two multiplying r14 is the key metal-inventory term. Bound-alcohol poisoning and cleavage release the previously coordinated A, explaining its positive r16 and r18 terms. Omitting those releases would violate elemental bookkeeping even if an unweighted catalyst count appeared plausible.

The code checks agreement between the explicit derivative and S times r. This cross-check can expose transcription errors because it compares two representations of the same model. It cannot establish that all these proposed elementary steps occur in solution. The model's chemical legitimacy still depends on structures, transition states and experiments that are outside an algebraic identity.

## 22. Exact metal, base, elemental and charge ledgers

Conservation is a left-nullspace property of the stoichiometric matrix. If a row vector l satisfies lS=0, then

$$
\frac{d}{dt}(\mathbf l\mathbf c)=\mathbf lS\mathbf r=0.
$$

The metal row assigns weight one to each monomeric catalyst state and weight two to D. The conserved inventory is therefore

$$
M_{tot}=C+CA+CHB+CH+CHI+CP+2D+CO+X.
$$

The sum of catalyst entity concentrations is not conserved during dimerization. Two monomers become one dimer, while two metal atoms remain two metal atoms. The required numerical tolerance applies to this correctly weighted inventory. The separate base-equivalent inventory is

$$
B_{tot}^{eq}=U+X.
$$

Base-assisted addition and dehydration cancel U from the stoichiometric column. Only paths 17–18 transfer one equivalent from U into the cleavage bookkeeping associated with X. Both ledger rows annihilate S exactly in integer arithmetic before floating-point integration begins.

The hydrogenated dimer contains twice the active-catalyst composition plus four H atoms under the stated hydrogenation hypothesis: two hydrides and two additional ligand protons. Cleavage removes a neutral fragment R–OtBu and makes the degraded catalyst's relative charge one lower than C. The neutral tBuOH reference is C4H10O; free tert-butoxide is C4H9O with charge −1. These composition and charge differences are enforced independently of species names.

The eliminated cleavage fragment need not add a nineteenth ODE because its concentration F follows exactly from X:

$$
F(t)=F(0)+X(t)-X(0).
$$

For element e, let a-e-i be its count in state i and a-e-F its count in the fragment. The effective dynamic ledger is

$$
L_{e,i}=a_{e,i}+a_{e,F}\delta_{i,X}.
$$

Then LeS=0. The actual elemental total adds the constant a-e-F times [F(0)−X(0)], which matters for a pre-existing degraded population but cancels from drift residuals. The same construction applies to charge. An inert potassium counterion contributes a constant. These exact ledgers verify the declared reaction model; they do not imply that real solution speciation has been exhaustively represented.

## 23. Deriving every nonzero analytical Jacobian entry

Define the rate derivative D with reaction rows and species columns:

$$
D_{jk}=\frac{\partial r_j}{\partial c_k}.
$$

Since S is constant, differentiation commutes with its finite sum:

$$
J_{ik}=\frac{\partial\dot c_i}{\partial c_k}=\sum_jS_{ij}D_{jk},\qquad J=SD.
$$

The production code differentiates each polynomial directly. For example, differentiating f1·C·A with respect to C holds A fixed and yields f1·A. Differentiating U·(f9·B·N−b9·H) with respect to U retains the entire parenthesis, while its derivative with respect to B retains U. These product-rule factors are crucial for the off-cycle sensitivities.

The complete nonzero rate-derivative entries are below. Each item `species:value` means the derivative with respect to that species. Every entry absent from this table is exactly zero.

| Rate | Nonzero derivatives |
|---:|---|
| 1 | C:f1·A; A:f1·C; CA:−b1 |
| 2 | CA:f2; CHB:−b2 |
| 3 | CHB:f3; CH:−b3·B; B:−b3·CH |
| 4 | CH:f4·I; I:f4·CH; CHI:−b4 |
| 5 | CHI:f5; CP:−b5 |
| 6 | CP:f6; C:−b6·P; P:−b6·C |
| 7 | B:f7·N; N:f7·B; H:−b7 |
| 8 | H:f8; I:−b8·W; W:−b8·I |
| 9 | B:U·f9·N; N:U·f9·B; H:−U·b9; U:f9·B·N−b9·H |
| 10 | H:U·f10; I:−U·b10·W; W:−U·b10·I; U:f10·H−b10·I·W |
| 11 | H:a·f11; I:−a·b11·W; W:−a·b11·I |
| 12 | H:a²·f12; I:−a²·b12·W; W:−a²·b12·I |
| 13 | H:a³·f13; I:−a³·b13·W; W:−a³·b13·I |
| 14 | CH:2·f14·CH; D:−b14 |
| 15 | C:f15·B; B:f15·C |
| 16 | CA:f16·B; B:f16·CA |
| 17 | C:f17·U; U:f17·C |
| 18 | CA:f18·U; U:f18·CA |

Together with the explicit S construction or all eighteen ODEs, this table determines every element of the eighteen-by-eighteen ODE Jacobian without displaying a mostly zero matrix. For example J-C,C=−f1·A−b6·P−f15·B−f17·U, and J-CH,CH includes −4f14·CH because two CH are consumed per dimer event.

No rate divided by concentration is used. Consequently the derivatives remain defined at zero concentration. Central finite differences are an independent test, not the implementation. The forty-condition benchmark reports a maximum absolute central-difference discrepancy of 9.630873876×10−11, below the requested 10−5 threshold. This numerical agreement validates the derivative of the declared rates, not the scientific provenance of any rate coefficient.

## 24. Radau integration, stiffness and independently checkable limits

The master system can become stiff when coordination, conversion and deactivation act on very different time scales. The implementation calls SciPy's Radau method with the analytic Jacobian, relative tolerance 10−9, absolute tolerance 10−12 M and dense output. The solver supports callable Jacobians and adaptive time stepping. [SciPy solve_ivp documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).

The requested five hundred times are a reporting grid. They are not a claim that exactly five hundred internal steps were accepted. Both the adaptive accepted solution values and the five hundred dense-output samples are audited. Every value must be finite; metal and base-equivalent drifts must remain below 10−10 M; concentrations below −10−10 M cause failure. The solver does not clip a negative trajectory back into compliance.

An analytical limiting case tests more than local Jacobian entries. Turn off every path except irreversible dimerization, with initial hydrogenated catalyst concentration h0. The CH equation becomes

$$
\frac{dh}{dt}=-2f_{14}h^2.
$$

Separate variables and integrate:

$$
\int_{h_0}^{h(t)}h^{-2}\,dh=-2f_{14}\int_0^t dt.
$$

$$
\frac{1}{h(t)}=\frac{1}{h_0}+2f_{14}t,\qquad h(t)=\frac{h_0}{1+2f_{14}h_0t}.
$$

The dimer concentration follows from h+2D conservation. This test detects both the stoichiometric factor and the quadratic concentration dependence. It complements randomized derivative comparisons and exact ledger checks.

Finite precision still matters. An analytic Jacobian does not make the numerical initial-value solution exact, and local error tolerances do not guarantee a particular global error for every stiff problem. The external mass and positivity audits therefore remain necessary. The observed benchmark's very small conservation residuals demonstrate numerical consistency for the sampled fixtures, without establishing physical rate accuracy.

The model is a finite batch with irreversible sinks. It is inappropriate to infer a steady state simply because a plot appears flat over a short window. A small product-release rate can instead reflect substrate depletion, catalyst sequestration or catalyst loss. Stability or bifurcation claims would require explicitly defined steady equations, a parameter continuation and nontrivial eigenvalue analysis after accounting for conserved inventories.

## 25. Eyring coefficients and the thermochemical provenance gate

For a reaction side with molecularity m, use standard chemical potentials at the same temperature and reference concentration. The transition-state activation free energy is

$$
\Delta G_j^{\ddagger,\circ}=G_{TS,j}^{\circ}-\sum_i\nu_{ij}^{react}G_i^\circ.
$$

The standard-concentration Eyring coefficient is

$$
k_j=\frac{k_BT}{h}(c^\circ)^{1-m}\exp[-\Delta G_j^{\ddagger,\circ}/(RT)].
$$

Multiplication by m concentration factors then gives M/s. Reservoir substitution follows the molecularity derivation in Section 20. Transmission coefficients or diffusion limits would require additional assumptions; none is inferred from an electronic association energy.

If both directions share one transition state and consistent endpoint energies, their ratio follows from subtraction of the two activation free energies. Parallel condensation routes then share their appropriate equilibrium ratios automatically. Arbitrarily perturbing only a forward coefficient would change equilibrium, even if the intended question were only a change to the transition-state height. This is why the control calculation scales both directions together.

Physical rate construction requires all twenty state/reference records and eighteen matched transition-state records. Required attributes include temperature-matched finite one-molar free energies, a common method/solvent/reference protocol, accepted stationarity and reviewed chemical identity, matching source-file hashes, full atom and charge closure, and appropriate transition-state mode/connectivity evidence. The project also records a −1800 to −300 cm−1 alcohol-dehydrogenation mode interval as an acceptance rule; it is not a universal imaginary-frequency law.

Missing pathways are not filled with reaction free energies or coefficients from a different catalyst. Negative or nonfinite standard-state barriers are rejected for review of association and diffusion treatment rather than clipped. Computed rate arrays are regenerated from their supplied thermochemical snapshot and compared against every stored coefficient and the stated temperature. Fixtures require explicit opt-in and cannot carry a computed snapshot.

Hash consistency proves the identity of the supplied bytes. It does not independently prove the scientific truth of a manually asserted energy or review flag. The validator checks supplied attestations and integrity; external chemical-identity and calculation review remains necessary. This boundary is explicit so a provenance schema is not mistaken for an automatic scientific oracle.

## 26. Exact parameter tangents and finite-time rate control

Define a dimensionless perturbation of transition state j at fixed intermediates, temperature and reservoir activity:

$$
\theta_j=-\delta G_{TS,j}/(RT).
$$

Both forward and reverse coefficients of that path acquire exp(theta-j), preserving their ratio. Let the sensitivity matrix Z have entries Z-i-j=partial c-i/partial theta-j. Since the initial state is fixed, Z(0)=0. Differentiating the governing ODE by the chain rule gives

$$
\dot Z=J(\mathbf c(t))Z+S\operatorname{diag}(\mathbf r(\mathbf c(t))).
$$

The second term follows because each net rate changes proportionally to itself when both directions scale. There are 324 tangent variables. Column-major vectorization produces

$$
\operatorname{vec}_F(\dot Z)=(I_{18}\otimes J)\operatorname{vec}_F(Z)+\operatorname{vec}_F(S\operatorname{diag}(\mathbf r)).
$$

The tangent solver therefore has an exact sparse block-diagonal Jacobian. It integrates along the accepted concentration trajectory, whose model signature prevents accidental mixing with another rate array or reservoir activity. Total rate derivatives are

$$
Q=DZ+\operatorname{diag}(\mathbf r).
$$

Write v=r6 for net desired-product release and w=r15+r16 for modeled benzene production. Instantaneous TOF is v/M-total. The finite-time logarithmic rate control is

$$
X_{rate,j}=Q_{6j}/v.
$$

For the instantaneous two-branch selectivity s=v/(v+w), differentiate ln v minus ln(v+w):

$$
X_{selectivity,j}=\frac{Q_{6j}}{v}-\frac{Q_{6j}+Q_{15j}+Q_{16j}}{v+w}.
$$

Cumulative net product has its own control,

$$
X_{product,j}=\frac{Z_{Pj}}{P(t)-P(0)}.
$$

This also controls average TOF over a fixed positive time. Values requiring a positive flux or product increment at or below the exported 10−20 floor are NaN; negative desired flux has no real logarithmic control. Negative byproduct flux also invalidates the stated selectivity. Missing values must not become zeros in tables or plots.

These are batch-transient derivatives motivated by the transition-state perturbation framework discussed by [Campbell](https://pubs.acs.org/doi/10.1021/acscatal.7b00115). They do not establish steady-state Campbell DRC, and no sum-to-one relation is imposed. The two-flux selectivity excludes bound material and cleavage fragments; it is not a complete experimental carbon selectivity.

## 27. Solver-fixture TOF, control rankings and longevity diagnostics

The completed verification campaign covers eight temperatures—360, 370, 380, 383.15, 390, 400, 410 and 420 K—and five initial free-base equivalents—0.01, 0.025, 0.05, 0.10 and 0.20. Every condition has five hundred reporting points over 3600 seconds. Neutral alcohol activity is fixed at 0.1 as an arbitrary fixture assumption. The resulting forty trajectories and their controls are software tests, with zero physical catalyst predictions.

![Synthetic-rate solver verification with TOF and finite-time control diagnostics](../examples/plots/phase3_solver_verification.png)

**Figure 2. Software fixture only.** Assigned coefficients produce the displayed TOF profiles and finite-time control rankings. These curves verify the ODE and tangent implementations; they are not measured kinetics, computed catalytic performance or evidence for a particular metal ranking.

| Verification quantity | Observed maximum | Acceptance interpretation |
|---|---:|---|
| Jacobian central-difference error | 9.63087×10−11 | Below 10−5 requested bound |
| Metal inventory drift / M | 3.64292×10−17 | Below 10−10 M |
| Base-equivalent drift / M | 4.99600×10−16 | Below 10−10 M |
| Metal sensitivity residual / M | 5.72391×10−17 | Conserved tangent ledger |
| Physical catalyst predictions | 0 | Missing accepted TS network |

At 360 K, changing fixture free-base equivalents from 0.01 to 0.20 changes the terminal fixture TOF from 0.00277494 to 0.00286615 s−1. This numerical illustration is reported only to make the archived grid inspectable. Its sign and size cannot be interpreted as the response of a real catalyst to total added tBuOK because the coefficients and free-base mapping are not physical inputs.

The complete terminal baseline controls below use 383.15 K, 0.05 initial free-base equivalents, assumed a=0.1 and t=3600 s. Every entry is dimensionless and uses the arbitrary rate fixture. Path numbers match Sections 20–23. Scientific catalyst-control rankings remain unavailable.

| Path | Instantaneous rate control | Branch-selectivity control | Cumulative-product control |
|---:|---:|---:|---:|
| 1 | +0.000480414 | +1.56049×10−8 | +0.000365467 |
| 2 | +0.005690930 | +1.18107×10−7 | +0.004427289 |
| 3 | +0.000758250 | +1.04421×10−9 | +0.000749558 |
| 4 | +0.397297003 | −6.81859×10−8 | +0.402873363 |
| 5 | +0.534418278 | −9.18188×10−8 | +0.541870724 |
| 6 | +0.007557573 | −1.31914×10−9 | +0.007686408 |
| 7 | +0.004355047 | +6.31320×10−8 | +0.005078550 |
| 8 | +0.008935764 | +9.84437×10−9 | +0.010886576 |
| 9 | +0.000653097 | +9.46751×10−9 | +0.000761702 |
| 10 | +0.011166962 | +1.23024×10−8 | +0.013606749 |
| 11 | +0.003350912 | +3.69164×10−9 | +0.004082466 |
| 12 | +0.000446788 | +4.92219×10−10 | +0.000544329 |
| 13 | +0.000022339 | +2.46109×10−11 | +0.000027216 |
| 14 | +0.000364505 | +6.86669×10−10 | −0.010873456 |
| 15 | −0.000002785 | −3.94631×10−8 | −0.000002085 |
| 16 | −0.000010157 | −9.65566×10−8 | −0.000008374 |
| 17 | −0.000402091 | +6.76768×10−11 | −0.000212181 |
| 18 | −0.001450854 | +2.44184×10−10 | −0.000794151 |

Within this terminal fixture, imine hydrogenation and imine coordination have the largest positive rate controls. Dimerization's instantaneous and cumulative controls have different signs, illustrating why a terminal flux response is not identical to total product accumulated over the preceding trajectory. The near-zero selectivity values also reflect this particular fixture and observable definition; they are not evidence of experimentally perfect selectivity.

Dimerization can sequester active metal without destroying its elemental inventory; reversible release can return it. CO poisoning and arm cleavage instead remove catalyst from the modeled productive cycle. A longevity diagnostic must therefore distinguish the population available to the on-cycle routes from total conserved metal. The fixture's active-metal fraction is useful for checking this definition, but cannot be called a measured lifetime.

Control rankings can change with time because the state distribution changes. Early aldehyde formation, later imine accumulation and eventual loss of active species can alter which perturbation changes product flux most strongly. Such behavior follows from the coupled ODE and is not automatically a bifurcation. Numerical continuation of defined stationary solutions, or experimentally constrained transient modeling, would be required for that stronger claim. The archived CSV retains all eighteen pathway controls and undefined-value semantics at all reporting times.

## 28. Part IV — Geometric types in a native PyTorch EGNN

The model implements message passing directly with native PyTorch tensors, shared multilayer perceptrons and indexed segment sums. It uses neither PyG nor DGL. Its architectural starting point is the scalar-feature/vector-coordinate construction of [Satorras, Hoogeboom and Welling](https://proceedings.mlr.press/v139/satorras21a.html). The equations and proofs here describe the present bounded-displacement implementation, including its masks, graph contexts, independent heads and matched-pair training.

For each atom i, position x-i has three Cartesian components and is stored in angstrom. Hidden features h-i have F scalar channels. Scalar edge attributes and graph context supply additional invariant information. A feature channel index is not a spatial coordinate index; mixing those two concepts is a common way to introduce a transformation error.

An arbitrary Euclidean transformation has an orthogonal matrix Q and translation t:

$$
x_i^{\prime a}=Q^a{}_b x_i^b+t^a,\qquad Q^TQ=I.
$$

Repeated spatial indices are summed over the three coordinates. Proper rotations have determinant +1. Orthogonal matrices with determinant −1 include reflection. The scalar prediction requirement and the coordinate-output requirement are different:

$$
f(QX+t)=f(X).
$$

$$
F(QX+t)=QF(X)+t.
$$

The first is invariance. The second is equivariance. A scalar activation energy should satisfy the first in an isolated achiral environment; a coordinate array should satisfy the second. The network can satisfy these laws for every weight value, including random or poorly trained weights. Therefore exact symmetry is a property of the representation and implementation, not a statement about chemical accuracy.

Input node features comprise nine element indicators and metal, donor and proton-site flags. Five graph inputs comprise charge, nominal unpaired occupancy divided by four and three state indicators. Nominal GFN2 occupation is not a measured spin splitting. These scalar contexts allow the model to distinguish declared conditions without inventing physical spin-gap labels. Barrier and crossing-gap heads retain independent readiness flags and remain untrained in this campaign.

## 29. Invariant tensor contraction and scalar message construction

Set the relative vector between atoms i and j to

$$
r_{ij}^a=x_i^a-x_j^a.
$$

Translation cancels exactly. Under the transformed coordinates,

$$
r_{ij}^{\prime a}=Q^a{}_b r_{ij}^b.
$$

Contract two vector indices with the Euclidean metric to obtain squared distance:

$$
d_{ij}^2=\delta_{ab}r_{ij}^ar_{ij}^b.
$$

The transformed contraction follows by substitution, then orthogonality:

$$
d_{ij}^{\prime2}=\delta_{ab}Q^a{}_cQ^b{}_d r_{ij}^cr_{ij}^d.
$$

$$
d_{ij}^{\prime2}=(Q^TQ)_{cd}r_{ij}^cr_{ij}^d=\delta_{cd}r_{ij}^cr_{ij}^d=d_{ij}^2.
$$

No determinant enters this identity, so it proves reflection as well as rotation invariance. It is the explicit tensor contraction that supplies the symmetry property; a generic neural network receiving three raw coordinate components would not gain that property merely by being called geometric.

The radial input compresses large separations smoothly:

$$
\rho_{ij}=\ln(1+d_{ij}^2/\ell^2),\qquad \ell=2\ \mathrm{angstrom}.
$$

A shared scalar MLP then constructs messages:

$$
m_{ij}^{\alpha}=\phi_m^{\alpha}(h_i,h_j,\rho_{ij},e_{ij}).
$$

All arguments are invariants, so each output channel is invariant. Linear–SiLU–Linear–SiLU operations mix scalar channels only. Their learned coefficients have no spatial direction, and therefore need not transform when the molecule is rotated.

This proof has a clear assumption: public edge and graph attributes must themselves be scalar invariants. Passing a laboratory-frame x coordinate as an allegedly scalar edge feature would violate that assumption. Likewise, supplied graph connectivity must transform consistently with atom relabeling. Complete graphs satisfy the connectivity condition automatically; sparse graphs require correct edge remapping. Symmetry is guaranteed for well-typed inputs, while the API tests reject malformed edges and cross-graph connections.

The use of squared distances also avoids division by zero for coincident atoms. Their physical plausibility can be checked elsewhere, but their presence should not make the mathematical layer emit an undefined normalized direction.

## 30. Scalar-to-vector lifting and a bounded equivariant displacement

Each incoming message produces a scalar edge coefficient. Let D-i be the larger of one and the number of valid incoming neighbors. The implemented coefficient is

$$
a_{ij}=\eta\frac{\tanh[\phi_x(m_{ij})]}{\sqrt{1+d_{ij}^2/a_0^2}},\qquad \eta=0.05,\quad a_0=1\ \mathrm{angstrom}.
$$

Because the coefficient depends on invariant arguments, it is invariant. Multiplying it by a relative vector lifts the scalar into a polar-vector contribution:

$$
u_i^a=D_i^{-1}\sum_{j\in\mathcal N(i)}r_{ij}^a a_{ij},\qquad x_i^{+,a}=x_i^a+u_i^a.
$$

Substitute the transformation of the relative vector while leaving the coefficient unchanged:

$$
u_i^{\prime a}=D_i^{-1}\sum_jQ^a{}_b r_{ij}^b a_{ij}.
$$

Q is constant with respect to the finite neighbor sum, so

$$
u_i^{\prime a}=Q^a{}_b u_i^b.
$$

Adding the transformed starting coordinate then proves

$$
x_i^{\prime+,a}=Q^a{}_b x_i^b+t^a+Q^a{}_b u_i^b=Q^a{}_b x_i^{+,b}+t^a.
$$

There is no learned vector bias. Such a fixed laboratory-frame vector would introduce a preferred direction and break the proof. The vector direction comes entirely from geometry; the neural network chooses only invariant scalar weights.

The bounded displacement also follows directly. Since absolute tanh is at most one,

$$
\|r_{ij}a_{ij}\|\le\eta\frac{d_{ij}}{\sqrt{1+d_{ij}^2/a_0^2}}\le\eta a_0.
$$

The triangle inequality and division by D-i give a per-layer node displacement no larger than eta times a0. This bounds numerical movement without adding an external directional restraint. For coincident atoms the relative vector is exactly zero and the denominator remains finite, so the contribution vanishes smoothly.

These updated coordinates are latent model coordinates. They have no geometry-optimization supervision and must not be exported as chemically optimized structures merely because their transformation law is correct. Equivariance preserves the relationship between transformed inputs and outputs; it does not establish that the displacement follows a force or lowers a physical energy.

## 31. Induction, permutation invariance, batching and gradient covariance

The hidden update uses a normalized invariant message aggregate:

$$
\overline{m}_i=D_i^{-1}\sum_jm_{ij},\qquad h_i^+=h_i+\phi_h(h_i,\overline{m}_i).
$$

Every argument is scalar invariant, so h-plus is invariant. Sections 29–30 prove that one layer preserves both required types: scalar hidden channels and equivariant coordinates. The initial element/role embedding has the scalar property. Mathematical induction therefore establishes the property for any finite stack of these layers.

For graph g, sum and mean pooling yield

$$
s_g=\sum_{i:b_i=g}h_i,\qquad \mu_g=\frac{s_g}{\max(1,N_g)}.
$$

Each scalar output head applies its own MLP to the sum, mean and invariant graph context. Every argument remains invariant, completing the proof for graph-level output. A missing-label flag affects whether an output is scientifically available; it does not alter this algebraic property.

Now let a permutation pi consistently relabel atoms and edges. Shared MLPs do not contain absolute atom-index inputs, so the relabeled message equals its original counterpart. A finite sum is unchanged by relabeling its terms. Hidden features and coordinates therefore follow the atom permutation, while graph scalars remain unchanged. Native index-add performs the segment sums; floating-point order can change last-bit rounding, which is measured rather than assumed zero.

Packed graphs reject cross-graph edges. Each graph therefore depends only on its own features and coordinates, giving agreement between separate and batched evaluation modulo rounding. Masked padding contributes no messages, hidden state, pooled count or readout. Empty graphs carry a false graph mask; the trained prediction API reports NaN, not a fabricated observed zero.

Gradient covariance follows from the chain rule. Since X-prime=QX+t and X=Q-transpose(X-prime−t),

$$
\frac{\partial f}{\partial x_i^{\prime a}}=\frac{\partial x_i^b}{\partial x_i^{\prime a}}\frac{\partial f}{\partial x_i^b}=Q^a{}_b\frac{\partial f}{\partial x_i^b}.
$$

Thus position gradients rotate as polar vectors. The tests differentiate actual network output under rotation. No force labels supervise this auxiliary model, so the covariance test does not establish physical force accuracy.

## 32. Data ingestion, matched compositions and family-level splitting

The prior campaign supplied 812 actual ALPB(toluene) conformer optimization records. The new audit checks native convergence, retained donor coordination, the complete previous identity gate and a fresh full ligand/ancillary identity review. Geometry hashes, native-output identity, zero return code, normal termination and agreement between parsed native energy and stored eV energy are checked before a record becomes eligible.

The command-level audit checks GFN2, electronic temperature 300 K, accuracy 0.5, charge and nominal occupation. This rejects an easy source of silent label mixing: two rows can have similar filenames while representing different protocols or an altered chemical structure. The audit retains 610 structures and records 202 exclusions with their reasons. These are real optimization records, but they were not individually characterized with Hessians.

There are seventy matched composition/protocol groups, each fixing catalyst, state, formula, atomic-number order, charge, nominal occupancy and protocol identity. Each selects the lowest original conformer index, breaking ties by structure ID without inspecting energy. Every remaining eligible structure is paired with its group's fixed reference, excluding self-pairs. The resulting dataset contains 540 actual electronic-energy differences:

$$
y_{ar}=E_{\mathrm{GFN2/ALPB}}(X_a)-E_{\mathrm{GFN2/ALPB}}(X_r).
$$

Composition, charge and protocol matching make this subtraction interpretable. It is a conformer-relative electronic-model target, not an activation free energy. Nominal occupancy changes from spin-independent native GFN2 do not create an independent spin-gap dataset.

| Backbone and substituent family | Pair count | Split |
|---|---:|---|
| bipyridine_pnnoh / Ph | 52 | Training |
| bipyridine_pnnoh / iPr | 90 | Training |
| pyridine_pnn / Ph | 98 | Training |
| pyridine_pnn / iPr | 134 | Training |
| macho_pnp / iPr | 80 | Validation |
| macho_pnp / Ph | 86 | Test |

Every backbone/substituent family stays in one split across metals, states and conformers. This avoids direct family leakage between candidate/reference pairs. It is not a leave-one-backbone-out evaluation because different substituent families may share a backbone across validation and test. Only six families are present, making the held-out family sample statistically weak. The split protects a defined comparison without claiming broad chemical generalization.

Dataset manifests, rejected sources, fixed-reference selections and geometry files are preserved under `data/phase3/egnn/`. Training and evaluation use the same audited pair list; no missing barrier or spin-gap label is replaced by a numerical zero.

## 33. Pair-potential gauge, normalization and the masked training objective

A shared scalar auxiliary head f-theta assigns a latent potential to each graph. The physical API for this auxiliary task uses only the matched difference

$$
\widehat y_{ar}=f_\theta(X_a)-f_\theta(X_r).
$$

The difference is invariant because each scalar is invariant. It also obeys exact algebraic constraints: swapping candidate and reference negates the prediction, and pairing a geometry with itself yields zero. Adding any shared composition-dependent constant to f-theta leaves the pair prediction unchanged. Absolute potential values are therefore unidentifiable under this training objective and are not exposed as absolute molecular energies.

Only observed training labels determine the normalization mean and scale. The actual values are mu=−0.12956554 eV and sigma=0.42901018 eV. Both prediction and target are normalized in the loss:

$$
z_{ar}^{pred}=\frac{\widehat y_{ar}-\mu}{\sigma},\qquad z_{ar}^{target}=\frac{y_{ar}-\mu}{\sigma}.
$$

Their difference shows that the mean cancels:

$$
z_{ar}^{pred}-z_{ar}^{target}=\frac{\widehat y_{ar}-y_{ar}}{\sigma}.
$$

Inference returns the raw potential difference. Adding the mean to a raw difference would break same-pair zero and antisymmetry, so the API does not do it. Held-out labels never determine normalization; tests alter them and confirm unchanged training statistics.

For observation mask M, the loss is the mean squared normalized residual over observed entries:

$$
\mathcal L=\frac{\sum_{(p,k):M_{pk}=1}(z_{pk}^{pred}-z_{pk}^{target})^2}{\#\{(p,k):M_{pk}=1\}}.
$$

The implementation selects observed entries before arithmetic. Multiplying an undefined residual by zero would still produce NaN. An all-missing batch produces a differentiable empty sum and zero gradient. Unsupported normalizer columns carry a false support mask; placeholder mean/scale values are never target labels.

Two 24-channel layers, twelve node inputs and five graph contexts give 12,675 parameters. Training uses AdamW, learning rate 0.002, weight decay 10−5, norm clipping at ten, two pairs per batch and seed 20260913. Barrier and crossing-gap heads are frozen and unready. The shared representation learns from the auxiliary target alone; that cannot turn an untrained head into a scientific prediction.

## 34. Observed auxiliary benchmark and exact symmetry tests

The CPU run completed nine epochs in 85.032 seconds with two computation threads and peak process RSS 392,376,320 bytes, approximately 374.20 MiB. Eight epochs without validation improvement triggered stopping. The selected checkpoint is epoch zero after its full 374-pair training pass; it is not an untrained initialization. No additional tuning followed inspection of test labels.

| Split | Pairs | Model MAE / eV | Model RMSE / eV | Equal-reference MAE / eV |
|---|---:|---:|---:|---:|
| Training | 374 | 0.28666436 | 0.44369871 | 0.28643789 |
| Validation | 80 | 0.26161685 | 0.52399179 | 0.25939262 |
| Test | 86 | 0.19519977 | 0.38456639 | 0.19399894 |

The equal-reference baseline predicts a zero energy difference. It is a prediction used for comparison, not an invented zero training target. The model is worse than this baseline on the held-out test family by 0.00120083 eV MAE. No predictive improvement is demonstrated. Training execution, retained checkpoints and symmetry compliance are completed results; useful catalyst screening accuracy is not.

![Actual auxiliary EGNN learning and held-out baseline comparison](../examples/plots/phase3_egnn_auxiliary.png)

**Figure 3. Observed auxiliary benchmark.** The network was trained on audited conformer-energy differences. Barrier and MECP-gap heads remain untrained. The held-out baseline comparison must accompany any reported error so that a small-looking number is not mistaken for an improvement.

The actual checkpoint was reloaded in float64 and tested at Euler angles (0.371, −0.893, 1.247) radians with translation (3.17, −2.41, 0.63) angstrom. Maximum scalar discrepancy across rotation, reflection, permutation and batching was 1.7763568394×10−15 eV, below the requested 10−6 eV tolerance. Maximum coordinate-equivariance discrepancy was 1.7763568394×10−15 angstrom. Per-case values and checkpoint hashes are in `symmetry_verification.json`.

These near-machine-precision results are expected from the proved architecture and support correct tensor implementation. They do not reduce the held-out error or supply missing physical labels. With complete graphs, edge count is the sum of N-g times (N-g−1), so dense message work scales approximately as L·E·H² and storage as E·H. The bounded local run does not establish large-scale speed. Reflection-invariant geometry alone cannot distinguish enantiomers without an appropriate chiral environment or additional symmetry-consistent information, and this model contains no explicit spin–orbit coupling tensor.

## 35. Reproduction, integrity checks and the scope of verification

The primary implementation files are `quantum/spin_mecp.py`, `generators/solvation_clusters.py`, `kinetics/master_kinetics.py` and `models/pincer_egnn.py` under `src/pincer_catmech/`. The master runner is `scripts/run_advanced_campaign.py`; specialized runners retain the native outputs and summaries used here. The repository's configured existing environment is reused instead of reinstalling scientific packages or overwriting prior campaign environments.

The microsolvation runner can be invoked from the repository with `python scripts/run_solvation_clusters.py --seeds 3 --threads 2`. In the configured Windows setup the existing wrapper supplies the interpreter. It runs one physical xTB process at a time, bounds numerical-library threading and protects scratch work with a worker lock. Existing native directories are reused only after deterministic geometry, output-checksum and solvent-protocol verification. A new physical campaign should use distinct scratch and output locations.

The portable solvation evidence archive contains 430 files with a byte-level SHA256 manifest. It includes native input/output, optimized coordinates, gradients, restart data, frozen subsets and full Hessian evidence. Original mappings and pre-extension thermochemistry metadata are retained when later diagnostics are added. The three-alcohol failed Hessian classification and repair are both present. This enables a reviewer to distinguish a parser change from a new native calculation.

Focused verification includes 23 solvation tests, 68 independent kinetics/Jacobian tests after review and 16 EGNN tests described in the module receipts. Spin tests include analytical seams, KKT residuals, penalty finite differences, rigid-motion covariance, failed SCF propagation and local minimum versus saddle curvature. Full integration-suite counts belong to the final campaign verification receipt rather than being inferred by summing potentially overlapping focused suites.

Scientific and software validation have different stopping criteria. A finite-difference derivative test asks whether the implemented derivative matches the model. A Hessian asks whether a geometry is stationary and how local curvature behaves under a Hamiltonian. A provenance check asks which artifacts produced a number. An experimental comparison asks whether that model predicts observations. Passing one does not substitute for the others.

The source records and structured summaries are authoritative for full precision and final run status. Printed tables round values for readability. Git deployment, source identity and export checks should refer to the final verified commit and artifact manifest; this document does not declare a push successful merely because a deployment script exists.

### Final integration and portable evidence

The completed full regression run recorded 556 tests, with zero failures or errors, in `data/phase3/verification/pytest.xml`. The tested-source manifest binds every Python module, script and test plus package configuration. Numerical fixtures, mocked process tests and adversarial artifact tests remain software evidence. The physical stage count is taken independently from native receipts.

Summary refresh now binds EGNN checkpoint, dataset, family split, normalization and supported-head flags together. Three exact historical training-source files are preserved in `data/phase3/egnn/source_at_training/`, with separate hashes for current verification code. `python scripts/train_phase3_egnn.py --verify-only` verifies and summarizes the retained run without rebuilding data or training again. The new physical-kinetics intake checks the named design, composition, charge, solvent, grid and base inventory before archiving certified inputs; a valid model for another catalyst cannot be relabeled as a requested design. These gates do not create missing chemical certificates.

The spin publication contains 794 scientific files. The 93 large numbered Psi4 DF-SCF B-matrix and DIIS cache files, totaling 18,136,488,536 bytes, remain unchanged on the originating computer and are omitted from Git/the release ZIP. Their exact path, byte count and SHA256 are listed in `data/phase3/spin/evidence_publication.json`. Inputs, native outputs, gradients and state responses remain published. File codes 97 and 64 are identified in the [official PSIO file map](https://psi4.github.io/psi4docs/master/autodoc_psifiles.html) and the installed Psi4 1.11 header. The release is not a complete scratch backup.

The proton-wire prospective driver, finalizer and replay sequence are documented in `docs/PROTON_WIRE_REPRODUCIBILITY.md`; archived execution-source hashes retain their historical meaning. `python scripts/run_advanced_campaign.py --stage all` audits completed native receipts and reruns all forty software grids. `bash scripts/deploy_phase3.sh` performs the complete verification and no-force publication gate. The final release exporter independently checks the pushed Git tree, every ZIP member and all six reviewed PDFs.

## 36. Decision gates, limitations and primary-source map

The completed campaign supports four concrete conclusions. The crossing solver implements and tests a mathematically defined constrained problem while retaining physical backend failures. Neutral microsolvation produces nine relaxed clusters and three selected characterized minima, with attractive model association energies but positive standard association free energies at the baseline temperature. The eighteen-state kinetic system has an explicit analytical Jacobian, exact material ledgers and verified transient tangents on declared fixtures. The native equivariant network passes its symmetry checks and completes real auxiliary training, but does not improve the simple held-out baseline.

The next chemical decision requires actual accepted transition states and consistent state/reference thermochemistry. Until that dataset exists, no physical TOF ranking or rate-control ranking is available. Open cluster minima cannot fill the dehydration barrier. Small-model MECP diagnostics cannot fill catalyst spin labels. Native occupation-dependent GFN2 outputs cannot fill a validated spin-polarized gap target. The machine-learning readiness API enforces these distinctions rather than returning unsupported numerical predictions.

The following primary-source map identifies the methodological basis of the report. The derivations above are explicit reconstructions for this implementation, and the numerical campaign values originate in repository artifacts rather than in these papers.

| Topic | Primary source | Use in this report |
|---|---|---|
| Crossing surfaces | [Harvey et al., 1998](https://doi.org/10.1007/s002140050309) | Projected crossing-gradient context |
| Constrained seam optimization | [Farazdel and Dupuis, 1991](https://doi.org/10.1002/jcc.540120219) | Historical constrained formulation |
| GFN2-xTB | [Bannwarth et al., 2019](https://doi.org/10.1021/acs.jctc.8b01176) | Electronic-model source |
| ALPB | [Ehlert et al., 2021](https://doi.org/10.1021/acs.jctc.1c00471) | Implicit-solvation source |
| Spin-polarized xTB | [spGFN publication](https://doi.org/10.1002/jcc.27185) | Distinct spin-dependent extension |
| Low-mode thermochemistry | [Grimme, 2012](https://doi.org/10.1002/chem.201200497) | qRRHO motivation |
| Rate control | [Campbell, 2017](https://doi.org/10.1021/acscatal.7b00115) | TS-energy perturbation interpretation |
| Equivariant networks | [Satorras et al., 2021](https://proceedings.mlr.press/v139/satorras21a.html) | Scalar/vector architecture source |

Official executable documentation is linked where used: xTB solvent and full-Hessian conventions, Psi4 SCF and gradients, and SciPy integration interfaces. No source is cited as proof that the present pincer catalyst follows the hypothesized mechanism. A future predictive campaign should close each evidence gate in sequence, preserve failed attempts, compare simple baselines and report uncertainty at the level supported by its actual data.
