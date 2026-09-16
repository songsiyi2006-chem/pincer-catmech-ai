# Thermodynamically consistent kinetic foundation

This module establishes mathematical consistency gates. It does not establish a new catalytic mechanism, certify a physical rate, or demonstrate readiness for Nature, Science, Cell, or any other journal. All examples in `tests/test_thermodynamic_network.py` are explicitly synthetic mathematical fixtures. The previous 18-species implementation and all archived chemical calculations remain unchanged.

## State variables and exact material ledgers

Let species concentrations be $c_i\geq0$ in mol/L, with one declared concentration standard $c^\circ>0$. For an ideal mixture $a_i=c_i/c^\circ$ is dimensionless. Reaction $r$ has nonnegative integer reactant and product vectors $\alpha_r$ and $\beta_r$, with $\nu_r=\beta_r-\alpha_r$. Matrices have shape `(species, reactions)`.

Each species has an explicit elemental inventory and charge. The constructor checks $B\nu=0$ and $q^T\nu=0$ exactly, using integer arithmetic with overflow protection. An omitted proton, potassium ion, bound counterion or dissociated ligand is not repaired by an implicit reservoir. Actual active species, protonation states, ion pairing and stoichiometry must be supported by experiment and suitable calculations; satisfying an integer ledger alone does not establish their identity.

## Shared energies and reversible flux

Supply species standard chemical potentials $\mu_i^\circ$ and one common transition-state standard chemical potential $G_r^{\ddagger\circ}$ per reversible reaction, all in kcal/mol at the same temperature, solvent model and standard state. The implemented transition-state-theory model has unit transmission coefficient:

$$J_r=c^\circ\frac{k_BT}{h}\left[e^{-(G_r^{\ddagger\circ}-\alpha_r^T\mu^\circ)/(RT)}\prod_i a_i^{\alpha_{ir}}-e^{-(G_r^{\ddagger\circ}-\beta_r^T\mu^\circ)/(RT)}\prod_i a_i^{\beta_{ir}}\right].$$

Here $R=8.31446261815324/4184$ kcal/(mol K); $k_B$ and $h$ use SI units so $k_BT/h$ has units s⁻¹. Consequently $J_r$ has units M/s, and $\dot c=\nu J$. The forward and reverse frequency constants multiplying dimensionless activities have units s⁻¹. Conventional constants multiplying concentration monomials instead have units $M^{1-\sum_i\alpha_{ir}}s^{-1}$ and $M^{1-\sum_i\beta_{ir}}s^{-1}$; the implementation exposes these molarity powers explicitly. A bimolecular coefficient cannot be silently treated as a first-order constant. Symmetry, pathway degeneracy and any justified transmission factors require consistent treatment in the model parameterization; the present core does not invent them.

The input energy arrays carry no physical certificate. Arbitrary supplied values can satisfy every mathematical identity below while remaining chemically wrong. Species and transition-state energies must come from consistent free energies, including justified conformational, thermal and standard-state treatment, before interpreting calculated rates. The IUPAC definition describes the assumptions underlying transition-state theory; normal termination of a quantum calculation is not such validation. [IUPAC transition-state theory](https://goldbook.iupac.org/terms/view/T06470)

## Detailed balance and cycle consistency

Here $k_r^\pm$ specifically denotes the s⁻¹ frequency constants in $J_r^+/c^\circ=k_r^+\prod_i a_i^{\alpha_{ir}}$ and $J_r^-/c^\circ=k_r^-\prod_i a_i^{\beta_{ir}}$. Their ratio is dimensionless. It is not the unnormalized ratio of conventional concentration-power coefficients, whose units differ when the two molecularities differ. Because both frequency constants share the same transition-state energy,

$$\ln(k_r^+/k_r^-)=-\nu_r^T\mu^\circ/(RT)=\ln K_r.$$

$K_r$ is dimensionless when expressed in activities. For any stoichiometric cycle vector $z$ satisfying $\nu z=0$, $z^T\ln K=0$, equivalently $\prod_r K_r^{z_r}=1$. These are the Wegscheider compatibility conditions. They hold algebraically because the equilibrium constants derive from the same species energies; independently fitting unrelated forward and reverse constants can violate them.

At a positive detailed-balance state, $\nu^T\mu=0$ with $\mu=\mu^\circ+RT\ln a$, each reaction flux vanishes. The implementation uses logarithms and a shared transition-state offset to reduce cancellation near equilibrium. Floating-point tests use explicit numerical tolerances; they do not constitute exact-arithmetic proofs.

Changing an elemental or charge energy reference requires changing transition energies too. If $\delta_i=b_i^T\lambda+q_i\lambda_q$, replace $\mu_i^\circ$ by $\mu_i^\circ+\delta_i$ and $G_r^{\ddagger\circ}$ by $G_r^{\ddagger\circ}+\alpha_r^T\delta$. Conservation ensures $\alpha_r^T\delta=\beta_r^T\delta$, leaving barriers, fluxes, affinities and equilibrium constants invariant. Adding an arbitrary identical constant to every species is generally not invariant for reactions changing total molecule number. Similarly, changing $c^\circ$ without the corresponding standard-potential transformation changes the model; the API does not pretend otherwise.

## Ideal closed-system dissipation and its boundaries

For a closed, ideal, isothermal mixture in the positive concentration interior, define the concentration free-energy density (up to irrelevant conserved reference terms)

$$g(c)=\sum_i c_i\{\mu_i^\circ+RT[\ln(c_i/c^\circ)-1]\},\qquad A_r=-\nu_r^T\mu.$$

Its concentration gradient is $\mu$. Consequently

$$\dot g=\mu^T\nu J=-\sum_r J_rA_r\leq0.$$

Writing forward and reverse nonnegative fluxes as $x_r,y_r$, positive-interior kinetics gives $A_r=RT\ln(x_r/y_r)$. Thus $J_rA_r=RT(x_r-y_r)\ln(x_r/y_r)\geq0$. Chemical entropy production is $\sum_rJ_rA_r/T$, in kcal/(L K s). This is the isothermal ideal closed-network diagnostic implemented by `ideal_dissipation`; it is not a proof of globally unique equilibrium or convergence for all networks. Conserved totals, disconnected stoichiometric classes and boundary equilibria still matter. Thermodynamic treatment of open chemical networks requires the corresponding work and exchange terms. [Rao and Esposito, Physical Review X 6, 041064](https://doi.org/10.1103/PhysRevX.6.041064)

With material exchange $u$, $\dot c=\nu J+u$ and $\dot g=-J^TA+\mu^Tu$; the closed-system sign claim no longer applies. Chemostats, evaporation, sampling, gas exchange and feeds must be represented explicitly before applying a physical energy balance. They are not implemented automatically in this core.

`flux_from_activities` and `rhs` accept nonnegative activities/concentrations. They implement $0^0=1$ according to stoichiometric monomials: an absent species contributes no factor when its exponent is zero. At $c_i=0$, every flux direction consuming that species vanishes, so $\dot c_i\geq0$. Tests cover this boundary property. `chemical_potentials`, `affinities`, `ideal_free_energy_density` and `ideal_dissipation` require strictly positive inputs; no floor or clipping is used to manufacture finite chemical potentials. The module supplies a vector field, **not an ODE solver** or a proof that an arbitrary numerical time integrator preserves positivity.

For a nonideal mixture $a_i=\gamma_i(c)c_i/c^\circ$, activity-based fluxes may still be evaluated, but the ideal $g$ is not automatically a Lyapunov function. A compatible excess free energy must provide $RT\ln\gamma_i=\partial g^{\rm ex}/\partial c_i$ under the chosen state constraints, together with integrability and stability conditions. Arbitrary independently chosen $\gamma_i$ do not establish such a potential. The implementation intentionally has no generic nonideal dissipation certificate.

## Coarse graining and spin-changing paths

Eliminating fast intermediates requires demonstrated timescale separation over the relevant conditions, a specified quasi-equilibrium or quasi-steady-state approximation, and validation against the uneliminated network. An equilibrium partition-function free energy $G_{\rm eff}$ alone does not establish effective kinetics: the reduced model must reproduce the relevant exit fluxes and timescales while the internal states equilibrate sufficiently rapidly. Projection can introduce memory, altered effective rate laws, or hidden entropy production; simply deleting a state does not preserve detailed balance or parameter meaning. Reversible dormancy and irreversible loss also require different material and observation models.

An MECP is not an ordinary first-order transition state on a single adiabatic surface. Its energy cannot be inserted into this Eyring expression as a physical rate without an appropriate nonadiabatic treatment, including relevant spin-orbit coupling, state character, crossing dynamics and competing paths. Any reduction to an effective rate requires independent justification. None is supplied here.

## Local observational rank is not global identifiability

`observational_identifiability` takes an already computed mean-observation Jacobian $D_{kj}=\partial m_k/\partial\theta_j$, where $m_k=\mathbb{E}[y_k]$, positive parameter scales $s_j$, and positive observation/noise scales $\sigma_k$. It performs SVD on $S=\operatorname{diag}(\sigma)^{-1}D\operatorname{diag}(s)$. These are diagonal row and column scalings, making the reported singular values comparable under consistent changes of observable and parameter units. Only $s_j=\theta_j>0$ at the evaluation point gives $D\operatorname{diag}(s)=\partial m/\partial\log\theta$; arbitrary positive scales do not imply a logarithmic parameterization. The function returns numerical rank, scaled null directions and a threshold-dependent local condition diagnostic.

The API neither accepts a full noise covariance nor computes or validates a whitening transformation. Correlated errors require externally justified positive-definite covariance $\Sigma$ and a map $W$ with $W\Sigma W^T=I$, consistently applied to observations, mean predictions and their Jacobian. A caller may supply $WD$ with unit row scales after independently checking this transformation. For local coordinates $d\theta=\operatorname{diag}(s)d\xi$, $F_\xi=S^TS$ is the true mean-parameter Fisher information under a Gaussian observation model with known, parameter-independent covariance: either $\Sigma=\operatorname{diag}(\sigma^2)$ or the correctly whitened correlated model. Otherwise it is a local weighted-least-squares/Gauss–Newton information approximation, not a general Fisher identity; covariance derivatives or a non-Gaussian likelihood require additional analysis. The SVD function does not fit a likelihood or itself return a Fisher-information calculation.

No parameter derivatives or experiments are inferred by this function. Full local column rank at one point does not establish global structural identifiability, statistical precision or practical identifiability across the domain. Rank deficiency can reveal locally confounded parameter combinations, but resolving it requires informative measurements or a revised model. A derivative-free network alone supplies no observed parameter identifiability. All included rank examples are synthetic fixtures.

For a second ligand backbone, freeze the rule mapping structure and independently obtained calculations to model parameters, including allowed calibration, before viewing the prediction outcomes. Different structures need not share identical numerical microscopic rates. Re-fitting to the held-out trajectory or endpoint and calling it a blind prediction would invalidate that claim.

## API and evidence status

`ThermodynamicNetwork` accepts `species`, `reactions`, `species_mu0`, `transition_mu0`, `alpha`, `beta`, `compositions`, `charges`, `temperature_K`, and optional `standard_concentration_M=1`. Arrays and elemental mappings are copied into immutable storage. Invalid units are prevented by the API contract and documented dimensional convention; arrays alone do not encode instrument units, so external data conversion must be verified before construction.

Main operations are `flux_from_activities`, `rhs`, `log_equilibrium_constants`, `concentration_rate_constants`, `chemical_potentials`, `affinities`, `ideal_free_energy_density`, `ideal_dissipation`, and `shifted_energy_reference`. `physical_rates_validated` is always false. This module neither enters the existing production certificate gate nor produces catalyst rankings, experimental TOFs, fitted lifetime predictions or new physical labels. Its contribution is a testable mathematical foundation for future evidence-based model construction, not a standalone novelty claim.
