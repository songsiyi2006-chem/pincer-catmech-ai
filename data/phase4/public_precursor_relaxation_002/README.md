# Public Mn-PNP precursor: native xTB preparation diagnostic

Two crystallographically independent molecules of **the same precursor** were
optimized with GFN2-xTB/ALPB(toluene), followed by fresh gradients and complete
Cartesian Hessians. See `summary.json` for exact inputs, resource limits,
results and hashes; each molecule has a complete `native_manifest.json`.
Both endpoints satisfy the model minimum checks. This does not establish two
distinct conformational basins, ground spin, solution-active identity, DFT
accuracy, a transition state, or a physical catalytic rate. No reaction free
energy was computed. The 300 K electronic smearing parameter is not the
383.15 K reaction-temperature context.

The earlier `public_precursor_relaxation_001` record was stopped by a resource
precheck and started **zero** native calculations. This separate run started
six calls: two optimizations, two fresh gradients and two Hessians. Hessian
internal evaluations are not separately counted as independent structures.

## Attribution and data license

Input and relaxed coordinates derive from Elangovan et al., *Selective
Catalytic Hydrogenations of Nitriles, Ketones, and Aldehydes by Well-Defined
Manganese Pincer Complexes*, JACS (2016),
[10.1021/jacs.6b03709](https://doi.org/10.1021/jacs.6b03709), supporting data
[10.1021/jacs.6b03709.s002](https://doi.org/10.1021/jacs.6b03709.s002).
The original dataset and coordinate derivatives retain
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/), not the
repository's MIT software license. Preserve the detailed
[source attribution](../public_structure/import_v002/THIRD_PARTY_LICENSE.md).
Changes in this directory additionally include native quantum-model geometry
optimization and the associated gradient/Hessian evaluations. These new
calculations are this project's model results, not new experimental
observations by the original authors. Original CIF coordinates remain frozen
in `public_structure/import_v002`.
