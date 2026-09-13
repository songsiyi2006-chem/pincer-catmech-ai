# Neutral micro-solvated proton-wire search: physical negative result

Two new 43-atom hemiaminal/tBuOH seeds and one immutable continuation were evaluated with actual GFN2-xTB, ALPB(toluene), gsolv, charge 0 and no unpaired electrons. One xTB worker and one thread were used throughout. The physical search finished at 15:26:19 UTC on 13 September 2026, before its original 15:26:39 UTC deadline.

**Accepted transition states: 0. Activation free-energy barriers: 0.**

| Physical path run | FIRE/NEB steps | Final NEB fmax, eV/angstrom | Required fmax | Outcome |
| --- | ---: | ---: | ---: | --- |
| Initial seed 02 | 100 | 0.375882 | 0.07 | Unconverged |
| Initial seed 01 | 100 | 1.347282 | 0.07 | Unconverged |
| Immutable seed-02 continuation | 200 additional | 1.206630 | 0.07 | Unconverged |

The two initial runs each passed free reactant and product optimization, exact mapped covalent graphs, C-N bond-length gates, and full native Hessian minimum certification. The continuation repeated both endpoint Hessians and all force checks. These are four initial endpoint certifications plus two repeat certifications, not six independently discovered minima. Because the bands remained unconverged, saddle refinement, a saddle Hessian, mapped imaginary-mode validation and opposite endpoint descents were not reached. No acceptance threshold was relaxed.

The new product preparation fits accepted imine, water and tBuOH component reference geometries and therefore begins with the correct C=N bonding. Every atom retains its original index. Original N-H (atom 17, zero-based) becomes the new alcohol O-H; original tBuOH O-H (atom 42) becomes the second water proton. Original water/leaving oxygen is atom 0; shuttle oxygen is atom 32. The seed-02 geometric product seed contained one H-H close-contact graph flag, preserved in its assembly metadata; unconstrained physical optimization removed that contact and passed the complete product-basin tests. Preparation coordinates are never claimed to be stationary-point evidence.

The model contains C, H, N and O only. It omits a metal catalyst, potassium and tert-butoxide. The experimental boundary of 0.05 equivalents total tBuOK remains a speciation input unresolved by this neutral cluster. Closed six- or eight-membered proton-wire transition states, a base kinetic order and a bulk reaction rate have not been established. The thermodynamic baseline is 383.15 K, while xTB electronic smearing remains 300 K. Endpoint qRRHO tables include 383.15 K and retain the existing documented 1 M standard-state convention and temperature-independent ALPB excess-solvation approximation.

`sampled_band_profiles.csv` contains only actual sampled potential energies and mapped bond distances. A maximum along an unconverged band is not an activation free energy or an accepted saddle barrier. `endpoint_thermochemistry.csv` contains endpoint thermochemistry only. `attempts.csv`, `endpoint_stationarity.csv`, `summary.json` and `verification.json` expose the numerical acceptance boundaries directly.

All native stdout, launch records, gradients, XYZ structures, full Hessians and electronic restart/cache files are retained in the published numbered ZIP archives. `native_manifest.json` qualifies each member by archive name, member path, SHA256 and byte count. Unchanged original native members appear only in the first ZIP; new continuation members appear in subsequent ZIPs. The temporary aggregate is preserved outside the repository in owned scratch. Historical source snapshots and results remain available. Four atom-mapping/model-rejection tests passed; archive and restart checks are documented separately in `verification.json`.
