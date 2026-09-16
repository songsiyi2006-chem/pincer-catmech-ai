# Offline ORCA diagnostic bundle

Status: prepared_not_submitted. No ORCA calculation is certified by preparation.

Protocol: r2SCAN-3c / CPCM(toluene) / TightSCF / single point; ORCA 6.x.
383.15 K is a reaction reference, not a computed thermochemical correction.
No TS, MECP, solution active species, or kinetic certificate is created.

Run `python bridge.py verify --bundle .` after transfer.
Site inputs: HPC_ACCOUNT, HPC_PARTITION, ORCA_EXE, HPC_PYTHON, HPC_SCRATCH.
A licensed site installation and administrator-reviewed MPI environment are required.
Use `python bridge.py submit --bundle . --mode pilot` to inspect the command.
Only adding `--execute` calls sbatch; start with one pilot and inspect native output.
Production requires a bound human pilot_review.json; see docs/HPC_BRIDGE_ZH.md in the repository.

## Per-structure origin

- MnPNP_complex1_Mn1_M1: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
- MnPNP_complex1_Mn1_M3: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
- MnPNP_complex1_Mn1_M5: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
- MnPNP_complex1_Mn2_M1: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
- MnPNP_complex1_Mn2_M3: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
- MnPNP_complex1_Mn2_M5: public_crystal_derived_candidate; public crystal-derived candidate; solution identity and ground spin unverified.
  Source metadata (not independently authenticated by this bridge): {"source_url": "https://ndownloader.figshare.com/files/5467289", "source_identifier": "doi:10.1021/jacs.6b03709.s002", "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"}
