# PYTHIA6 standalone generation status (2026-03-13)

This note records the standalone PYTHIA6 setup used to generate a truth-level
`K/pi` versus `dN_ch/deta` reference for comparison to the DELPHI analysis.

## Version
- Official final PYTHIA 6 release: `6.428`
- Banner string from the executable:
  - `This is PYTHIA version 6.428`
  - `Last date of change: 5 Sep 2013`
- Source file:
  - `external/pythia6428-src/pythia6428.f`
- Download source:
  - `https://pythia.org/download/pythia6/pythia6428.f`

## Physics configuration
The standalone generator is configured as:
- beams: `e+ e-`
- center-of-mass energy: `sqrt(s) = 91.2 GeV`
- hard process: `gamma*/Z` production
  - `MSEL = 0`
  - `MSUB(1) = 1`
- allowed boson decays: hadronic only
  - all `Z -> lepton` channels disabled by looping over the `Z` decay table and
    keeping only quark final states

This should therefore be labeled as:
- `e+e- -> gamma*/Z -> hadrons at sqrt(s)=91.2 GeV`

## Tune / parameter choice
- No explicit PYTHIA6 tune is applied in the standalone driver.
- The run therefore uses the default PYTHIA 6.428 parameter set for this process.
- In figures and summaries this is labeled as:
  - `Standalone PYTHIA 6.428 default`

## Local driver
- Fortran driver:
  - `tools/pythia6_generate_truth_kpi_dndeta.f90`
- Build/run wrapper:
  - `tools/run_pythia6_generate_truth.sh`

## Observable definition
The standalone truth observable matches the DELPHI comparison convention:
- `dN_ch/deta`:
  - charged stable final-state particles with `|eta| < 0.5`
- `K/pi`:
  - charged kaons over charged pions
  - `0.4 < pT < 5.0 GeV/c`

## Output products
- raw event-summary text:
  - `result/20260313/pythia6_truth/pythia6_zpole_truth_200k.txt`
- binned standalone plot:
  - `result/20260313/pythia6_truth/Pythia6Truth_KtoPi_vs_dNdEta.pdf`
- binned summary text:
  - `result/20260313/pythia6_truth/Pythia6Truth_KtoPi_vs_dNdEta.txt`
