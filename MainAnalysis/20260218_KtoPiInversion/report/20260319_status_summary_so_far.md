# 2026-03-19 consolidated status summary

This note summarizes the main analysis and documentation work completed up to the current state of the repository.

## 1. Core DELPHI analysis chain

- Active production uses:
  - data: `merged_data_v2.5.root`
  - nominal MC: `merged_pythia_v2.5.root`
  - cross-check MC: `merged_mc_v2.4.root`
- Event selection in the active chain now uses the archived event bit:
  - `PassAll == 1`
- The main `K/pi` analysis chain has been regenerated and propagated to:
  - `Nch_tag`
  - `dN_ch/deta`
  - thrust-axis `dN_ch/dy`

## 2. Thrust-axis `dN_ch/dy` branch

- Added a full parallel `K/pi` branch versus `dN_ch/dy` in `|y_T| < 0.5`
- `y_T` is defined relative to the thrust axis:
  - `y_T = 0.5 * ln((E + p·n_thrust)/(E - p·n_thrust))`
- Added:
  - raw/corrected histograms
  - response matrices
  - unfolding macro and driver
  - systematics finalizer
  - control plots
  - generator comparisons
- Replaced the old coarse `dN/dy` activity binning with a variable-bin scheme:
  - `[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 15.5, 20.5, 25.5, 30.5]`
- Tail merging is now auto-selected in unfolding:
  - `keepBinsOverride = -1`
  - `keepBinsUsed = 14`
- Added dedicated toy-coverage calibration for this branch and switched the quoted statistical bars to toy-calibrated RMSE.

## 3. Standalone generator program

- Installed and/or built local standalone generator studies for:
  - `PYTHIA 8.317`
  - `PYTHIA 8.317 + Ropewalk`
  - `PYTHIA 8.315 + Dire`
  - `PYTHIA 6.428`
  - `HERWIG 7.3.0`
  - `SHERPA 3.0.3`
  - exploratory `X-SCAPE/JETSCAPE` colorless and hybrid controls
- Main standalone statistics were raised to `400k` events for:
  - PYTHIA8 baseline
  - PYTHIA8 Ropewalk
  - PYTHIA8 Dire
  - PYTHIA6
  - HERWIG
  - SHERPA
- Added generator comparison figures for:
  - `K/pi` vs `dN_ch/deta`
  - `K/pi` vs thrust-axis `dN_ch/dy`
  - normalized `dN_ch/deta` shapes
  - normalized thrust-axis `dN_ch/dy` shapes

## 4. Weak-decay daughter veto study

- Studied the DELPHI-style truth definition that excludes charged daughters of ancestors with `c tau > 1 cm`
- Explicitly tested impact on:
  - `K/pi`
  - `dN_ch/deta`
  - thrust-axis `dN_ch/dy`
- Main conclusion:
  - the active DELPHI merged `Gen` record already behaves like the vetoed truth definition
  - therefore the DELPHI merged-truth study shows zero change when the veto is re-applied
- Standalone PYTHIA8 control showed the veto matters physically in a fully inclusive stable-particle record:
  - about `11%` of charged `|eta|<0.5` activity removed
  - about `11%` of counted pions removed
  - kaons nearly unchanged
  - inclusive `K/pi` moved from about `0.135` to about `0.152`

## 5. Standalone-generator alignment to the DELPHI truth policy

- Updated all standalone generator comparisons so that the default counting convention matches the DELPHI-like truth policy:
  - DELPHI-like charged-particle definition
  - weak-decay-daughter veto applied by default
- This was applied to:
  - `K/pi`
  - `dN_ch/deta`
  - thrust-axis `dN_ch/dy`
- Also fixed generator-side technical issues found during this update:
  - replaced the broken generic Herwig HepMC read path with a direct Herwig ASCII parser
  - fixed the parallel `dN/dy` truth-conversion build race

## 6. PID scale-factor correction

- The active resonance-driven PID inputs remain:
  - `SF_K = 0.9101` from `phi -> K+K-`
  - `SF_pi = 0.9571` from `K* -> Kpi`
  - `SF_K / SF_pi = 0.9509`
- The implementation bug was the direction of the correction.
- Old incorrect convention:
  - `R_after = R_before * (SF_K / SF_pi)`
- Current corrected convention:
  - `R_after = R_before / (SF_K / SF_pi)`
- This correction has now been propagated through:
  - `Nch_tag`
  - `dN_ch/deta`
  - thrust-axis `dN_ch/dy`
  - inclusive `e+e-` world comparison
  - before/after PID-SF validation plots
- Updated inclusive corrected DELPHI value:
  - `K/pi = 0.155069 +- 0.000276 (stat) +- 0.010583 (PID SF)`

## 7. Analysis note and paper draft

- The analysis note has been updated repeatedly and now includes:
  - rebinned thrust-axis `dN_ch/dy` branch
  - standalone generator comparisons
  - weak-daughter-veto study
  - corrected PID scale-factor convention
- The paper draft has been updated consistently for the parts already mirrored there.
- Current synchronized repo heads:
  - main repo: `d7d5c57`
  - analysis note Overleaf repo: `d80484e`
  - paper draft Overleaf repo: `8971ec6`

## 8. Status reports and supporting documentation

- Created a detailed standalone-generator slide deck:
  - `report/20260314_generator_standalone_status.pdf`
- Created a delta update slide deck relative to the last reports:
  - `report/20260319_status_update_since_last_report.pdf`
- Supporting markdowns now exist for:
  - v2.5 migration
  - PassAll event-selection update
  - thrust-axis `dN/dy` branch
  - standalone generator installations
  - X-SCAPE/JETSCAPE setup and controls
  - long-lived-daughter veto study
  - standalone generator veto follow-up

## 9. Current interpretation

- The main DELPHI physics conclusion is unchanged:
  - no strong multiplicity-driven rise indicative of a heavy-ion-like strangeness-enhancement signal is seen in this `K/pi` observable
- The thrust-axis `dN_ch/dy` branch is now a useful `e+e-`-native cross-check
- The standalone generator comparisons are now on a more defensible footing because they follow the DELPHI-like weak-daughter-vetoed counting convention
- Main remaining caution:
  - the current `X-SCAPE/JETSCAPE` hybrid result remains exploratory and strongly hadronization-model dependent

## 10. Delivered artifacts worth checking first

- Analysis note:
  - `overleaf_repo_git/main.pdf`
- Paper draft:
  - `paper/overleaf_paper_draft/main.pdf`
- New delta slides:
  - `report/20260319_status_update_since_last_report.pdf`
- Standalone generator slides:
  - `report/20260314_generator_standalone_status.pdf`
