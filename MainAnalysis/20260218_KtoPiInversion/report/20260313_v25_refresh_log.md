# v2.5 Refresh Log

Date: 2026-03-13
Area: `MainAnalysis/20260218_KtoPiInversion`

## Scope
This update switches the active nominal analysis from the v2.4 nominal MC/data pair to the newly available v2.5 nominal sample set and propagates the change through the current analysis note and paper draft.

## Sample configuration used
Active nominal inputs after this refresh:
- Nominal MC: `sample/Strangeness/merged_pythia_v2.5.root`
- Data: `sample/Strangeness/merged_data_v2.5.root`

Cross-check MC retained for now:
- Cross-check MC: `sample/Strangeness/merged_mc_v2.4.root`

Reason:
- `merged_pythia_v2.5.root` and `merged_data_v2.5.root` are available.
- `merged_mc_v2.5.root` was not found in the current sample area, so the PYTHIA 5.7 / JETSET 7.4 cross-check remains on v2.4.

Located v2.5 files:
- `/mnt/data2/data2/chenyi/Strangeness/merged_data_v2.5.root`
- `/mnt/data2/data2/chenyi/Strangeness/merged_pythia_v2.5.root`

## Code / configuration updates
Updated active defaults and drivers to v2.5 nominal:
- `KtoPiAnalysis.cpp`
- `BuildDNdEtaResponse.cpp`
- `analysis.sh`
- `analysis_mt.sh`
- `run_systematics_checks.py`
- `run_systematics_dndeta_unfolding.py`
- `run_variations_20260306.sh`
- `make_basic_control_plots.C`
- `make_kpid_pipid_correlation.C`
- `make_kpid_pipid_correlation_plots.py`
- `run_v5_duplicate_pid.sh`

Updated user-facing labels so the nominal branch is no longer hardcoded as "PYTHIA 8 v2.4":
- `make_kpi_vs_dndeta_comparison_plots.py`
- `make_ptopi_vs_dndeta_comparison_plots.py`

Documentation text updated:
- `overleaf_repo_git/main.tex`
- `paper/overleaf_paper_draft/sections/analysis.tex`
- `paper/overleaf_paper_draft/sections/results_and_discussion.tex`
- `paper/overleaf_paper_draft/sections/summary.tex`
- `report/20260312_complete_status_report_from_overleaf_note.tex`
- `report/build_status_report_native_pptx.py`

## Rerun summary
### Full nominal chain rerun
Executed the active nominal production with v2.5:
- `output/KtoPi-MC-Gen-Closure.root`
- `output/KtoPi-MC-Reco-Closure.root`
- `output/KtoPi-MC-Gen-Nominal.root`
- `output/KtoPi-MC-Reco-Nominal.root`
- `output/KtoPi-Data-Reco-Nominal.root`

The v2.5 nominal pass completed cleanly. No NaN/corruption issue reappeared in the core corrected spectra.

### Targeted quoted systematics refresh
Refreshed the currently quoted binning variations for the active systematic model:
- `output/systematics_20260306/nominal_data.root`
- `output/systematics_20260306/nominal_mc.root`
- `output/systematics_20260306/npt10_data.root`
- `output/systematics_20260306/npt10_mc.root`
- `output/systematics_20260306/npt14_data.root`
- `output/systematics_20260306/npt14_mc.root`

### Unfolding products refreshed
Refreshed nominal and binning-variation unfolding products used by the note figures:
- `output/NtagUnfolding_BayesSVD.root`
- `output/PtoPiUnfolding_BayesSVD.root`
- `output/systematics_20260306_dndeta/nominal_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/npt10_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/npt14_unfold_dndeta.root`
- `output/systematics_20260306_dndeta/nominal_unfold_dndeta_ptopi.root`
- `output/systematics_20260306_dndeta/npt10_unfold_dndeta_ptopi.root`
- `output/systematics_20260306_dndeta/npt14_unfold_dndeta_ptopi.root`

### Summary builders refreshed
Rebuilt the current quoted summaries:
- `output/systematics_20260306_v3/systematics_table_v3.txt`
- `output/systematics_20260306_v3/Systematics_ByNtag_Table_v3.tex`
- `output/systematics_20260306_v3/KtoPi_DoubleRatio_with_Systematics.pdf`
- `output/systematics_20260306_v3/Systematics_v3_Components_All.pdf`
- `output/systematics_20260306_v3/Systematics_Components_vs_Ntag.pdf`
- `output/systematics_20260306_dndeta/systematics_dndeta_summary.root`
- `output/systematics_20260306_dndeta/systematics_dndeta_table.txt`
- `output/systematics_20260306_dndeta/KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf`
- `output/systematics_20260306_dndeta/Systematics_Components_vs_dNdEta.pdf`

## Figure regeneration
Regenerated the note-facing figure payload, including:
- basic control plots
- PID correlation plots
- closure / corrected-ratio overlays
- Ntag and dN/deta unfolding diagnostics
- Bayesian iteration scans
- reviewer follow-up validation plots
- K/pi PID scale-factor before/after plots
- K/pi and p/pi final result plots
- DELPHI vs ALICE comparison plots
- DELPHI vs historical e+e- comparison plots
- inclusive K/pi vs sqrt(s) plot

Representative refreshed assets:
- `overleaf_repo_git/StrangenessV3/Generated/KtoPi_DoubleRatio_with_Systematics.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_WorldEE.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/ControlPlot_NchTag_DataMC.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/ControlPlot_DNdEtaEta05_DataMC.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/ControlPlot_TrackKinematics_DataMC.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KPID_vs_PiPID_Correlation_MCData.pdf`

Paper draft figures refreshed:
- `paper/overleaf_paper_draft/figures/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf`
- `paper/overleaf_paper_draft/figures/KtoPi_vs_dNdEta_DELPHI_vs_WorldEE.pdf`

## Documentation rebuild status
Rebuilt successfully after syncing the v2.5 figure payload:
- Analysis note: `overleaf_repo_git/main.pdf`
- Paper draft: `paper/overleaf_paper_draft/main.pdf`

## Overleaf push status
Analysis note repo:
- commit pushed: `3bd04cd`

Paper draft repo:
- commit pushed: `61e4c34`

## Important caveat
The current quoted analysis is now v2.5 for the nominal MC/data chain, but the PYTHIA 5.7 / JETSET 7.4 cross-check still points to `merged_mc_v2.4.root` because a `merged_mc_v2.5.root` file was not available in the sample area at the time of this refresh.

## Non-quoted legacy variation files
The older acceptance/extrapolation variation files in `output/systematics_20260306/` were not exhaustively regenerated in this refresh because they are not part of the currently quoted final systematic budget. The active quoted summary uses:
- binning
- residual non-closure
- unfolding
- PID scale-factor propagation

## Bottom line
The active v6 chain in the analysis note and paper draft now uses the v2.5 nominal sample pair:
- `merged_pythia_v2.5.root`
- `merged_data_v2.5.root`

and the refreshed figures and compiled PDFs have been propagated to both Overleaf projects.
