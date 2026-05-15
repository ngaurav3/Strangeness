# 2026-03-19 v7 Analysis-Note Plot Audit

## Scope
I audited every plot file referenced by `overleaf_repo_git/main.tex` through `\includegraphics` and checked whether it is:
- updated to the current v7 Yi-style note state,
- intentionally kept as a legacy/reference plot,
- or inconsistent with the rewritten note text.

I also rebuilt the note-facing Yi comparison figures where needed and verified the figure text directly with `pdftotext`.

## Main conclusion
The note is now internally consistent at the figure level for the promoted v7 result.

The only real inconsistency I found was in the two Yi comparison figure pairs:
- `YiIndependent_DNdEta_Comparison.pdf`
- `YiIndependent_DNdEta_ClosureComparison.pdf`
- `YiIndependent_DNdY_Comparison.pdf`
- `YiIndependent_DNdY_ClosureComparison.pdf`

Their embedded legend text was still using the old wording:
- `Nominal result`
- `Nominal syst.`
- `Yi-style independent`
- `Yi / nominal`
- `Nominal closure`

That was inconsistent with the rewritten note text, which now presents the result as `current v7` compared to `legacy v6`.

I regenerated those four figures and re-synced them into `overleaf_repo_git/StrangenessV3/Figures/`.

Verified new embedded figure text:
- `Legacy v6 result`
- `Legacy v6 syst.`
- `Current v7 Yi-style`
- `v7 / v6`
- `Legacy v6 closure`
- `Current v7 closure`

## Plots updated and verified as v7-facing
These are the main plots that now match the current v7 note narrative and central values.

### Beam-axis `dN_ch/deta` result plots
- `StrangenessV3/Figures/KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_WorldEE.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_WorldEE_NoPIDSFOnly.pdf`
- `StrangenessV3/Figures/KtoPi_PIDScaleFactor_BeforeAfter_dNdEta.pdf`

### Thrust-axis `dN_ch/dy` result plots
- `StrangenessV3/Figures/KtoPi_DoubleRatio_vs_dNdY_with_Systematics.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdY_DataMC_with_Systematics.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdY_DELPHI_vs_ALICE.pdf`
- `StrangenessV3/Figures/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf`

### Legacy-v6 to current-v7 comparison figures
- `StrangenessV3/Figures/YiIndependent_DNdEta_Comparison.pdf`
- `StrangenessV3/Figures/YiIndependent_DNdEta_ClosureComparison.pdf`
- `StrangenessV3/Figures/YiIndependent_DNdY_Comparison.pdf`
- `StrangenessV3/Figures/YiIndependent_DNdY_ClosureComparison.pdf`

## Plots intentionally kept as legacy/reference
These were not part of the v7 central-value promotion and were left unchanged on purpose.

### Legacy `N_tag` result and validation package
- `StrangenessV3/Generated/KtoPi_DoubleRatio_with_Systematics.pdf`
- `StrangenessV3/Figures/NtagUnfolding_ResponseMatrix.pdf`
- `StrangenessV3/Figures/NtagUnfolding_ResponseMatrix_BySpecies.pdf`
- `StrangenessV3/Figures/NtagUnfolding_MethodDifference.pdf`
- `StrangenessV3/Figures/NtagUnfolding_MCClosure_BayesVsSVD.pdf`
- `StrangenessV3/Figures/NtagUnfolding_DataMC_BayesVsSVD.pdf`
- `StrangenessV3/Figures/NtagUnfolding_RefoldingValidation.pdf`
- `StrangenessV3/Figures/KtoPi_PIDScaleFactor_BeforeAfter_Ntag.pdf`
- `StrangenessV3/Generated/Systematics_v3_Components_All.pdf`
- `StrangenessV3/Generated/Systematics_ByNtag_Table_v3.tex`

These are still used in the note as legacy v6 reference material. That is intentional and now stated in the text.

### Legacy/systematics validation retained for reference
- `StrangenessV3/Figures/DNdEtaUnfolding_RefoldingValidation.pdf`
- `StrangenessV3/Figures/DNdEtaUnfolding_RatioRefoldingValidation.pdf`
- `StrangenessV3/Figures/DNdEtaKeepBinsScan_Summary.pdf`
- `StrangenessV3/Figures/DNdEtaKeepBinsScan_SpeciesClosure.pdf`
- `StrangenessV3/Figures/UnfoldingToyCoverage_Summary.pdf`
- `StrangenessV3/Figures/DNdYToyCoverage_Summary.pdf`
- `StrangenessV3/Figures/DNdYUnfolding_ResponseMatrix.pdf`
- `StrangenessV3/Figures/DNdYUnfolding_MethodDifference.pdf`
- `StrangenessV3/Figures/DNdYUnfolding_MCClosure_BayesVsSVD.pdf`
- `StrangenessV3/Figures/DNdYUnfolding_DataMC_BayesVsSVD.pdf`
- `StrangenessV3/Figures/Systematics_Components_vs_dNdEta.pdf`

These figures still document the validated legacy uncertainty and stability package. The note now explains that the current v7 `dN_ch/deta` and `dN_ch/dy` bands inherit the validated legacy-v6 systematic envelope, rather than coming from a fully separate frozen v7 variation sweep.

### Standalone generator truth-shape plots
- `StrangenessV3/Figures/Generator_dNdEta_Comparison.pdf`
- `StrangenessV3/Figures/Generator_dNdY_Comparison.pdf`

These were not regenerated in the v7 promotion because they are standalone truth-level activity-shape comparisons and do not depend on the choice of v6 vs v7 DELPHI data correction chain.

### Other control / appendix plots intentionally unchanged
- control plots (`ControlPlot_*`)
- closure plots (`KtoPiClosure_Overlay`, `KPtClosure_Overlay`, `PiPtClosure_Overlay`, `PPtClosure_Overlay`, `Closure_vs_Theta_KPiP_Panels`)
- weak-decay-veto appendix figures
- `MCEfficiency` figures
- `DataMCSF_v5.4` figures
- `IndependentCrossCheck` embedded material

These are not part of the v7 central-value promotion itself.

## Additional audit findings
- No `p/pi` result plots remain referenced from `main.tex`.
- The note title/date and the main `Results`, `Systematic Uncertainties`, `Conclusions`, and appendix Yi sections are now consistent with the v7 narrative.
- The current note correctly distinguishes:
  - current quoted v7 activity-based results,
  - retained legacy v6 plots,
  - and standalone truth/generator context plots.

## Decision
The note is in a consistent state after the Yi comparison figure refresh.

What is now true:
- the central `dN_ch/deta` and `dN_ch/dy` result figures are v7,
- the main-text comparison figures now use `v7` vs `v6` language consistently,
- legacy validation/systematics plots remain present but are explicitly described as legacy/reference material,
- no `p/pi` result block remains in the note.

## Files changed during this audit
- regenerated and resynced:
  - `overleaf_repo_git/StrangenessV3/Figures/YiIndependent_DNdEta_Comparison.pdf`
  - `overleaf_repo_git/StrangenessV3/Figures/YiIndependent_DNdEta_ClosureComparison.pdf`
  - `overleaf_repo_git/StrangenessV3/Figures/YiIndependent_DNdY_Comparison.pdf`
  - `overleaf_repo_git/StrangenessV3/Figures/YiIndependent_DNdY_ClosureComparison.pdf`

## Remaining limitation
The current v7 activity-result systematic bands are still inherited from the validated legacy-v6 envelope. That is documented correctly in the note, but it remains the main technical limitation of the present v7 status.
