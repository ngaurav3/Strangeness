# Analysis Note Consistency Status (2026-03-09)

**Author:** Yen-Jie Lee and OpenAI

## Scope
- Review the current analysis note source and its build path.
- Cross-check multiplicity definitions, PID-scale-factor documentation, systematic-uncertainty narrative, and figure references against the active v3 code.
- Record what has already been completed and what remains to be edited.

## Note Source Located
The current PDF-building sources were identified at:
- `overleaf_repo_git/main.tex`
- `overleaf_repo_git/StrangenessV3/analysis_note_v3_source.tex`
- `result/20260306/notes/20260306_analysis_note_v3.tex`

The active compiled note is built from:
- `overleaf_repo_git/main.tex`

## Build and Workflow Files Confirmed
The current v3 workflow files present in the analysis area are:
- `analysis.sh`
- `run_variations_20260306.sh`
- `runNtagUnfolding_BayesSVD.C`
- `runDNdEtaUnfolding_BayesSVD.C`
- `build_systematics_v3_with_unfolding.py`
- `finalize_dndeta_systematics.py`

## What Has Been Checked So Far
### 1. Section 3.1 figure-reference audit
- Located the paragraph starting with `Discussion of Figs.~1--2.` in `overleaf_repo_git/main.tex`.
- Verified that Figure 1 is the workflow schematic, so the paragraph is off by one and should refer to Figures 2--3 instead.

### 2. Multiplicity-definition audit
Reviewed the active v3 code path in:
- `KtoPiAnalysis.cpp`

Current implemented reco-level `NchTag` definition in the active v3 code:
- `RecoGoodTrack == 1`
- `RecoCharge != 0`
- `pT > NtagPtMin` with nominal `NtagPtMin = 0.2 GeV/c`
- optional `|eta| < 0.5` only when `UseCentralEtaNtag=true`
- no PID-tag requirement enters `NchTag`

Current implemented truth-level `NchTagTrue` definition in the active v3 code:
- stable charged particles only (`GenStatus == 1`, charged PDG)
- `pT > NtagPtMin`
- optional `|eta| < 0.5` only when `UseCentralEtaNtag=true`

Current implemented truth-level `dNch/deta(|eta|<0.5)` definition:
- stable charged particles only
- `|eta| < 0.5`
- no extra `pT` threshold

### 3. Legacy/public prototype conflict identified
An older prototype exists at:
- `MainAnalysis/20260213_KtoPi/KtoPiAnalysis.cpp`

This older prototype defines `NchTag` differently:
- it counts PID-tagged tracks rather than all good charged tracks above threshold.

Conclusion:
- the current note should stay aligned to the active v3 production chain in `20260218_KtoPiInversion/`, not the older prototype.
- if the older prototype is kept in the repo, the note should explicitly state that it is not the code path used for the quoted v3 results.

### 4. PID scale-factor consistency check
Observed inconsistency to resolve in the note:
- Section 4 currently contains a standalone `phi -> K+K-` scale-factor discussion.
- Section 6 applies final `SF_K / SF_pi` values in the main result chain.

Numerical values currently hard-coded in the final-result/systematics scripts:
- `SF_K = 0.893777 +- 0.043210`
- `SF_pi = 1.019310 +- 0.093383`

Additional context found in the embedded PID note material:
- `phi -> K+K-` summary uses the same kaon central value with a decomposed uncertainty
- `K0_S -> pi+pi-` summary exists in the older embedded PID note material

Implication:
- the main note needs an explicit summary table that distinguishes:
  - channel-specific calibration studies,
  - the final v3-applied scale factors,
  - and which entries are used vs not used.

### 5. Systematics narrative check
Current note issue identified:
- Section 6 opening says the budget has three components.
- The equation and the figures already include PID scale-factor propagation as a fourth term.

Required documentation fix:
- rewrite the opening sentence to four components:
  - binning
  - residual non-closure
  - unfolding
  - PID SF propagation

### 6. Figure-caption and axis-label audit
Already confirmed as needing harmonization:
- the note text mixes `N_{tag}^{ch}`, `N_{tag}^{ch}`, `Ntag`, and `Nch_tag`
- several plot generators still use `N_{tag}^{ch}` in rendered axis titles
- the requested preferred notation is `Nch_tag` in text and `N_{ch}^{tag}` in LaTeX/math

Main plot generators identified for updates:
- `plotClosure.C`
- `plotResultCor.C`
- `plotClosureSpectra.C`
- `runNtagUnfolding_BayesSVD.C`
- `runNtagUnfolding_PoverPi_BayesSVD.C`
- `plotSpeciesResponseMatrices.C`
- `study_bayes_iterations_ntag.py`
- `build_systematics_v3_with_unfolding.py`
- `make_kpi_pid_sf_before_after_plots.py`
- `make_ptopi_results.py`
- `make_systematics_component_ratio_plots.py`

### 7. Reproducibility appendix audit
Current appendix is too prose-heavy for the requested style.
The repo already contains the core runnable files needed for a short shell snippet:
- `analysis.sh`
- `run_variations_20260306.sh`
- post-processing scripts in the same directory

## Work Already Completed Before This Consistency Pass
The note-facing figure style pass has already been applied and pushed earlier:
- multi-panel ratio pad cleanup
- stronger ratio reference lines
- response-matrix palette/margin cleanup
- iteration-scan legend cleanup
- final comparison-plot legend cleanup

Overleaf was last pushed with these figure-style updates at:
- `overleaf_repo_git` commit `972ed30`

## Remaining Note Edits Identified
The following edits are still pending in the note source:
- fix Section 3.1 figure references from Figures 1--2 to Figures 2--3
- add a boxed multiplicity-definition block after the baseline configuration subsection
- standardize notation to `Nch_tag` / `N_{ch}^{tag}` throughout the note and figure captions
- add a Section 4 summary table for PID scale factors used in v3
- reconcile Section 6 text, equations, and table with PID SF propagation
- tighten the reproducibility appendix into a short command snippet with commit hash
- clean the implementation-consistency audit wording and remove dangling phrasing
- improve standalone readability of key captions with correction-stage and `pT`-range information

## Recommended Next Execution Order
1. Edit `overleaf_repo_git/main.tex` first.
2. Regenerate figure PDFs after axis-label changes.
3. Sync mirrored note source copies.
4. Rebuild `overleaf_repo_git/main.pdf`.
5. Run a final consistency sweep between the note and:
   - `KtoPiAnalysis.cpp`
   - `analysis.sh`
   - `run_variations_20260306.sh`
   - unfolding/systematics post-processing scripts.
