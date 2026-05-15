2026-03-19
Iteration 14 Response

# Summary

Addressed the reviewer concern by making the analysis-note workflow figure set explicitly v7-driven in the main body of the note.

# Changes Made

1. Replaced the legacy v6 workflow schematic with a workflow-v7 schematic.
   - The new schematic now follows the quoted Yi-style order:
     - reco observed tags in flat `mu`
     - fake correction
     - flat reco-`mu` to true-`mu` unfolding
     - truth-cell `3x3` inversion
     - truth-cell matched/gen efficiency correction
     - projection to final `K/pi` plus PID-SF and inherited systematics

2. Removed legacy v6 figure blocks from the main note body.
   - Removed the legacy pre-unfolding closure package.
   - Removed the legacy `N_tag` unfolding / refolding / Bayes-iteration / toy-coverage figure package.
   - Removed the legacy thrust-axis Bayes/SVD diagnostic package.
   - Removed the exploratory `K/pi` response-matrix / binning-study block that was not part of workflow v7.

3. Removed the legacy `N_tag` result figure from the main Results section.
   - The main Results section now focuses on the quoted v7 observables only:
     - `K/pi` vs `dN_ch/deta`
     - `K/pi` vs thrust-axis `dN_ch/dy`

4. Updated the systematic-summary figure pair.
   - Left panel remains `Systematics_Components_vs_dNdEta.pdf`.
   - Right panel is now `Systematics_Components_vs_dNdY.pdf`.
   - The caption now states explicitly that both panels are inherited legacy-v6 fractional envelopes mapped onto the quoted v7 central values.

5. Removed remaining legacy `N_tag` figure usage from the Systematics section.
   - Removed the legacy `N_tag` systematics-component figure.
   - Removed the `Systematics_ByNtag_Table_v3.tex` inclusion from the main note.
   - Reduced the PID-SF before/after figure to the workflow-v7 `dN_ch/deta` panel only.

6. Updated the surrounding method text.
   - The main workflow description now presents v7 as the active quoted method.
   - Legacy v6 is described only as the source of the inherited systematic envelope and as the comparison baseline used in the explicit v7-vs-v6 figures.

# Build Status

- Rebuilt successfully:
  - `overleaf_repo_git/main.pdf`

# Scope Note

This response updates the workflow-defining figure set in the main analysis note.
Two classes of figures remain by design:

- the resonance-fit PID appendix figures in `DataMCSF_v5.4/`
- the step-by-step independent-cross-check material in `IndependentCrossCheck/`

These are source-documentation / appendix figures rather than legacy-v6 central-result figures, so they were not removed as part of this workflow-v7 figure cleanup.

# Bottom Line

The main analysis-note figure set is now workflow-v7 consistent:
- the quoted result figures are v7,
- the main workflow schematic is v7,
- the remaining main-body systematic-summary figures are tied to the quoted v7 observables,
- and the legacy v6 workflow figure package has been removed from the main note body.
