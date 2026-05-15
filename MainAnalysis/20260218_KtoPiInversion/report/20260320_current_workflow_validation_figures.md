# 2026-03-20 current-workflow validation figures

## Scope
Added current-workflow validation material back into the analysis note for the two quoted observables:
- `K/pi` vs `dN_ch/deta`
- `K/pi` vs thrust-axis `dN_ch/dy`

## What was added
Two new response-matrix figures were produced from the current flat-`mu` inputs:
- `CurrentWorkflow_DNdEta_ResponseMatrix.pdf`
- `CurrentWorkflow_DNdY_ResponseMatrix.pdf`

Each figure shows an activity-projected view of the current flat-`mu` response for the three observed tag categories:
- `K`-tag
- `pi`-tag
- `p`-tag

The projection is done by summing over the coarse momentum and `|cos(theta)|` groups used in the current solver. Each true-activity row is normalized to unity so the migration pattern along the activity axis is visible.

## Existing closure figures reused
The current-workflow closure figures already existed and were reused in the main note body:
- `YiIndependent_DNdEta_ClosureComparison.pdf`
- `YiIndependent_DNdY_ClosureComparison.pdf`

These are the closure comparisons for the current branch against the legacy factorized branch.

## Note changes
Updated:
- `overleaf_repo_git/main.tex`

Main-text changes:
- in the `dN_ch/deta` results section, replaced the duplicated current-versus-legacy comparison block with a validation block containing:
  - latest-workflow response matrix
  - latest-workflow closure comparison
- in the thrust-axis `dN_ch/dy` results section, added a new validation block containing:
  - latest-workflow response matrix
  - latest-workflow closure comparison

## Figure assets synced to the note repo
- `overleaf_repo_git/StrangenessV3/Figures/CurrentWorkflow_DNdEta_ResponseMatrix.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/CurrentWorkflow_DNdY_ResponseMatrix.pdf`

## Build status
- `overleaf_repo_git/main.pdf` rebuilt successfully after the update.
