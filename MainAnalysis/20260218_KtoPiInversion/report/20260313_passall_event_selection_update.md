# 2026-03-13 PassAll event-selection update

## What changed
- The active nominal event preselection in the main analysis code was switched from recomputed cuts to the archived event-selection bit `PassAll` stored in the DELPHI open-data trees.
- The main reco/gen analysis now uses `PassAll == 1` by default.
- The old recomputed selection was kept only as a fallback option in the main executable through `UsePassAllSelection=false` for debugging and compatibility studies.

## Why
- The tree already carries the archived preselection bits `PassNch`, `PassThrust`, `PassTotalE`, and `PassAll`.
- Using `PassAll` makes the analysis follow the sample-defined event selection directly instead of assuming that a local recomputation from `sum(RecoE)`, `Nch`, and `ThrustZ` is exactly equivalent.

## Affected code paths
- Main analysis executable:
  - `KtoPiAnalysis.cpp`
- Standalone dN/deta response builder:
  - `BuildDNdEtaResponse.cpp`
- Basic data/MC control-plot macro:
  - `make_basic_control_plots.C`
- p/K/pi closure-vs-theta diagnostic:
  - `make_closure_vs_pt_theta.C`
- eta-closure diagnostic:
  - `make_closure_vs_eta.py`

## Implementation detail
- New nominal behavior in `KtoPiAnalysis.cpp`:
  - `UsePassAllSelection=true` by default
  - event accepted iff `M->PassAll == 1`
- Legacy fallback remains available:
  - `UsePassAllSelection=false`
  - then the code recomputes the old cuts from `sumRecoE / EcmRef > 0.5`, `Nch >= 7`, and `30 deg < acos(ThrustZ) < 150 deg`

## Important status note
- This step changes the nominal event-selection definition in code.
- The full v2.5 production outputs and note figures have not yet been regenerated under the new `PassAll`-based selection in this step.
- Until that rerun is done, existing physics outputs should be regarded as produced with the previous recomputed-cut implementation.

## Recommendation
- Regenerate the nominal reco/gen outputs, unfolding products, control plots, and note figures under the `PassAll` selection before treating this as the new frozen baseline.
