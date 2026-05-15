# 2026-03-13 closure investigation vs pT and theta

## Scope

This study extends the existing MC closure checks for identified hadrons by examining species-level closure for charged kaons, pions, and protons as a function of:

- `pT`
- track polar angle `theta_track`

The nominal v2.5 sample and the nominal reco correction logic were used:

- input MC: `sample/Strangeness/merged_pythia_v2.5.root`
- closure references:
  - `output/KtoPi-MC-Gen-Closure.root`
  - `output/KtoPi-MC-Reco-Closure.root`

## Method

The theta-differential closure was built with the same main-chain ingredients as the nominal analysis:

- event selection:
  - `sum(RecoE)/91.2 > 0.5`
  - `Nch >= 7`
  - `30 < theta_thrust < 150 deg`
- track definition:
  - `RecoGoodTrack == 1`
  - charged tracks only
  - `0.4 <= pT < 5.0 GeV/c`
  - PID fiducial `0.15 < |cos(theta_track)| < 0.675`
- PID observation model:
  - exclusive `K / pi / p` categories
  - legacy tie rule `K > pi > p`

For `theta_track`, the same `(Nch_tag, pT)`-dependent correction matrix and reco/gen matching efficiencies used in the main chain were applied independently in each theta bin. This gives a faithful differential diagnostic without changing the production analysis outputs.

## Outputs

The new closure package is written to:

- `result/20260313/closure_pt_theta/Closure_vs_Pt_KPiP.pdf`
- `result/20260313/closure_pt_theta/Closure_vs_Theta_KPiP.pdf`
- `result/20260313/closure_pt_theta/Closure_vs_Pt_KPiP_Panels.pdf`
- `result/20260313/closure_pt_theta/Closure_vs_Theta_KPiP_Panels.pdf`
- `result/20260313/closure_pt_theta/Closure_vs_PtTheta_KPiP.root`
- `result/20260313/closure_pt_theta/Closure_PtTheta_Summary.md`

Main macro:

- `make_closure_vs_pt_theta.C`

## Quantitative summary

| Species | Axis | RMS(|closure-1|) | Worst bin center | Worst closure |
| --- | --- | ---: | ---: | ---: |
| K | pT | 0.0400 | 4.808 GeV/c | 1.080 |
| pi | pT | 0.0230 | 0.592 GeV/c | 1.055 |
| p | pT | 0.0834 | 0.975 GeV/c | 0.803 |
| K | theta | 0.0321 | 132.5 deg | 1.067 |
| pi | theta | 0.0306 | 57.5 deg | 1.053 |
| p | theta | 0.0763 | 57.5 deg | 0.888 |

## Detailed observations

### Kaons

`pT` closure:

- low `pT` starts slightly high: `1.030` at `0.59 GeV/c`
- mid `pT` is slightly low: about `0.964-0.980` from `1.0` to `2.5 GeV/c`
- high `pT` turns high again: `1.037-1.080` above `3.3 GeV/c`

`theta` closure:

- mostly flat at the percent level inside the accepted windows
- small edge structure is visible near the PID-fiducial boundaries:
  - `1.056` at `47.5 deg`
  - `1.041` at `77.5 deg`
  - `0.942` at `82.5 deg`
  - `0.966` at `97.5 deg`
  - `1.067` at `132.5 deg`

Interpretation:

- kaon closure is reasonably good overall
- the main theta dependence is localized near the edges of the accepted PID region rather than a broad angular trend

### Pions

`pT` closure:

- very stable overall
- worst bin is the first one: `1.055` at `0.59 GeV/c`
- all other bins stay close to unity, mostly within `1-3.5%`

`theta` closure:

- mild overclosure across most accepted bins
- typical values are `1.01-1.05`
- no strong localized failure at the fiducial edges

Interpretation:

- pion closure is the cleanest of the three species
- the theta trend is broad and mild, not dominated by a single geometric edge effect

### Protons

`pT` closure:

- this is the weakest species by a large margin
- worst bin is `0.803` at `0.975 GeV/c`
- most of the remaining `pT` bins are still low: roughly `0.90-0.95`
- the deficit is not limited to the tail; it is present through almost the full `pT` range

`theta` closure:

- also systematically low across both accepted theta lobes
- values are typically `0.89-0.96`
- worst bin is `0.888` at `57.5 deg`
- the deficit is visible across the accepted region, not only at the fiducial edges

Interpretation:

- the proton closure issue is not primarily a theta-edge effect
- it looks like a broader species-level normalization or response problem that survives both the `pT` and `theta` projections

## Main conclusion

The new differential checks sharpen the closure picture:

- `pi` closure is good in both `pT` and `theta`
- `K` closure is acceptable and the angular structure is mainly concentrated near the PID-fiducial boundaries
- `p` closure is the real outlier and remains low in both `pT` and `theta`

This strongly suggests that the remaining proton non-closure is not caused mainly by a geometric theta-acceptance mismatch. It is more consistent with a species-dependent issue in the proton correction chain itself, for example the proton tagging/mis-tagging matrix and/or the proton matching efficiencies.

## Recommended next checks

1. Compare proton closure with and without the PID fiducial to see whether the proton deficit is truly independent of the accepted theta windows.
2. Inspect the proton `(Nch_tag, pT)` cell matrices directly, especially the low-`pT` cells where the deficit is largest.
3. Compare the proton closure against the v2.4 nominal chain to determine whether the effect changed when switching to the v2.5 nominal MC.
