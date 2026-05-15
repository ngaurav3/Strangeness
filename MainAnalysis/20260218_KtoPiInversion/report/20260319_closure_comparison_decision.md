# 2026-03-19 closure comparison and method decision

## Scope
This note compares the closure and validation performance of the three relevant
`K/pi` vs `dN_ch/deta` correction branches that have been tested so far:

1. current nominal chain
2. approximate order-swapped cross-check
3. independent Yi-style `3D` cross-check in `mu = (dN_ch/deta, p, |cos(theta)|)`

The removed thrust-axis `dN_ch/dy` Yi-style appendix branch is intentionally not
part of this decision note.

## Definitions used for the comparison
### Nominal chain
1. fake / reco-match correction
2. reco-side `3x3` PID inversion
3. truth-side gen-matching correction
4. activity unfolding

### Approximate order-swapped branch
1. fake-correct observed reco tags
2. unfold observed tags in `dN_ch/deta`
3. apply truth-matched `3x3` inversion in `dN_ch/deta x pT`
4. apply truth-matched gen-efficiency correction

This remains a reduced factorized test, not the full Yi procedure.

### Independent Yi-style branch
1. fake-correct observed reco tags in flat `mu = (dN_ch/deta, p, |cos(theta)|)`
2. unfold observed tags in full flat `mu`
3. apply truth-cell `3x3` inversion in true `mu`
4. apply truth-cell matched/gen efficiency correction
5. integrate back to `K/pi` vs `dN_ch/deta`

## Core closure metrics
### Nominal branch
From `output/systematics_20260306_dndeta/nominal_unfold_dndeta.root`:
- closure RMS(`unfolded / truth - 1`): `0.022939`
- MC reco-space refolding RMS(`refolded / reco - 1`): `0.008957`
- data reco-space refolding RMS(`refolded / reco - 1`): `0.010087`

This is the current reference performance.

### Approximate order-swapped branch
From [20260319_order_swapped_dndeta_crosscheck.md](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/report/20260319_order_swapped_dndeta_crosscheck.md):
- closure RMS: `0.0691`
- MC reco-space refolding RMS: `0.0493`
- data reco-space refolding RMS: `0.0533`
- maximum `|shift| / nominal systematic`: `2.312`
- reordered closure by visible bin:
  - `0.986, 0.978, 0.994, 0.992, 1.024, 1.066, 1.104, 1.148`

Interpretation:
- clearly worse than nominal in both closure and refolding
- degradation is strongest in the merged tail bins
- the method shift exceeds the nominal quoted systematic from bin 6 onward

### Independent Yi-style `dN/deta` branch
From [20260319_yi_independent_dndeta_crosscheck.md](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/report/20260319_yi_independent_dndeta_crosscheck.md):
- closure RMS: `0.036203`
- maximum `|shift| / nominal systematic`: `1.758720`
- closure by visible bin:
  - `0.9957, 0.9796, 1.0063, 0.9849, 0.9876, 1.0014, 1.0073, 1.0978`
- first 7 visible bins agree with nominal within about `2%`
- last visible merged bin is about `12.4%` low relative to nominal

Important limitation:
- this branch currently has a strong closure test, but it does not yet have a
  full common-`mu` reco-space refolding package comparable to the nominal
  branch’s dedicated `refolded / reco` validation.
- So it is more independent than the approximate order-swapped test, but still
  less fully validated than the nominal chain.

## Direct comparison
### Closure quality ranking
Best to worst:
1. nominal: `0.0229`
2. independent Yi-style `dN/deta`: `0.0362`
3. approximate order-swapped: `0.0691`

### Refolding quality ranking
Only branches with a directly packaged reco-space refolding test can be ranked
strictly here:
1. nominal: `0.0090` (MC), `0.0101` (data)
2. approximate order-swapped: `0.0493` (MC), `0.0533` (data)

So the order-swapped branch is decisively worse than nominal on the most direct
unfolding-performance metric.

### Tail sensitivity
Both non-nominal branches agree reasonably well with the nominal result in the
populated low- and mid-activity bins.

The difference shows up in the final merged tail bin:
- approximate order-swapped: strongest instability and strongest self-nonclosure
- independent Yi-style: milder but still visible tail-only method sensitivity
- nominal: best overall closure and best refolding performance

That pattern matters. A method should not be promoted over the nominal result if
its largest disagreement appears exactly where its own validation is weakest.

## Informed decision
### Decision
Keep the current nominal `dN_ch/deta` branch as the quoted DELPHI result.

### Reasoning
1. It has the best overall closure.
2. It has by far the best reco-space refolding performance.
3. The approximate order-swapped branch is not competitive as a nominal method:
   it is worse on every validation metric that matters.
4. The independent Yi-style branch is technically meaningful and useful, but it
   does not outperform the nominal branch overall, and its only significant
   difference is confined to the merged tail bin.
5. The current significant method sensitivity is therefore a tail-bin issue, not
   evidence that the nominal low- and mid-activity result is unstable.

### Recommended treatment of the alternatives
- Approximate order-swapped branch:
  - keep only as an appendix stress test
  - do not use for quoted numbers
- Independent Yi-style `dN/deta` branch:
  - keep as a stronger appendix cross-check
  - use it to support the robustness of the nominal result in the first seven
    visible bins
  - do not replace the nominal chain unless a full common-`mu` refolding and
    systematics package is added and shown to be at least as stable as nominal

## Practical conclusion for the note
The current note should communicate:
- the nominal branch remains the quoted result
- alternative orderings have been tested
- low- and mid-activity conclusions are stable
- the merged highest-activity `dN_ch/deta` bin remains the method-sensitive
  region and should be interpreted with more caution than the rest of the curve
