# X-SCAPE hadronization control for `e+e- -> gamma*/Z -> hadrons`

This note records the first direct control check on the local X-SCAPE
`e+e-` workflow.

The question was simple:
- is the low X-SCAPE `K/pi` result coming from the `epemGun` hard process or
  the ASCII-to-ROOT conversion,
- or is it specific to the current `hybrid` hadronization choice?

## Setup

The same `epemGun` hard-process path was used in both runs:
- `eCM = 91.2`
- `Tune:ee = 7`
- `PartonLevel:FSR = off` inside `epemGun`
- vacuum `Matter` path with
  - `in_vac = 1`
  - `Q0 = 1.0`
  - `QS = 1.0`
  - `vir_factor = 0.25`
  - `recoil_on = 0`
  - `broadening_on = 0`

The only intentional change between the two runs was the hadronization module:
- `hybrid`
- `colorless`

Shared hadronization-side settings:
- `pythia_decays = on`
- `tau0Max = 10.0`
- `reco_hadrons_in_pythia = 1`
- `eCMforHadronization = 0`
- `had_postprop = 0.0`
- `part_prop = 0.0`

Wrapper:
- [tools/run_xscape_generate_truth.sh](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/run_xscape_generate_truth.sh)

Comparison figure:
- [XscapeHadronization_Comparison.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_controls/XscapeHadronization_Comparison.pdf)

Comparison table:
- [XscapeHadronization_Comparison.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_controls/XscapeHadronization_Comparison.txt)

## Samples used

Hybrid run:
- [xscape_epem_hybrid_zpole_truth_clean20k.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_clean20k.root)

Colorless run:
- [xscape_epem_colorless_zpole_truth_20k.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_colorless/xscape_epem_colorless_zpole_truth_20k.root)

Standalone per-run `K/pi` summaries:
- [XscapeTruth_KtoPi_vs_dNdEta.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/XscapeTruth_KtoPi_vs_dNdEta.txt)
- [XscapeTruth_KtoPi_vs_dNdEta.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_colorless/XscapeTruth_KtoPi_vs_dNdEta.txt)

## Inclusive result

Hybrid:
- events: `20000`
- mean `nChEta05`: `7.3651`
- inclusive `K/pi` in `0.4 < pT < 5.0 GeV/c`: `0.090350`

Colorless:
- events: `20000`
- mean `nChEta05`: `6.9429`
- inclusive `K/pi` in `0.4 < pT < 5.0 GeV/c`: `0.156402`

For context, the current standalone generator set is:
- `PYTHIA8 baseline`: `0.146336`
- `HERWIG`: `0.151443`
- `SHERPA`: `0.157318`

So:
- X-SCAPE `colorless` is in-family with the existing standalone baseline band
- X-SCAPE `hybrid` is not

## Binned result

`hybrid / colorless` by `dN_ch/deta` bin:

- `1.41`: `0.6227`
- `5.22`: `0.5717`
- `9.03`: `0.6077`
- `12.84`: `0.6326`
- `16.66`: `0.5014`
- `20.47`: `0.4774`
- `24.28`: `0.3921`
- `28.09`: `0.4597`

So the hybrid setup is suppressed relative to colorless in every plotted bin,
roughly:
- `38%` low at the first point
- `40%` to `50%` low through the middle
- up to `60%` low in the sparse high-activity tail

## Interpretation

This control resolves the first engineering question.

The low X-SCAPE `hybrid` `K/pi` curve is **not** caused by:
- the `epemGun` hard-process path itself
- the basic final-state ASCII parsing
- the event-level counting logic for `nChEta05`, `nKPt0405`, and `nPiPt0405`

Reason:
- keeping the same hard process and conversion chain but switching only the
  hadronization module from `hybrid` to `colorless` moves the result from
  `K/pi ~ 0.09` to `K/pi ~ 0.156`

That means the present discrepancy is specific to the current X-SCAPE
`hybrid` hadronization setup.

## Practical consequence

At the current stage:
- the X-SCAPE `hybrid` curve can be shown as an **exploratory** model result
- but it should not yet be presented as a settled apples-to-apples generator
  baseline without qualification

The safer interpretation is:
- `colorless` validates the local `e+e-` X-SCAPE plumbing
- `hybrid` currently produces a much lower strange-to-nonstrange ratio and
  needs additional physics-side scrutiny before being merged into the main
  generator compilation without caveat

## Next checks

The next useful controls are:

1. Repeat the hybrid run with higher statistics.
2. Vary `pythia_decays` and `tau0Max` to see how much of the suppression comes
   from decay handling.
3. Check whether `reco_hadrons_in_pythia = 0` changes the result materially.
4. Inspect the X-SCAPE hybrid hadron species content directly before the final
   ROOT conversion.
