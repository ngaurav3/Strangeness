# 2026-03-12 v5 duplicate-PID-candidate status

## Scope

Requested v5 change:

- replace the exclusive K > pi > p tie-breaking observation with an inclusive raw-spectrum definition
- fill every kaon candidate with `RecoPIDKaon >= 2` into the raw kaon spectra
- fill every pion candidate with `RecoPIDPion >= 2` into the raw pion spectra
- fill every proton candidate with `RecoPIDProton >= 2` into the raw proton spectra
- allow duplicate counting across K / pi / p for the same track

Implemented in code:

- `KtoPiAnalysis.cpp`
  - new `PIDObservationMode` / `AllowDuplicatePIDCandidates` option
  - `exclusive` remains the default
  - `inclusive` fills all species with PID score `>= 2`
- `run_v5_duplicate_pid.sh`
  - dedicated v5 driver with isolated output/result directories

## Validation outcome

The current DELPHI trees already behave as effectively exclusive inputs at the threshold used here.

Full-sample tree-level check:

- MC (`sample/Strangeness/merged_mc_v2.2.root`)
  - `max duplicate candidates per event = 0`
- data (`sample/Strangeness/merged_data_v2.2.root`)
  - `max duplicate candidates per event = 0`

Definition used for the duplicate check:

- reconstructed good charged tracks with
  - `RecoGoodTrack == 1`
  - `RecoCharge != 0`
- duplicate candidate means at least two of
  - `RecoPIDKaon >= 2`
  - `RecoPIDPion >= 2`
  - `RecoPIDProton >= 2`

Observed result:

- no event in the current MC or data sample contains any track with multiple species passing the `>= 2` PID threshold

## Smoke-test comparison

Bounded smoke test on 20k MC events:

- old chain: exclusive observation
- new chain: inclusive duplicate-allowed observation

Raw histogram comparison:

- `hK`: identical bin-by-bin
- `hPi`: identical bin-by-bin
- `hP`: identical bin-by-bin
- `hKPt`, `hPiPt`, `hPPt`: identical bin-by-bin
- `hKDNdEta`, `hPiDNdEta`, `hPDNdEta`: identical bin-by-bin

Maximum absolute bin difference in the tested raw spectra:

- `0.0` for every checked histogram

## Consequence

For the current DELPHI PID branches, the requested v5 duplicate-candidate observation is a strict no-op:

- raw observed spectra do not change
- corrected spectra do not change
- unfolding inputs do not change
- final `K/pi` and `p/pi` results do not change

Because the analysis is exactly unchanged on the current trees, a full v5 rerun was not executed. It would only regenerate the same plots and tables.

## Practical status

What is available now:

- code path for inclusive duplicate-allowed PID observation
- v5 driver script prepared for isolated reruns:
  - `run_v5_duplicate_pid.sh`

What was intentionally not done:

- no new v5 plot set was produced
- no analysis-note or Overleaf figures were changed

Reason:

- the input PID tags already contain no multi-species `>= 2` candidates, so the v5 observable definition is numerically identical to the current chain

## Recommendation

If we want a materially different v5 PID observation model, it needs to be defined from a different input than the current integer `RecoPIDKaon/Pion/Proton >= 2` thresholds, for example:

- using non-exclusive per-species probabilities directly
- lowering the threshold and accepting overlaps there
- redefining observed spectra from calibration weights rather than hard threshold tags
