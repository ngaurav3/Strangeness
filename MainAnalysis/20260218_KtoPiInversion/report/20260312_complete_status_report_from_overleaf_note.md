# Complete Status Report: DELPHI Strangeness Analysis

Author: Yen-Jie Lee and OpenAI  
Date: March 12, 2026

## Scope

This status report is synchronized to the latest Overleaf analysis note and summarizes:
- motivation and dataset
- observable definitions and selection
- correction and unfolding workflow
- validation status
- systematic uncertainty model
- final `K/pi` and `p/pi` results
- comparisons to ALICE and historical `e+e-` data

The PID scale-factor study is intentionally summarized only briefly here because Prof. Yi Chen's group will report that part in detail.

## Executive summary

- The main DELPHI result is the unfolded `K/pi` double ratio versus `Nch_tag` and versus `dNch/deta`.
- The nominal chain includes species correction, species-dependent unfolding, residual correction, and PID-SF propagation to data only.
- `K/pi` vs `Nch_tag` is the strongest and most mature result.
- `K/pi` vs `dNch/deta` now uses the same-variable reco/truth definition with `keepBins=8`; the ratio-level refolding is strong, but the merged species-level tail bin is still cautionary.
- The unfolding-check package now explicitly includes the keepBins scan, ratio-level refolding, pseudo-experiment coverage, migration/stress tests, and SVD `kReg` scan.
- `p/pi` is produced in parallel; the `dNch/deta` branch remains a status-level cross-check rather than a validation-quality headline result.

## Contents of the slide deck

1. Executive summary
2. Physics motivation
3. Data sample and observables
4. Baseline selection and PID definition
5. Analysis workflow
6. Pre-unfolding closure and data/MC status
7. `Nch_tag` unfolding setup and validation
8. `dNch/deta` unfolding setup and validation
9. `dNch/deta` ratio-level refolding and keepBins scan
10. Migration and stress-test diagnostics
11. Pseudo-experiment coverage and SVD `kReg` scan
12. PID scale factors: summary only
13. Quoted systematic uncertainty model
14. Main result: `K/pi` vs `Nch_tag`
15. Main result: `K/pi` vs `dNch/deta`
16. Comparison to ALICE and historical `e+e-` context
17. Inclusive `K/pi` cross-check versus world `e+e-`
18. Parallel `p/pi` results
19. Current strengths and caveats
20. Conclusion and next steps

## Key messages carried from the current note

### Main nominal result
- `K/pi` vs `Nch_tag` is the primary DELPHI multiplicity-differential output.
- Final points include:
  - species correction
  - Bayesian unfolding
  - residual non-closure correction
  - PID scale-factor update on data only

### `dNch/deta` branch
- The current branch uses a same-variable response with reco and truth both defined from charged multiplicity in `|eta| < 0.5`.
- The current working choice is `keepBins=8`, so the final visible point is a merged overflow bin.
- The active scripted default now matches that branch: `KEEPBINS_OVERRIDE = 8`.
- The note now distinguishes:
  - ratio-level `K/pi` refolding quality
  - species-level kaon/pion refolding quality
- The ratio-level refolding RMS is now about `0.0118` in MC and `0.0111` in data.
- The species-level final merged bin still remains around `0.33` for kaons and pions and should not be described as fully closed.

### `p/pi` branch
- `p/pi` vs `Nch_tag` is presented as a full parallel result.
- `p/pi` vs `dNch/deta` is retained, but the note explicitly states that it remains weaker than `K/pi`.

### Systematics
Quoted total includes four components:
- binning
- residual non-closure
- unfolding variation
- PID-SF propagation

Statistical bars shown in the current `K/pi` status plots are no longer the raw diagonal Bayes errors. They are the per-bin toy-calibrated absolute RMSE from the pseudo-experiment study because the direct diagonal propagation under-covers badly.
The DELPHI-vs-MC, DELPHI-vs-ALICE, and DELPHI-vs-world `dNch/deta` comparison figures now use that same toy-stat source.

Cross-checks retained but not included in the quoted total:
- acceptance variations
- extrapolation variation

## Files produced

- `report/20260312_complete_status_report_from_overleaf_note.tex`
- `report/20260312_complete_status_report_from_overleaf_note.pdf`
- `report/20260312_complete_status_report_from_overleaf_note.md`
- `report/20260312_unfolding_checks_status_note.md`

## Style

The slide deck uses the NTU-style Beamer theme derived earlier for the status-report series.
