# 2026-03-12 unfolding checks status note

Author: Yen-Jie Lee and OpenAI

## Scope
This note records the additional unfolding checks propagated into the updated status-report slides on March 12, 2026.

## What was added to the status deck
- explicit distinction between species-level and ratio-level `dN_{ch}/d\eta` refolding
- explicit `keepBins = 8, 9, 10` scan summary
- pseudo-experiment coverage study summary
- explicit SVD `kReg` scan summary
- updated interpretation of the current `dN_{ch}/d\eta` branch

## Current working `dN_{ch}/d\eta` branch
- reco and truth axes are both charged multiplicity in `|eta| < 0.5`
- active working choice: `keepBins = 8`
- last visible point is a merged overflow bin
- active scripted default: `run_systematics_dndeta_unfolding.py` now defaults to `KEEPBINS_OVERRIDE = 8`
- authoritative summary builder: `finalize_dndeta_systematics.py`

## Key validation numbers used in the slides
### Ratio-level refolding for `K/pi` vs `dN_{ch}/d\eta`
- MC RMS: `0.0118`
- data RMS: `0.0111`

### keepBins scan
- best tested point: `keepBins = 8`
- none of `keepBins = 8, 9, 10` fully passes the stricter species-level rule
- final merged-bin species-level refolding remains about `0.33` for kaons and pions at `keepBins = 8`

### Pseudo-experiment coverage
- raw diagonal Bayes statistical errors under-cover strongly
- current status plots therefore use toy-calibrated per-bin absolute RMSE as the statistical bar
- the `K/pi` `dN_{ch}/d\eta` comparison/context figures now use that same toy-stat source

### SVD cross-check
- SVD remains strongly `kReg`-dependent
- it is kept as a qualitative cross-check, not a quoted uncertainty input

## Presentation consequence
The updated status deck now presents the `dN_{ch}/d\eta` branch honestly:
- acceptable at ratio level for the current `K/pi` observable
- still cautionary at species level in the merged tail bin
- stronger than before, but not promoted to the same validation status as the nominal `N_{ch}^{tag}` result

## Files updated
- `report/20260312_complete_status_report_from_overleaf_note.tex`
- `report/20260312_complete_status_report_from_overleaf_note.pdf`
- `report/20260312_complete_status_report_from_overleaf_note.md`
- `report/20260312_unfolding_checks_status_note.md`
