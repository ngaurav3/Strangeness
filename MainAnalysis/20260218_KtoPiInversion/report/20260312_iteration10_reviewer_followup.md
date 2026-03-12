# 20260312 reviewer follow-up validation

Author: Yen-Jie Lee and OpenAI

## Scope
This note records the remaining unfolding-validation studies requested by the reviewer after iteration 9:
- explicit `keepBins` scan for the `dN_{ch}/d\eta` tail merge
- pseudo-experiment coverage study
- explicit SVD `kReg` scan with singular-value summary

## keepBins scan
Scanned values:
- `keepBins = 8, 9, 10`

Acceptance rule used in the scan:
1. ratio-level `K/pi` refolding remains acceptable
2. species-level final-bin refolding no longer collapses catastrophically
3. final merged-bin purity and stability improve beyond the marginal `0.20-0.25` level

- Best ratio-level refolding stability in this scan: `keepBins=8` with data RMS `0.0111` and MC RMS `0.0118`.
- Passing scan point(s): none under the current explicit rule.
- Best scanned final-bin species refolding remains limited: K data `0.327`, pi data `0.336`.
- Best scanned final-bin migration quality remains: K purity/stability `0.454/0.489`, pi purity/stability `0.466/0.480`.

Generated outputs:
- `result/20260312/reviewer_followup_validation/keepbins_scan/DNdEtaKeepBinsScan_Summary.pdf`
- `result/20260312/reviewer_followup_validation/keepbins_scan/DNdEtaKeepBinsScan_SpeciesClosure.pdf`

## Pseudo-experiment coverage
- `Ntag` toy coverage: mean |pull bias| `1.979`, mean pull width `0.348`, minimum coverage `|pull|<1` `0.000`, `|pull|<2` `0.000`.
- `dNdEta` toy coverage: mean |pull bias| `2.019`, mean pull width `0.518`, minimum coverage `|pull|<1` `0.000`, `|pull|<2` `0.000`.

Generated output:
- `result/20260312/reviewer_followup_validation/toy_coverage/UnfoldingToyCoverage_Summary.pdf`

## SVD regularization scan
- `Ntag` SVD scan: smallest closure RMS at `kReg=2` (0.4431); smallest Bayes-SVD method-difference RMS at `kReg=2` (0.0379); production cross-check value remains `kReg=8`.
- `dNdEta` SVD scan: smallest closure RMS at `kReg=2` (0.0341); smallest Bayes-SVD method-difference RMS at `kReg=2` (0.0302); production cross-check value remains `kReg=8`.

Generated output:
- `result/20260312/reviewer_followup_validation/svd_kreg_scan/SVD_kReg_Scan_Summary.pdf`

## Bottom line
- The note-side reviewer comments from iteration 8 and 9 are addressed.
- The `dN_{ch}/d\eta` tail problem is now explicitly scanned rather than only discussed, but the species-level final-bin issue is not fully closed by the present `keepBins` scan alone.
- The pseudo-experiment coverage and explicit SVD `kReg` scan now exist as concrete validation studies rather than open placeholders.
