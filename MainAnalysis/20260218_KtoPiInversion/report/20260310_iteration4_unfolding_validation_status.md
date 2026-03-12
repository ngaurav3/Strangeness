# 2026-03-10 iteration 4 unfolding-validation status

Superseded in part by:
- `report/20260312_iteration10_reviewer_followup.md`

The original iteration 4 package introduced the refolding, stress-test, and
migration-metric figures. The later 2026-03-12 follow-up completed the two
then-pending reviewer requests below (pseudo-experiment coverage and explicit
SVD `kReg` scan). The quantitative numbers listed here remain correct for the
iteration 4 figure set itself, but the current validation status should be read
together with the newer follow-up note.

Implemented in this pass:
- refolding validation for `Ntag` unfolding
- refolding validation for `dNch/deta` unfolding
- mismatched-shape closure stress test for `Ntag` and `dNch/deta`
- migration purity/stability summary for K and pi responses

Quantitative summary:
- `Ntag` refolding RMS(refold/reco - 1): K MC = 0.0435, K data = 0.1152, Pi MC = 0.0264, Pi data = 0.1725
- `dNch/deta` refolding RMS(refold/reco - 1): K MC = 0.4578, K data = 0.4560, Pi MC = 0.4549, Pi data = 0.4566
- `Ntag` mismatched-shape RMS(unfolded/injected - 1) = 0.0058
- `dNch/deta` mismatched-shape RMS(unfolded/injected - 1) = 0.0833

Originally pending reviewer items, completed later:
- pseudo-experiment / bootstrap coverage study
- explicit SVD `kReg` stability scan and singular-value summary

Generated outputs:
- `result/20260310/unfolding_validation/NtagUnfolding_RefoldingValidation.pdf`
- `result/20260310/unfolding_validation/DNdEtaUnfolding_RefoldingValidation.pdf`
- `result/20260310/unfolding_validation/UnfoldingStressTest_Combined.pdf`
- `result/20260310/unfolding_validation/MigrationMetrics_PurityStability.pdf`
