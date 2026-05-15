# 2026-03-10 iteration 4 unfolding-validation status

Implemented in this pass:
- refolding validation for `Ntag` unfolding
- refolding validation for `dNch/deta` unfolding
- mismatched-shape closure stress test for `Ntag` and `dNch/deta`
- migration purity/stability summary for K and pi responses

Quantitative summary:
- `Ntag` refolding RMS(refold/reco - 1): K MC = 0.0699, K data = 0.1456, Pi MC = 0.0319, Pi data = 0.2072
- `dNch/deta` refolding RMS(refold/reco - 1): K MC = 0.3157, K data = 0.3367, Pi MC = 0.3130, Pi data = 0.3357
- `Ntag` mismatched-shape RMS(unfolded/injected - 1) = 0.0636
- `dNch/deta` mismatched-shape RMS(unfolded/injected - 1) = 0.0474

Pending reviewer items:
- pseudo-experiment / bootstrap coverage study
- explicit SVD `kReg` stability scan and singular-value summary

Generated outputs:
- `result/20260310/unfolding_validation/NtagUnfolding_RefoldingValidation.pdf`
- `result/20260310/unfolding_validation/DNdEtaUnfolding_RefoldingValidation.pdf`
- `result/20260310/unfolding_validation/UnfoldingStressTest_Combined.pdf`
- `result/20260310/unfolding_validation/MigrationMetrics_PurityStability.pdf`
