# 2026-03-10 iteration 5 parameter-choice summary

Status:
- stricter reviewer-style parameter-choice study
- not yet adopted as the production parameter-setting rule in the main note

Decision rule:
- choose the smallest Bayes iteration that passes refolding, mismatched-shape closure, and neighboring-step stability

- `Ntag` chosen iteration: no passing iteration in 1-6
- `dNch/deta` chosen iteration: no passing iteration in 1-6

Interpretation:
- under this stricter rule, the current study does not identify a fully passing Bayes iteration for either axis
- the production note still retains `nIter=1` as the current nominal setting from the older monotonic/minimal-deformation scan
- this stricter study should therefore be read as an open validation item, not as a superseding production choice

Outputs:
- `result/20260310/bayes_parameter_selection/Ntag_BayesParameterSelection_Summary.pdf`
- `result/20260310/bayes_parameter_selection/dNdEta_BayesParameterSelection_Summary.pdf`
- `result/20260310/bayes_parameter_selection/bayes_parameter_selection_summary.csv`
