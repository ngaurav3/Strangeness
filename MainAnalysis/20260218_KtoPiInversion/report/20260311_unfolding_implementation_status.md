# 2026-03-11 Unfolding Implementation Status

## Summary

The active analysis in `MainAnalysis/20260218_KtoPiInversion` does **not** use RooUnfold.
The unfolding code is implemented locally in the analysis macros using ROOT matrix tools.

## Active implementation

### Nch_tag unfolding
- Macro: `runNtagUnfolding_BayesSVD.C`
- Bayesian unfolding is implemented in the local function `IterativeBayesUnfold1D`.
- SVD unfolding is implemented in the local function `SVDUnfold1D`.

Relevant code locations:
- `runNtagUnfolding_BayesSVD.C:84` for `IterativeBayesUnfold1D`
- `runNtagUnfolding_BayesSVD.C:165` for `SVDUnfold1D`

### dNch/deta unfolding
- Macro: `runDNdEtaUnfolding_BayesSVD.C`
- Bayesian unfolding is implemented in the local function `IterativeBayesUnfold1D`.
- SVD unfolding is implemented in the local function `SVDUnfold1D`.

Relevant code locations:
- `runDNdEtaUnfolding_BayesSVD.C:219` for `IterativeBayesUnfold1D`
- `runDNdEtaUnfolding_BayesSVD.C:298` for `SVDUnfold1D`

## What is used instead of RooUnfold

The current code uses ROOT linear algebra classes directly:
- `TMatrixD`
- `TVectorD`
- `TDecompSVD`

The macros do **not** include or instantiate:
- `RooUnfoldResponse`
- `RooUnfoldBayes`
- `RooUnfoldSvd`
- any `RooUnfold*.h` header

## ROOT environment

The active ROOT version in this workspace is:
- `6.34.04`

## Practical implication

When the note or slides refer to “Bayesian unfolding” or “SVD unfolding”, this currently means:
- a custom iterative Bayes implementation in the macro
- a custom SVD inversion based on `TDecompSVD`

It does **not** mean a RooUnfold package dependency.

## Recommendation for documentation

The analysis note and slides should avoid wording that implies RooUnfold was used unless the code is actually migrated to RooUnfold in the future.
A precise wording is:

> The unfolding is implemented in local ROOT macros using an iterative Bayes method and an SVD-based matrix inversion.

## Status

- RooUnfold dependency in active chain: **no**
- Custom Bayes/SVD implementation in active chain: **yes**
- ROOT version: **6.34.04**
