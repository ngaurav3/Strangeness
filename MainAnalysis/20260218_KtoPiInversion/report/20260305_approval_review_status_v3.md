# DELPHI Strangeness Analysis Approval Review (v3)

**Author:** Yen-Jie Lee and OpenAI  
**Date:** 2026-03-05

---

## 1) Scope and Approval Goal

- Deliver a complete approval-style review of the DELPHI open-data `K/pi` analysis.
- Cover full chain: motivation, data, selections, correction, closure, unfolding, systematics, cross-checks, final results, and interpretation.
- Freeze current status as **v3** baseline.

---

## 2) Physics Motivation

- Study strangeness enhancement in the smallest clean system (`e+e-` at Z pole).
- Provide complementary information to pp, pA, AA multiplicity-dependent strangeness measurements.
- Test whether multiplicity-ordered hadrochemistry can arise without hadronic initial-state geometry.

![ALICE context](../../../report/assets/ALICE_strangeness_ratios_vs_mult_2017.pdf)

---

## 3) Dataset and Observable

- Collision system: `e+e-` at `sqrt(s) ~ 91.2 GeV` (DELPHI open data).
- Inputs:
  - `sample/Strangeness/merged_mc_v2.2.root`
  - `sample/Strangeness/merged_data_v2.2.root`
- Main observable:
  - `K/pi` vs `Ntag^ch`
  - final double ratio: `(K/pi)_Data / (K/pi)_MC`

---

## 4) Selection and PID Definition

- Event/track baseline:
  - `RecoGoodTrack == 1`
  - `MinNch = 7`
  - `30 deg < theta < 150 deg`
  - `0.2 < pT < 5.0 GeV/c`
- PID tagging:
  - `RecoPIDKaon >= 2`
  - `RecoPIDPion >= 2`
  - `RecoPIDProton >= 2`
- Exclusive category assigned by highest PID tag.

---

## 5) Correction Chain

For each `(Ntag, pT)` cell:

1. reco-matching correction (`RecoEfficiency*`)
2. `3x3` species unfolding (`RecoEfficiency*As*`)
3. gen-matching correction (`RecoGenEfficiency*`)
4. integrate pT and form `K/pi` vs `Ntag`

---

## 6) Produced Core Outputs

- `KtoPi-MC-Gen-Closure.root`
- `KtoPi-MC-Reco-Closure.root`
- `KtoPi-MC-Gen-Nominal.root`
- `KtoPi-MC-Reco-Nominal.root`
- `KtoPi-Data-Reco-Nominal.root`

---

## 7) Closure and Baseline Comparison

![closure/data-mc](../top_plots/KtoPiClosureVsData_Overlay.pdf)

- Mid-`Ntag` bins are stable.
- Residual non-closure slope remains in `Ntag`.
- Tail bins are statistically sparse.

---

## 8) Unfolding Upgrade in v3

- Added explicit `Ntag` unfolding for corrected K and pi yields.
- Bayes (`nIter=1`) used as nominal.
- SVD (`kReg=8`) kept as cross-check only.
- Response and method diagnostics produced and validated.

---

## 9) Unfolding Performance

![unfold closure](../top_plots/NtagUnfolding_MCClosure_BayesVsSVD.pdf)

![unfold datamc](../top_plots/NtagUnfolding_DataMC_BayesVsSVD.pdf)

- Bayes improves high-tail stability relative to no unfolding.
- SVD used for validation but not in uncertainty sum.

---

## 10) Systematics Definition (Current)

Components in final budget:

- acceptance envelope
- extrapolation
- binning
- closure residual (`|Closure_Bayes-1|`)
- unfolding uncertainty = `max(|prior variation|, |iteration variation|)`

Total: quadrature sum of all components.

---

## 11) Ntag Final Result + Systematics

![ntag final](../top_plots/KtoPi_DoubleRatio_with_Systematics.pdf)

- Plot convention updated:
  - **systematic-only** band
  - **statistical-only** vertical bars

---

## 12) dN/deta-based Result and Systematics

![dndeta components](../top_plots/Systematics_Components_vs_dNdEta.pdf)

![dndeta final](../top_plots/KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf)

- Same component logic propagated to dN/deta result.
- Final plot uses systematic-only band + statistical bars.

---

## 13) Cross-checks

- Bayes vs SVD method-difference QA retained in outputs.
- Central-eta `Ntag` test (`UseCentralEtaNtag=true`) for dN/deta chain performed.
- Result: closure did **not** improve in current setup; default `Ntag` retained.

---

## 14) Interpretation

- Data/MC double ratio is close to unity in controlled bins, with measurable structure.
- Residual closure/systematics indicate sensitivity to migration and acceptance in sparse regions.
- Current v3 supports robust status-level physics statements and identifies publication-grade improvements.

---

## 15) Conclusion and Approval Recommendation

- Full analysis chain is implemented, reproducible, and documented.
- Unfolding integrated in correction strategy and reflected in systematics.
- Final plots and note are updated to current v3 conventions.
- Recommendation:
  - approve v3 as analysis baseline,
  - proceed with publication-grade finalization (covariance, sparse-tail policy, response-by-variation).
