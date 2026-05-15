# DELPHI Strangeness Analysis Approval Review (v3)

**Author:** Yen-Jie Lee and OpenAI  
**Date:** 2026-03-06

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
- In A+A, enhanced strange-hadron production is a classic QGP-sensitive observable.
- In high-multiplicity pp and p+A, heavy-ion-like strangeness trends are also observed.
- Key open question: are these trends driven by final-state hadronization alone, or by additional nuclear/medium effects?

![ALICE context](assets/ALICE_strangeness_ratios_vs_mult_2017.pdf)

---

## 3) Why `e+e-` Here

- `e+e-` provides a clean reference:
  - no hadronic remnant underlying event
  - well-defined hard scale
  - final state dominated by parton shower + hadronization
- This isolates hadronization-driven multiplicity effects from nuclear geometry/initial-state effects.
- LEP and B-factory identified-hadron measurements provide strong historical constraints on fragmentation models.

---

## 4) Physics Observable and Comparison Program

- Primary observable in this note:
  - `R_{K/pi}(Nch) = (K+ + K-) / (pi+ + pi-)` vs event activity
- In implementation we use `Ntag^ch` as activity estimator and report:
  - `(K/pi)_Data / (K/pi)_MC`
- Goal of final interpretation:
  - compare `e+e-` trend to pp, p+A, asymmetric A+B, and A+A results
  - identify which components scale with multiplicity and which require additional medium effects

---

## 5) Dataset and Observable

- Collision system: `e+e-` at `sqrt(s) ~ 91.2 GeV` (DELPHI open data).
- Inputs:
  - `sample/Strangeness/merged_mc_v2.2.root`
  - `sample/Strangeness/merged_data_v2.2.root`
- Main observable:
  - `K/pi` vs `Ntag^ch`
  - final double ratio: `(K/pi)_Data / (K/pi)_MC`

---

## 6) Selection and PID Definition

- Event/track baseline:
  - `RecoGoodTrack == 1`
  - `MinNch = 7`
  - `30 deg < theta < 150 deg`
  - `Ntag` from all charged tracks with `pT > 0.2 GeV/c`
  - identified-hadron spectra for `K/pi`: `0.4 < pT < 5.0 GeV/c`
- PID tagging:
  - `RecoPIDKaon >= 2`
  - `RecoPIDPion >= 2`
  - `RecoPIDProton >= 2`
- Exclusive category assigned by highest PID tag (for K/pi/p spectra).

---

## 7) Correction Chain

For each `(Ntag, pT)` cell:

1. reco-matching correction (`RecoEfficiency*`)
2. `3x3` species unfolding (`RecoEfficiency*As*`)
3. gen-matching correction (`RecoGenEfficiency*`)
4. integrate pT and form `K/pi` vs `Ntag`

---

## 8) Produced Core Outputs

- `KtoPi-MC-Gen-Closure.root`
- `KtoPi-MC-Reco-Closure.root`
- `KtoPi-MC-Gen-Nominal.root`
- `KtoPi-MC-Reco-Nominal.root`
- `KtoPi-Data-Reco-Nominal.root`

---

## 9) Closure and Baseline Comparison

![closure/data-mc](../top_plots/KtoPiClosureVsData_Overlay.pdf)

- Mid-`Ntag` bins are stable.
- Residual non-closure slope remains in `Ntag`.
- Tail bins are statistically sparse.

---

## 9b) Closure Test vs `eta`: Per-species Panels

![closure eta panels](../top_plots/Closure_vs_Eta_KPiP_Panels.pdf)

- Kaon (left), pion (middle), proton (right).
- Each panel shows generator/corrected yields (top) and closure ratio (bottom).

---

## 9c) Closure Test vs `eta`: Ratio Compilation

![closure eta overlay](../top_plots/Closure_vs_Eta_KPiP.pdf)

- Overlay closure ratios vs `eta` for `pi`, `K`, `p`.
- Y-axis range is fixed to `0-2`.

---

## 10) Unfolding Upgrade in v3

- Added explicit `Ntag` unfolding for corrected K and pi yields.
- Bayes (`nIter=1`) used as nominal.
- SVD (`kReg=8`) kept as cross-check only.
- Response and method diagnostics produced and validated.

---

## 11) Unfolding Performance

![unfold closure](../top_plots/NtagUnfolding_MCClosure_BayesVsSVD.pdf)

![unfold datamc](../top_plots/NtagUnfolding_DataMC_BayesVsSVD.pdf)

- Bayes improves high-tail stability relative to no unfolding.
- SVD used for validation but not in uncertainty sum.

---

## 12) Systematics Definition (Current)

Components in final budget:

- binning
- residual non-closure correction applied bin-by-bin from MC closure
- residual non-closure systematic = `50%` of correction magnitude
- unfolding uncertainty = `max(|prior variation|, |iteration variation|)`
- acceptance variation kept as cross-check only (excluded from final quadrature)

Total: quadrature sum of all components.

---

## 12a) Systematics uncertainties: Binning

- What is done: vary `pT` bin granularity (`NPtBins=10,14` around nominal 12).
- What we see:
- generally subdominant.
- median `~8e-4` (`Ntag`) and `~4e-4` (`dNch/deta`).

![binning ratio](../top_plots/Systematics_Binning_RatioToNominal.pdf)

---

## 12b) Systematics uncertainties: Residual Non-closure Correction

- What is done: apply bin-wise residual correction from MC closure.
- Residual systematic assigned as `50%` of correction magnitude.
- What we see:
- proxy scale (`0.5 x closure`) median `~0.013` (`Ntag`) and `~0.019` (`dNch/deta`).

![residual correction ratio](../top_plots/Systematics_ResidualCorrection_RatioToNominal.pdf)

---

## 12c) Systematics uncertainties: Unfolding

- What is done:
- nominal unfolding is Bayesian with `nIter=1` and nominal MC prior.
- prior variation: replace prior with flat prior (`hKPriorFlat`, `hPiPriorFlat`) at fixed `nIter=1`.
- iteration variation: keep nominal prior and vary to `nIter=2` (`nIterVar=max(2,nIter+1)`).
- build per-bin shifts on final double ratio:
  - `Delta_prior = DR(priorVar) - DR(nominal)`
  - `Delta_iter = DR(iterVar) - DR(nominal)`
- quoted unfolding uncertainty per bin is `max(|Delta_prior|, |Delta_iter|)`.
- What we see:
- `Ntag`: median `~0.010`, up to `~0.057`.
- `dNch/deta`: median `~0.019`, up to `~0.040`.

![unfolding ratio](../top_plots/Systematics_Unfolding_RatioToNominal.pdf)

---

## 12d) Unfolding Systematics Plot

![unfolding ratio large](../top_plots/Systematics_Unfolding_RatioToNominal.pdf)

---

## 13) Ntag Final Result + Systematics

![ntag final](../top_plots/KtoPi_DoubleRatio_with_Systematics.pdf)

- Plot convention updated:
  - **systematic-only** band
  - **statistical-only** vertical bars

---

## 14) dN/deta-based Result and Systematics

![dndeta components](../top_plots/Systematics_Components_vs_dNdEta.pdf)

![dndeta final](../top_plots/KtoPi_DoubleRatio_vs_dNdEta_with_Systematics.pdf)

- `dNch/deta` axis uses all stable charged particles in `|eta|<0.5` at generator level (no extra `pT` threshold).
- Same component logic propagated to dN/deta result.
- Final plot uses systematic-only band + statistical bars.

---

## 14b) Systematics Breakdown (Both Axes)

![sys dndeta](../top_plots/Systematics_Components_vs_dNdEta.pdf)

![sys ntag](../top_plots/Systematics_v3_Components_All.pdf)

- Explicit component breakdown is now shown for:
- `K/pi` vs `dNch/deta`
- `K/pi` vs `Ntag`
- Final components: binning, residual-correction, unfolding.
- Acceptance is shown as cross-check only.

---

## 15) New: `K/pi` vs `dNch/deta` (Data and MC)

![dndeta data mc](../top_plots/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf)

- Black points: unfolded DELPHI data (`K/pi`) with stat bars.
- Blue band: DELPHI data systematics propagated from accepted final components (binning, residual-correction, unfolding); acceptance shown as cross-check only.
- Red line: DELPHI MC gen-level `K/pi`.
- Open blue markers/band: unfolded MC cross-check with the same systematic propagation.

---

## 15b) ALICE vs DELPHI Compilation Plot

![delphi alice compilation](../top_plots/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf)

- DELPHI unfolded data + systematics and DELPHI MC gen-level curve.
- ALICE pp and Pb--Pb compilation points overlaid for external comparison.

---

## 15c) New: `p/pi` ALICE vs DELPHI Compilation

![ptopi delphi alice compilation](../top_plots/PtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf)

- DELPHI unfolded data + systematics and DELPHI MC gen-level curve.
- ALICE pp and Pb--Pb `p/pi` references overlaid for external comparison.

---

## 15d) Residual-Correction Validation

![kpi residual validation](../top_plots/ResidualCorrectionValidation_KtoPi_Ntag.pdf)

![ptopi residual validation](../top_plots/ResidualCorrectionValidation_PtoPi_Ntag.pdf)

- Left: `K/pi`, right: `p/pi`.
- Curves show: raw double ratio, closure, and corrected value (`raw/closure`).
- Residual uncertainty definition used in final budget: `0.5 * |1/closure - 1|`.

---

## 17) Cross-checks

- Bayes vs SVD method-difference QA retained in outputs.
- Central-eta `Ntag` test (`UseCentralEtaNtag=true`) for dN/deta chain performed.
- Result: closure did **not** improve in current setup; default `Ntag` retained.

---

## 18) Interpretation

- Data/MC double ratio is close to unity in controlled bins, with measurable structure.
- Residual closure/systematics indicate sensitivity to migration and acceptance in sparse regions.
- Current v3 supports robust status-level physics statements and identifies publication-grade improvements.

---

## 19) Conclusion and Approval Recommendation

- Full analysis chain is implemented, reproducible, and documented.
- Unfolding integrated in correction strategy and reflected in systematics.
- Final plots and note are updated to current v3 conventions.
- Recommendation:
  - approve v3 as analysis baseline,
  - proceed with publication-grade finalization (covariance, sparse-tail policy, response-by-variation).
