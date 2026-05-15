# 2026-03-14 Long-lived weak-decay daughter veto study

## Question

The DELPHI measurement excludes charged particles that come from decays with
`c tau > 1 cm`, e.g. daughters of `K_S^0`, `Lambda`, `Xi`, and similar weakly
decaying hadrons. The requested study was:

- compare `K/pi` vs `dN_ch/deta` with and without those daughters
- apply the comparison separately to:
  - the `K/pi` numerator/denominator
  - the `dN_ch/deta` activity definition

The four variants used were:

1. include daughters in both `K/pi` and `dN/deta`
2. veto daughters in `K/pi` only
3. veto daughters in `dN/deta` only
4. veto daughters in both

## DELPHI merged truth result

Study inputs:

- MC file: `sample/Strangeness/merged_pythia_v2.5.root`
- event selection: `PassAll == 1`
- `K/pi` definition:
  - charged `K^\pm` and `pi^\pm`
  - `0.4 < p_T < 5.0 GeV/c`
  - same PID-fiducial momentum cut used in the nominal truth chain:
    - `0.15 < |cos(theta)| < 0.675`
- `dN_ch/deta` definition:
  - charged truth particles in `|eta| < 0.5`
  - no additional activity `p_T` cut

Output files:

- `result/20260314/longlived_daughter_veto/longlived_daughter_veto_dndeta.root`
- `result/20260314/longlived_daughter_veto/longlived_daughter_veto_dndeta.txt`
- `result/20260314/longlived_daughter_veto/LongLivedDaughterVeto_KtoPi_vs_dNdEta.pdf`
- `result/20260314/longlived_daughter_veto/LongLivedDaughterVeto_dNdEta_Distribution.pdf`

Observed result:

- processed events: `1,334,110`
- selected `PassAll` events: `1,101,081`
- total charged particles in `|eta|<0.5`:
  - include: `9,514,225`
  - veto: `9,514,225`
- total charged kaons in the `K/pi` definition:
  - include: `962,712`
  - veto: `962,712`
- total charged pions in the `K/pi` definition:
  - include: `6,465,205`
  - veto: `6,465,205`

Inclusive `K/pi`:

- include/include: `0.148907`
- veto in `K/pi` only: `0.148907`
- veto in `dN/deta` only: `0.148907`
- veto in both: `0.148907`

All eight `dN/deta` bins are identical within exact floating-point equality in
the current merged-truth study.

## Why the merged-truth result is zero

Two explicit checks were performed.

### 1. Full selected-event ancestor scan

Using the stored `GenParent` chain on the full selected sample:

- selected `PassAll` events: `1,101,081`
- charged status-1 `Gen` particles scanned: `20,379,344`
- charged status-1 particles with any stored
  `K_S^0 / K_L^0 / Lambda / Xi / Omega` ancestor: `0`

This shows that the active `Gen` record used by the analysis does **not**
contain charged final particles with those weak ancestors in the stored
truth ancestry.

### 2. Status-1 charged species content of the selected `Gen` record

In the selected `PassAll` truth record, the charged hadron content relevant for
this analysis is effectively:

- `pi^\pm`: `16,394,974`
- `K^\pm`: `2,538,364`
- `p/\bar p`: `892,571`

No `K_S^0`, `Lambda`, `Xi`, or `Omega` appear as status-1 objects in the truth
list used by the main chain.

## Interpretation of the DELPHI tree bookkeeping

The main conclusion is:

- the current DELPHI merged `Gen` record already behaves like the
  **post-veto** definition
- i.e. the truth record used in the unfolding support does not expose the
  long-lived weak-decay daughters as counted charged final particles

There are explicit `KShort*` branches in the tree, but those daughter links do
not point into the `Gen` arrays. They point into the separate `Sim` record.
That means:

- `Gen` is the relevant truth definition for the current nominal analysis
- `Sim` is not a drop-in replacement for a direct “with daughters” variation
  because its particle coding is different and is not used by the main chain

## Standalone PYTHIA8 control

To quantify the size of the effect in a fully inclusive stable-particle event
record, a standalone control was run with:

- generator: `PYTHIA 8.317`
- tune: `Tune:ee = 7` (`Monash2013-ee`)
- process: `e+e- -> gamma*/Z -> hadrons` at `sqrt(s)=91.2 GeV`
- generated events: `400,000`

Output files:

- `result/20260314/pythia8_truth/pythia8_longlived_daughter_veto_400k.root`
- `result/20260314/pythia8_truth/pythia8_longlived_daughter_veto_400k.txt`
- `result/20260314/pythia8_truth/Pythia8_LongLivedDaughterVeto_KtoPi_vs_dNdEta.pdf`
- `result/20260314/pythia8_truth/Pythia8_LongLivedDaughterVeto_dNdEta_Distribution.pdf`

Observed result:

- total charged particles in `|eta|<0.5`:
  - include: `3,365,181`
  - veto: `2,993,721`
  - removed fraction: `11.0%`
- charged kaons in the `K/pi` definition:
  - include: `310,350`
  - veto: `310,268`
  - removed fraction: `0.026%`
- charged pions in the `K/pi` definition:
  - include: `2,297,447`
  - veto: `2,045,993`
  - removed fraction: `10.9%`

Inclusive `K/pi`:

- inclusive stable definition: `0.135085`
- weak-daughter-vetoed definition: `0.151647`

Largest removed parent contributions:

- for `dN/deta`:
  - `K_S^0`: `258,898`
  - `Lambda^0`: `75,035`
  - `Sigma^+`: `17,067`
  - `Sigma^-`: `15,462`
  - `Xi^-`: `4,705`
- for removed pions:
  - overwhelmingly `K_S^0`, then `Lambda` and hyperons
- for removed kaons:
  - negligible, dominated by a tiny `Omega^-` contribution

## Take-home message

The study separates two distinct facts.

1. **Current DELPHI merged truth**
   - the active `Gen` record already corresponds to the long-lived-daughter
     vetoed definition
   - therefore the current DELPHI `K/pi` vs `dN/deta` result is unchanged by
     applying the veto again

2. **Inclusive stable-particle control**
   - if one starts from a fully inclusive stable-particle generator record,
     the veto matters a lot
   - it removes about `11%` of the charged midrapidity activity and about
     `11%` of the counted charged pions in the `K/pi` definition
   - because kaons are almost unaffected while pions decrease strongly,
     the vetoed `K/pi` ratio is significantly larger than the inclusive one

So the correct interpretation is:

- the DELPHI analysis is already using the intended weak-decay-daughter-vetoed
  truth definition
- the reason the current tree-level DELPHI study gives exactly zero difference
  is not that the veto is unimportant in principle
- it is that the veto has already been built into the truth bookkeeping of the
  merged `Gen` record used by the analysis
