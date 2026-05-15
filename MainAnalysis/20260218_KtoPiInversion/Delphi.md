# Delphi Notes

This file is the canonical repo-local summary for the DELPHI `K/pi` analysis in
this workspace. It merges the previous `agent.md` and `DELPHI.md` content into
one place.

## 1. Scope

This repository is the main working area for the DELPHI strangeness analysis at
the `Z` pole. The central measurement is a fiducial charged `K/pi` ratio as a
function of event activity, with two main activity definitions:

- beam-axis `dN_ch/deta` in `|eta| < 0.5`
- thrust-axis `dN_ch/dy` in `|y_T| < 0.5`

The current note state uses an independent flat-`mu` correction workflow in the
full tuple:

- `mu = (activity, p, |cos(theta)|)`

rather than the older factorized workflow.

## 2. Git Repositories

This workspace involves multiple repositories and checkouts.

### Main code repository

- git root:
  - `/raid5/data/yjlee/strangeness/Strangeness`
- active analysis subdirectory:
  - `/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion`
- remote:
  - `origin = git@github.com:yenjie/Strangeness.git`

This is the main source-control home for:

- analysis code
- scripts
- tools
- standalone-generator code
- local reports in `report/`
- review comments and responses in `review/`

### Analysis-note repository

- local checkout:
  - `overleaf_repo_git/`
- remote:
  - Overleaf analysis-note project

This repo tracks:

- [overleaf_repo_git/main.tex](overleaf_repo_git/main.tex)
- note figures under `overleaf_repo_git/StrangenessV3/`
- note-side markdowns under `overleaf_repo_git/report/`

### Paper-draft repository

- local checkout:
  - `paper/overleaf_paper_draft/`
- remote:
  - Overleaf paper-draft project

This repo tracks:

- paper sections under `paper/overleaf_paper_draft/sections/`
- paper figures under `paper/overleaf_paper_draft/figures/`
- compiled paper PDF when intentionally tracked

### External generator-framework repositories

These are nested git checkouts under `external/`. They are not the main
DELPHI-analysis source-control homes, but they were used for exploratory
generator studies in this workspace.

#### JETSCAPE

- local checkout:
  - `external/JETSCAPE-src/`
- remote:
  - `origin = https://github.com/JETSCAPE/JETSCAPE.git`

Role in this workspace:

- upstream framework checkout used for exploratory JETSCAPE-based studies
- relevant to the `X-SCAPE/JETSCAPE` comparison work

#### X-SCAPE

- local checkout:
  - `external/X-SCAPE-src/`
- remote:
  - `origin = https://github.com/JETSCAPE/X-SCAPE.git`

Role in this workspace:

- upstream framework checkout used for the local X-SCAPE studies
- source of the exploratory `X-SCAPE/JETSCAPE` colorless and hybrid
  hadronization controls
- not the source-control home of the DELPHI analysis itself

### Upstream reduced-tree production repository

This is an important upstream reference even though it is not a local checkout
documented here as a working tree:

- reduced-tree production repo:
  - `https://github.com/FHead/Strangeness/tree/main/MainAnalysis/20260208_ReducedTreeProduction`

Use it when questions depend on how the reduced DELPHI trees were produced.

### Practical implication

Most substantial updates may require separate commits and pushes for:

1. main code repo
2. analysis-note repo
3. paper-draft repo

The external `JETSCAPE` and `X-SCAPE` clones should usually be treated as
upstream reference checkouts rather than repos to modify for routine DELPHI
analysis updates.

## 3. Repository Layout

Top-level areas worth knowing first:

- [KtoPiAnalysis.cpp](KtoPiAnalysis.cpp)
  - main event loop, reco/gen bookkeeping, core histograms, response inputs
- [include/TruthCountingPolicy.h](include/TruthCountingPolicy.h)
  - shared truth-level counting policy for charged activity and fiducial
    identified hadrons
- [analysis.sh](analysis.sh)
  - simple nominal driver for closure and nominal production
- [run_yi_independent_dndeta.py](run_yi_independent_dndeta.py)
  - current independent beam-axis correction chain
- [run_yi_independent_dndy.py](run_yi_independent_dndy.py)
  - current independent thrust-axis correction chain
- [finalize_dndeta_systematics.py](finalize_dndeta_systematics.py)
  - beam-axis result/systematics post-processing
- [finalize_dndy_systematics.py](finalize_dndy_systematics.py)
  - thrust-axis result/systematics post-processing
- `tools/`
  - standalone truth builders, generator converters, cross-check utilities
- `result/`
  - generated study outputs and control plots
- `output/`
  - main ROOT outputs from the DELPHI analysis chain
- [overleaf_repo_git/main.tex](overleaf_repo_git/main.tex)
  - analysis note source
- `paper/overleaf_paper_draft/`
  - paper draft source
- `report/`
  - status notes, study summaries, installation notes, internal markdowns
- `review/`
  - reviewer comments and responses
- [worklog.md](worklog.md)
  - prompt-oriented working history for this thread

## 4. Main Inputs

Active nominal input files:

- data: `sample/Strangeness/merged_data_v2.6.root`
- nominal MC: `sample/Strangeness/merged_pythia8ext_v2.6.root`

Important historical/cross-check inputs also exist, but the current note is
centered on the current `v2.6` data and `pythia8ext_v2.6` nominal MC samples.

## 5. Event Selection

The active event preselection uses the archived event bit:

- `PassAll == 1`

This replaced the earlier recomputed hadronic-event cuts. The analysis note and
code were aligned to that choice.

## 6. What the Quoted Observable Is

The quoted result is not a whole-event `4pi` inclusive `K/pi`.

It is a fiducial identified-hadron ratio:

- identified kaons and pions are restricted to:
  - `0.4 < p_T < 5.0 GeV/c`
  - `0.15 < |cos(theta_track)| < 0.675`

This is the standard DELPHI PID fiducial used in the analysis.

The event-activity variable is separate and event-level:

- `dN_ch/deta` branch:
  - charged reco tracks with `RecoGoodTrack==1` and `RecoCharge!=0`
  - counted in `|eta| < 0.5`
  - no PID requirement
  - no extra `p_T` threshold
- thrust-axis `dN_ch/dy` branch:
  - charged reco tracks with `RecoGoodTrack==1` and `RecoCharge!=0`
  - counted in `|y_T| < 0.5`
  - no PID requirement

So the quoted observable is:

- fiducial `K/pi` vs event-level activity

not:

- whole-event inclusive `K/pi`

## 7. Current Correction Workflow

The active workflow is:

1. Build observed reco tagged yields in flat reco `mu`
2. Apply reco fake/matching correction to the observed tags
3. Unfold observed tags from reco `mu` to true `mu`
4. Apply a truth-cell `3x3` PID/tag inversion
5. Apply truth-cell matched/gen efficiency correction
6. Project to the 1D activity axis
7. Form `K/pi`
8. Apply data-driven PID scale factor to data
9. Display the inherited total systematic envelope

This is the workflow described in the current note body.

## 8. Truth Counting Policy

The shared truth bookkeeping is centralized in
[include/TruthCountingPolicy.h](include/TruthCountingPolicy.h).

Important rules:

- charged activity uses a DELPHI-like charged-particle list, broader than just
  `pi/K/p`
- charged kaons and pions used in the ratio are counted only inside the fiducial
  `p_T` and angular window
- charged daughters of ancestors with `ctau > 1 cm` are vetoed in the
  DELPHI-like truth definition

That weak-daughter veto is meant to remove daughters of long-lived weak parents
such as:

- `K_S^0`
- `Lambda`
- `Xi`
- `Omega`

Important caveat:

- charm and bottom decays are not explicitly removed
- there is no dedicated `bbbar` / `ccbar` event veto
- the current result is therefore for inclusive hadronic `Z` decays, not a
  light-flavor-only sample

## 9. Heavy Flavor Caveat

The current DELPHI analysis does not remove kaons or pions from `D/B` decays.
The only ancestry veto in the truth policy is the long-lived weak-daughter veto
with `ctau > 1 cm`. Since `D` and `B` hadrons are much shorter lived than that
threshold, their decay products remain part of the measured fiducial `K/pi`.

This means:

- the quoted DELPHI `K/pi` result is inclusive in heavy flavor
- any heavy-flavor enhancement of kaon production is part of the current
  observable unless a dedicated tagged-veto/subtraction study is introduced

## 10. Activity Axes in Use

### Beam-axis branch

- observable: `K/pi` vs `dN_ch/deta`
- activity window: `|eta| < 0.5`
- final high-activity tail is merged into an explicit last visible bin

### Thrust-axis branch

- observable: `K/pi` vs `dN_ch/dy`
- `y_T = 0.5 * ln((E + p.n_thrust) / (E - p.n_thrust))`
- activity window: `|y_T| < 0.5`
- current branch uses rebinned variable activity edges:
  - `[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 15.5, 20.5, 25.5, 30.5]`

## 11. PID Scale Factors

The data-driven PID correction currently used for the `K/pi` result is:

- `SF_K = 0.9101` from `phi -> K+K-`
- `SF_pi = 0.9571` from `K* -> Kpi`

The active convention is:

- `R_after = R_before / (SF_K / SF_pi)`

This was previously implemented with the wrong sign convention and later fixed.

## 12. Standalone Generator Program

This repo also contains a large standalone generator comparison program.

Installed/built studies include:

- `PYTHIA 8.317`
- `PYTHIA 8.317 + Ropewalk`
- `PYTHIA 8.315 + Dire`
- `PYTHIA 6.428`
- `HERWIG 7.3.0`
- `SHERPA 3.0.3`
- exploratory `X-SCAPE/JETSCAPE` colorless and hybrid controls

These are used to build truth-level comparisons for:

- `K/pi` vs `dN_ch/deta`
- `K/pi` vs thrust-axis `dN_ch/dy`
- normalized activity-shape comparisons

The standalone curves were updated to follow the DELPHI-like truth definition,
including the weak-daughter veto.

## 13. Current Note / Paper State

Main document locations:

- note source: [overleaf_repo_git/main.tex](overleaf_repo_git/main.tex)
- paper source:
  - [paper/overleaf_paper_draft/sections/analysis.tex](paper/overleaf_paper_draft/sections/analysis.tex)
  - [paper/overleaf_paper_draft/sections/results_and_discussion.tex](paper/overleaf_paper_draft/sections/results_and_discussion.tex)

The analysis note is the most complete technical source of truth. The paper
draft mirrors the main physics story, but not every internal validation plot.

## 14. Internal Documentation Worth Reading

Good first internal notes:

- [report/20260319_status_summary_so_far.md](report/20260319_status_summary_so_far.md)
- [report/20260314_dndy_thrust_branch.md](report/20260314_dndy_thrust_branch.md)
- [report/20260314_longlived_daughter_veto_study.md](report/20260314_longlived_daughter_veto_study.md)
- [report/20260314_standalone_generator_veto_followup.md](report/20260314_standalone_generator_veto_followup.md)
- [report/20260319_yi_independent_dndeta_crosscheck.md](report/20260319_yi_independent_dndeta_crosscheck.md)
- [report/20260319_closure_comparison_decision.md](report/20260319_closure_comparison_decision.md)

## 15. Problem-Solving Style

Use this style when working in this DELPHI analysis repository.

### 15.1 Start from the code, not from assumptions

- Verify behavior in the active implementation before answering.
- Prefer code references over memory.
- If note text, paper text, and code disagree, treat the code as the first thing
  to audit and then align the documents explicitly.

### 15.2 Keep nominal, legacy, and cross-check workflows separate

- Be explicit about which workflow is currently quoted.
- Do not silently mix legacy and current results.
- When presenting an alternative method, label it clearly as:
  - nominal
  - legacy
  - cross-check
  - exploratory

### 15.3 Preserve DELPHI-specific definitions

- Do not casually reinterpret the observable.
- Keep the distinction explicit between:
  - fiducial identified `K/pi`
  - event-level activity estimator
  - whole-event inclusive quantities
- State the active fiducial window when relevant:
  - `0.4 < p_T < 5.0 GeV/c`
  - `0.15 < |cos(theta_track)| < 0.675`

### 15.4 Prefer validation over aesthetics

- Before changing the quoted workflow, compare:
  - closure
  - refolding
  - method dependence
  - systematic impact
- If a cross-check looks visually better but performs worse globally, document
  the tradeoff and keep the more defensible nominal choice.

### 15.5 Keep repository scopes separate

Work on the main code repo, note repo, paper repo, and upstream/external repos
deliberately:

- commit/push each repo separately
- do not assume note/paper changes imply code changes, or vice versa
- avoid mixing unrelated local dirt into commits

### 15.6 Document nontrivial changes immediately

- Add or update markdown in `report/` for meaningful analysis changes.
- Mirror note-relevant markdown to `overleaf_repo_git/report/` when useful.
- If a change affects interpretation, note the caveat explicitly.

### 15.7 Be careful with generated assets

- When figures are updated, confirm the correct files are rebuilt.
- Check that the note or paper actually references the refreshed assets.
- If a user asks to push, verify whether there are real tracked changes first.

### 15.8 Keep the scientific caveats visible

Important recurring caveats in this repo include:

- weak-decay-daughter veto versus inclusive stable-particle counting
- heavy-flavor feed-down not being explicitly removed
- inherited versus fully regenerated systematic envelopes
- exploratory status of some generator curves

Do not hide these caveats when they materially affect interpretation.

### 15.9 Use concise, direct communication

- Answer the exact question first.
- Then give the operational or physics caveat if needed.
- Keep wording concrete and auditable.

### 15.10 When in doubt, leave a trace

- If a workflow choice, cross-check result, or method comparison required real
  reasoning, capture it in markdown so the next pass does not need to reconstruct
  it from memory.

## 16. Work Log Rule

- Record all visible user prompts and the corresponding work in
  [worklog.md](worklog.md).
- Treat `worklog.md` as a prompt-oriented project history, not as a shell log.

## 17. Practical Starting Points

If you are new to the repo, the most efficient entry path is:

1. Read [overleaf_repo_git/main.tex](overleaf_repo_git/main.tex)
2. Read [KtoPiAnalysis.cpp](KtoPiAnalysis.cpp)
3. Read [include/TruthCountingPolicy.h](include/TruthCountingPolicy.h)
4. Read [analysis.sh](analysis.sh)
5. Read the beam-axis and thrust-axis drivers:
   - [run_yi_independent_dndeta.py](run_yi_independent_dndeta.py)
   - [run_yi_independent_dndy.py](run_yi_independent_dndy.py)
6. Then use the `report/` markdowns for study history and caveats

## 18. Current Physics Interpretation

The working interpretation in the note is:

- DELPHI does not show a strong multiplicity-driven rise in this fiducial
  `K/pi` observable
- the thrust-axis `dN_ch/dy` branch is a useful `e+e-`-native cross-check
- standalone generator comparisons are now more defensible because they use a
  DELPHI-like counting convention

Important caveats still active:

- heavy-flavor feed-down is included
- the current systematic band for the present independent workflow is still
  inherited from the validated older envelope in some places, rather than from a
  fully fresh standalone variation campaign
- `X-SCAPE/JETSCAPE` hybrid curves remain exploratory and hadronization-model
  sensitive

## 19. Maintenance Note

This file is meant to be a living summary. If the active workflow, fiducial
definition, truth policy, repository layout, or note/paper narrative changes,
update this file in the same pass.
