# 2026-03-14 standalone generator weak-daughter veto follow-up

## Starting point

This note is the follow-up to:

- `report/20260314_longlived_daughter_veto_study.md`
- `report/20260314_dndy_thrust_branch.md`

The earlier veto study established two things:

1. the active DELPHI merged `Gen` record already behaves like the intended
   weak-decay-daughter-vetoed truth definition
2. a fully inclusive standalone stable-particle record changes significantly if
   daughters of `c tau > 1 cm` ancestors are removed

The requested next step was therefore:

- make the standalone generator studies use the DELPHI-like veto definition by
  default
- apply it consistently to:
  - `K/pi`
  - `dN_ch/deta`
  - thrust-axis `dN_ch/dy`
- rebuild the comparison plots
- update the analysis note and paper draft

## Scope completed

The standalone generator comparison chain is now aligned to the DELPHI-like
charged-particle definition for:

- `PYTHIA 8.317`
- `PYTHIA 8.317 + Ropewalk`
- `PYTHIA 8.315 + Dire`
- `PYTHIA 6.428`
- `HERWIG 7.3.0`
- `SHERPA 3.0.3`
- `X-SCAPE/JETSCAPE` control curves

For all of the standalone generators above, the default event-level counting
used in the comparison plots now excludes charged daughters of ancestors with
`c tau > 1 cm`.

## Observable definition after the update

### `K/pi`

- charged `K^\pm / pi^\pm`
- `0.4 < p_T < 5.0 GeV/c`
- weak-decay daughters from ancestors with `c tau > 1 cm` removed

### `dN_ch/deta`

- charged final particles in `|eta| < 0.5`
- no additional activity `p_T` cut
- weak-decay daughters from ancestors with `c tau > 1 cm` removed

### thrust-axis `dN_ch/dy`

- charged final particles in `|y_T| < 0.5`
- `y_T = 0.5 ln[(E + p \cdot n_{thrust}) / (E - p \cdot n_{thrust})]`
- no additional activity `p_T` cut
- weak-decay daughters from ancestors with `c tau > 1 cm` removed

Important implementation detail:

- for `dN/dy`, the veto is applied to the particles counted in `dN/dy` and
  `K/pi`
- the thrust axis is still reconstructed from the full stored final-state
  sample, not from the vetoed subset

## Source changes made

### Standalone truth builders

Updated:

- `tools/pythia8_generate_truth_kpi_dndeta.cc`
- `tools/hepmc3_to_truth_root.cc`
- `tools/pythia6_generate_truth_kpi_dndeta.f90`
- `tools/pythia6_final_state_to_truth_root.py`
- `tools/xscape_final_state_to_truth_root.py`
- `tools/build_truth_kpi_vs_dndy_from_tree.cc`

New / replacement utility:

- `tools/herwig_ascii_to_truth_root.py`

Wrapper fixes:

- `tools/run_herwig_generate_truth.sh`
- `tools/run_truth_dndy_from_tree.sh`

### What changed in the standalone truth outputs

The event-level branches in the standalone truth trees now mean the vetoed
definition by default:

- `nChEta05`
- `nPiPt0405`
- `nKPt0405`
- `kPiPt0405`

Inclusive bookkeeping is still retained in parallel branches:

- `nChargedInclusive`
- `nChEta05Inclusive`
- `nPiPt0405Inclusive`
- `nKPt0405Inclusive`
- `kPiPt0405Inclusive`

Per-particle bookkeeping now includes:

- `isWeakDecayDaughter`

For the thrust-axis outputs, a metadata string is also written:

- `TNamed("weakDecayDaughterPolicy", "...")`

## Generator-specific notes

### PYTHIA8, SHERPA, and HepMC-based paths

The veto is applied explicitly by traversing stored ancestry and flagging
particles whose ancestry contains a long-lived weak parent.

### PYTHIA6

The standalone chain was upgraded so that the saved event content is rich
enough to support the same DELPHI-like veto logic downstream.

### HERWIG

The original HepMC3 reader path was not reliable on the local Herwig Ascii v3
output. The observed failure mode was:

- `ERROR::ReaderAscii: too few vertices were parsed`

This was replaced with a dedicated Herwig Ascii parser:

- `tools/herwig_ascii_to_truth_root.py`

This parser reconstructs enough production-vertex ancestry to apply the same
`c tau > 1 cm` ancestor veto used in the other standalone converters.

### X-SCAPE/JETSCAPE

This case is slightly different. The current X-SCAPE/JETSCAPE study starts from
`final_state_hadrons.dat`, which already contains hadrons like `K_S^0`,
`K_L^0`, `Lambda`, `Xi`, etc. directly rather than their charged final decay
daughters.

So for X-SCAPE/JETSCAPE:

- the DELPHI-like weak-daughter-veto behavior is effectively implicit
- the hadron-level input is already closer to the vetoed definition
- `isWeakDecayDaughter` stays `0` in the current local converter because the
  charged daughters are not present as explicit final-state records to be
  flagged

## Regeneration done

The updated standalone truth roots were regenerated for the current high-stat
samples:

- `result/20260314/pythia8_truth/pythia8_zpole_truth_400k.root`
- `result/20260314/pythia8_rope_truth/pythia8_zpole_truth_rope_400k.root`
- `result/20260314/pythia8_dire_truth/pythia8_zpole_truth_dire_400k.root`
- `result/20260314/pythia6_truth/pythia6_zpole_truth_400k.root`
- `result/20260314/herwig_truth/herwig_ee_zpole_truth_400k.root`
- `result/20260314/sherpa_truth/sherpa_ee_zpole_truth_400k.root`
- existing exploratory X-SCAPE/JETSCAPE controls reused:
  - `result/20260314/xscape_colorless/xscape_epem_colorless_zpole_truth_20k.root`
  - `result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_clean20k.root`

The thrust-axis conversion sweep was then rerun, producing:

- `result/20260314/pythia8_truth/pythia8_zpole_truth_400k_dndy.root`
- `result/20260314/pythia8_rope_truth/pythia8_zpole_truth_rope_400k_dndy.root`
- `result/20260314/pythia8_dire_truth/pythia8_zpole_truth_dire_400k_dndy.root`
- `result/20260314/pythia6_truth/pythia6_zpole_truth_400k_dndy.root`
- `result/20260314/herwig_truth/herwig_ee_zpole_truth_400k_dndy.root`
- `result/20260314/sherpa_truth/sherpa_ee_zpole_truth_400k_dndy.root`

These `dndy` outputs now carry an explicit weak-daughter-veto policy tag.

## Additional technical fixes made during the update

### 1. Herwig converter failure fixed

The local Herwig output could not be consumed robustly through the generic
HepMC3 `ReaderAscii` path, so the comparison chain would otherwise have
silently remained stale or broken for Herwig.

The new direct Ascii parser removes that failure mode.

### 2. Parallel `dN/dy` conversion build race fixed

Concurrent runs of `tools/run_truth_dndy_from_tree.sh` were rebuilding the same
binary in place, causing intermittent:

- `Text file busy`

The wrapper now builds to a temporary binary and atomically renames it into
place, so process-level parallel execution is stable.

## Plots rebuilt

### Beam-axis generator studies

- `result/20260313/generator_dndeta/Generator_dNdEta_Comparison.pdf`
- `result/20260306/top_plots/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf`
- `result/20260306/top_plots/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf`

### Thrust-axis generator studies

- `result/20260314/top_plots_dndy/Generator_dNdY_Comparison.pdf`
- `result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf`

## Analysis note update

Updated in:

- `overleaf_repo_git/main.tex`
- `overleaf_repo_git/appendix_standalone_generator_models.tex`

The note text now states explicitly that the standalone generator comparison
curves use the DELPHI-like weak-daughter-vetoed counting convention.

Copied updated figure assets into:

- `overleaf_repo_git/StrangenessV3/Figures/Generator_dNdEta_Comparison.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/Generator_dNdY_Comparison.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf`
- `overleaf_repo_git/StrangenessV3/Figures/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf`

Rebuilt note:

- `overleaf_repo_git/main.pdf`

Pushed note repo:

- `origin/master: 0b888c4 -> d6dfa49`
- commit: `Apply weak-decay veto to standalone generator plots`

## Paper draft update

Updated in:

- `paper/overleaf_paper_draft/sections/results_and_discussion.tex`

Copied updated figure assets into:

- `paper/overleaf_paper_draft/figures/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf`
- `paper/overleaf_paper_draft/figures/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf`

Rebuilt paper:

- `paper/overleaf_paper_draft/main.pdf`

Pushed paper repo:

- `origin/master: 1dd3f58 -> 6e1ccb2`
- commit: `Refresh standalone generator comparison figures`

## Current repo-state summary

### Already pushed

- analysis note Overleaf repo
- paper draft Overleaf repo

### Still local only

The standalone converter and wrapper source changes in the main analysis repo
have not been selectively committed yet in this follow-up. The reason is simple:

- the main repo still has unrelated pre-existing local dirt
- a clean source-side commit should stage only the standalone-veto migration
  files

So the documentation and note/paper outputs are already synchronized, while the
workspace-side generator-converter implementation still needs a selective main
repo commit if we want the source history itself to be cleanly archived.

## Take-home message

The standalone comparison now matches the DELPHI truth convention more honestly.

Concretely:

- the generator `K/pi` curves are no longer inflated or distorted by counting
  weak charged daughters that DELPHI intentionally excludes
- the activity estimators `dN_ch/deta` and thrust-axis `dN_ch/dy` are now based
  on the same DELPHI-like vetoed charged-particle definition
- the beam-axis and thrust-axis generator comparison figures in the note and
  paper are therefore on a materially cleaner footing than before

The main remaining bookkeeping task is not physics but source control:

- make a selective main-repo commit of the standalone converter updates

