# 2026-03-19 v7 method-text synchronization

## Scope
This note records the text-only synchronization pass performed after the v7 Yi-style result was promoted to the quoted DELPHI activity-based result in the analysis note and paper draft.

## Main changes
- Explicitly distinguished the two workflows:
  - legacy v6 factorized chain
  - current quoted v7 Yi-style flat-mu chain
- Rewrote the thrust-axis `dN_ch/dy` method text so the quoted result is described as the v7 Yi-style procedure, while the earlier Bayes/SVD same-variable package is described only as legacy validation material.
- Rewrote the paper `analysis` section so it no longer describes the legacy factorized chain as the active quoted result.
- Updated the paper `observable_definition` section to state clearly that the quoted activity-based results are now:
  - beam-axis `dN_ch/deta`
  - thrust-axis `dN_ch/dy`
  while `N_tag^ch` is retained only as a legacy reference branch.
- Updated `results_and_discussion` and `summary` in the paper so the method language matches the promoted v7 result.
- Updated the early workflow section of the analysis note so the factorized correction and unfolding description is labeled explicitly as the legacy v6 reference workflow.

## Files updated
- `overleaf_repo_git/main.tex`
- `paper/overleaf_paper_draft/sections/analysis.tex`
- `paper/overleaf_paper_draft/sections/observable_definition.tex`
- `paper/overleaf_paper_draft/sections/results_and_discussion.tex`
- `paper/overleaf_paper_draft/sections/summary.tex`

## Verification
- Rebuilt analysis note PDF successfully.
- Rebuilt paper draft PDF successfully.
- Grep audit confirmed that remaining references to factorized unfolding or same-variable Bayes/SVD treatment are now explicitly marked as legacy v6 context rather than active quoted-result methodology.
