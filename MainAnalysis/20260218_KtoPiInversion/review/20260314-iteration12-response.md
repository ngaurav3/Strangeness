# 2026-03-14 Iteration 12 Response

## Scope

This note records the implementation work done in response to
[20260314-iteration12.md](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/review/20260314-iteration12.md).

The reviewer concern was definition consistency across:

- DELPHI truth activity for `dN_ch/deta`
- DELPHI truth activity for thrust-axis `dN_ch/dy`
- standalone generator `K/pi` overlays
- standalone generator activity overlays

## Reviewer Findings Accepted

The reviewer findings were correct:

1. The charged-particle activity definition was not one common predicate across
   the active DELPHI truth path and all standalone generator converters.
2. The standalone `K/pi` curves did not yet apply the DELPHI PID fiducial used
   by the DELPHI truth branch.
3. The note wording had become stronger than the actual implementation.

## Implementation Strategy

The fix was to make the truth policy explicit and shared instead of letting it
drift separately in each standalone converter.

The shared policy now has three pieces:

- `is_counted_charged_for_activity`
- `is_counted_pion_for_ratio`
- `is_counted_kaon_for_ratio`

The activity definition follows the broader DELPHI `IsChargedPDG` choice rather
than the narrower light-charged-only list.

That means the counted charged-particle list is:

- `e`
- `mu`
- `tau`
- `pi`
- `K`
- `p`
- `Sigma`
- `Xi`
- `Omega`
- `D`
- `Ds`
- `B`
- `Bc`
- `W`

with the exact PDG implementation matching the DELPHI truth code.

For the ratio definition, the counted kaons and pions now follow the same truth
fiducial used by the DELPHI branch:

- charged species only
- `0.4 <= pT < 5.0 GeV/c`
- nominal DELPHI PID fiducial

The standalone weak-daughter veto remains:

- charged daughters of ancestors with `ctau > 1 cm` are removed from both the
  activity count and the counted kaon/pion yields

## Code Changes

Shared policy added:

- [include/TruthCountingPolicy.h](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/include/TruthCountingPolicy.h)
- [tools/truth_counting_policy.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/truth_counting_policy.py)

Active DELPHI truth path updated:

- [KtoPiAnalysis.cpp](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/KtoPiAnalysis.cpp)
- [BuildDNdEtaResponse.cpp](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/BuildDNdEtaResponse.cpp)

Standalone truth builders / converters updated:

- [tools/pythia8_generate_truth_kpi_dndeta.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/pythia8_generate_truth_kpi_dndeta.cc)
- [tools/hepmc3_to_truth_root.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/hepmc3_to_truth_root.cc)
- [tools/herwig_ascii_to_truth_root.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/herwig_ascii_to_truth_root.py)
- [tools/herwig_hepmc_to_truth_root.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/herwig_hepmc_to_truth_root.cc)
- [tools/pythia6_generate_truth_kpi_dndeta.f90](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/pythia6_generate_truth_kpi_dndeta.f90)
- [tools/pythia6_final_state_to_truth_root.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/pythia6_final_state_to_truth_root.py)
- [tools/xscape_final_state_to_truth_root.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/xscape_final_state_to_truth_root.py)
- [tools/build_truth_kpi_vs_dndy_from_tree.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/build_truth_kpi_vs_dndy_from_tree.cc)

## Charge-Branch Correction

One subtle follow-up fix was needed while addressing the review:

- the stored `charge` branch in some standalone truth trees had been using a
  sign shortcut based on the PDG sign

That is not correct for particles such as:

- `e-`
- `Sigma-`
- `Xi-`
- `Omega-`

The shared policy therefore also now defines the explicit PDG-to-charge mapping
used by the standalone converters.

This does not change the species or activity predicate itself, but it does make
the stored truth trees internally correct and keeps the thrust-axis `dN/dy`
conversion consistent.

## Regeneration Status

The standalone truth trees and the beam-axis / thrust-axis generator comparison
plots were regenerated from the updated shared-policy outputs after the code
changes above.

The affected note and paper comparison sections were then rebuilt so the text
and figures describe the same implemented definition.

The concrete regenerated comparison outputs are:

- beam-axis generator shape:
  - [result/20260313/generator_dndeta/Generator_dNdEta_Comparison.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260313/generator_dndeta/Generator_dNdEta_Comparison.pdf)
- beam-axis \(K/\pi\) overlays:
  - [result/20260306/top_plots/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260306/top_plots/KtoPi_vs_dNdEta_DELPHI_vs_ALICE.pdf)
  - [result/20260306/top_plots/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260306/top_plots/KtoPi_vs_dNdEta_DELPHI_vs_Generators.pdf)
  - [result/20260306/top_plots/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260306/top_plots/KtoPi_vs_dNdEta_DataMC_with_Systematics.pdf)
- thrust-axis generator shape:
  - [result/20260314/top_plots_dndy/Generator_dNdY_Comparison.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/top_plots_dndy/Generator_dNdY_Comparison.pdf)
- thrust-axis \(K/\pi\) overlays:
  - [result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DELPHI_vs_ALICE.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DELPHI_vs_ALICE.pdf)
  - [result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DELPHI_vs_Generators.pdf)
  - [result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DataMC_with_Systematics.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/top_plots_dndy/KtoPi_vs_dNdY_DataMC_with_Systematics.pdf)

The clean PYTHIA6 truth repack used for the thrust-axis regeneration was also
verified explicitly:

- [result/20260314/pythia6_truth/pythia6_zpole_truth_400k.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/pythia6_truth/pythia6_zpole_truth_400k.root)
  contains `400000` `Events` entries
- the derived thrust-axis summary file is:
  - [result/20260314/pythia6_truth/pythia6_zpole_truth_400k_dndy.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/pythia6_truth/pythia6_zpole_truth_400k_dndy.root)

## Documentation Update

The note and paper wording were tightened to describe the actual implementation:

- standalone activity counts use the same broader DELPHI charged-particle list
- standalone `K/pi` uses the DELPHI PID fiducial
- weak-daughter removal is explicit where ancestry is available
- X-SCAPE/JETSCAPE remains hadron-level input, so the weak-daughter part is
  effectively implicit there

The rebuilt document entry points are:

- note source:
  - [overleaf_repo_git/main.tex](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/overleaf_repo_git/main.tex)
- note PDF:
  - [overleaf_repo_git/main.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/overleaf_repo_git/main.pdf)
- paper source:
  - [paper/overleaf_paper_draft/sections/results_and_discussion.tex](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/paper/overleaf_paper_draft/sections/results_and_discussion.tex)
- paper PDF:
  - [paper/overleaf_paper_draft/main.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/paper/overleaf_paper_draft/main.pdf)

Document update / push status:

- analysis note Overleaf:
  - local commit `87dbe47`
  - remote update `68f3ca6 -> 87dbe47`
- paper draft Overleaf:
  - local commit `a4b1744`
  - remote update `6e1ccb2 -> a4b1744`

## Bottom Line

Iteration 12 changes the repo from “similar but not identical” truth definitions
to one shared DELPHI-style truth counting policy across the active DELPHI truth
path and the standalone generator overlays used in the analysis note and paper.
