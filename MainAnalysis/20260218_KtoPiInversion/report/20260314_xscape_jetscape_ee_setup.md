# X-SCAPE / JETSCAPE `e+e-` setup note

This note documents exactly what was needed locally to get a minimal
`e^+e^- -> gamma*/Z -> hadrons` truth-level workflow running with
X-SCAPE / JETSCAPE, and how that output was adapted into the same
`K/pi` vs `dN_ch/deta` analysis format used by the standalone
PYTHIA / HERWIG / SHERPA studies in this repo.

The goal here is reproducibility for another agent, not a polished physics
note.

## Scope and current status

- Upstream package used:
  - `X-SCAPE v2.1`
  - bundled `JETSCAPE v4.0.2`
- Local source trees:
  - [external/X-SCAPE-src](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-src)
  - [external/JETSCAPE-src](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/JETSCAPE-src)
- Local build tree:
  - [external/X-SCAPE-build](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-build)
- What is working now:
  - X-SCAPE builds locally
  - `runJetscape` runs with an `epemGun + hybrid hadronization` steering file
  - final-state hadrons are written as ASCII
  - that ASCII output is converted into a ROOT tree with the same event-level
    branches used by the other standalone generators
  - a standalone `K/pi` vs `dN_ch/deta` plot can be produced from that ROOT tree
- What is not finished as a production comparison:
  - no final X-SCAPE curve has been propagated into the main comparison figures
  - no large-statistics production sample has been completed and validated
  - the current state should be treated as an exploratory engineering path, not
    a finalized generator result
- First follow-up control note:
  - [20260314_xscape_hadronization_control.md](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/report/20260314_xscape_hadronization_control.md)
  - this compares the same `epemGun` setup with `hybrid` and `colorless`
    hadronization and shows that the present low `K/pi` value is specific to
    the `hybrid` choice, not to the basic `e+e-` plumbing

## Important point: what was already in upstream X-SCAPE

The main enabling feature for `e+e^-` was already present upstream.

X-SCAPE contains a dedicated `epemGun` initial-state module:
- [external/X-SCAPE-src/src/initialstate/epemGun.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-src/src/initialstate/epemGun.cc)
- [external/X-SCAPE-src/src/initialstate/epemGun.h](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-src/src/initialstate/epemGun.h)

That module already does the critical `e+e^-`-specific setup:
- `HadronLevel:all = off`
- `PartonLevel:FSR = off`
- `PDF:lepton = off`
- `WeakSingleBoson:ffbar2gmz=on`
- `23:onMode = off`
- `23:onIfAny = 1 2 3 4 5`
- `Beams:idA = 11`
- `Beams:idB = -11`
- beam energy read from XML via `<Hard><epemGun><eCM>`

So the key point is:
- I did **not** patch X-SCAPE source code to invent an `e+e^-` mode
- I used the existing `epemGun` path and built local tooling around it

## What was actually added locally

No upstream X-SCAPE source files were patched in this repo for the `e+e^-`
hard-process definition.

The local changes were wrapper / steering / conversion utilities:

### 1. Run wrapper

- [tools/run_xscape_generate_truth.sh](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/run_xscape_generate_truth.sh)

Purpose:
- writes a temporary user XML
- sets `LD_LIBRARY_PATH` for the local dependency stack
- runs `runJetscape` with explicit user and main XML paths
- converts the final-state hadron ASCII output into a ROOT tree

### 2. Final-state ASCII to ROOT converter

- [tools/xscape_final_state_to_truth_root.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/xscape_final_state_to_truth_root.py)

Purpose:
- reads `*_final_state_hadrons.dat`
- counts charged particles in `|eta| < 0.5`
- counts charged `pi`, `K`, and `p` in `0.4 < pT < 5.0 GeV/c`
- writes a ROOT `Events` tree with:
  - `nChEta05`
  - `nPiPt0405`
  - `nKPt0405`
  - `nPPt0405`
  - `kPiPt0405`

This makes X-SCAPE output look like the other standalone truth trees.

### 3. Standalone plotting script

- [make_xscape_truth_kpi_vs_dndeta.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/make_xscape_truth_kpi_vs_dndeta.py)

Purpose:
- reads the converted ROOT tree
- rebins onto the DELPHI `dN_ch/deta` comparison axis
- produces a standalone X-SCAPE `K/pi` vs `dN_ch/deta` plot

### 4. Minimal standalone XML used for the first direct test

- [tools/xscape_epem_hybrid_zpole.xml](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/xscape_epem_hybrid_zpole.xml)

This was used as the first direct smoke-test steering file before the wrapper
was added.

## Local dependency stack that had to be built

The main build blocker was not X-SCAPE logic; it was local dependency
resolution.

The following local packages were built because the required development
packages were not available in the current environment:

- Boost:
  - [external/boost_1_83_0](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/boost_1_83_0)
  - libraries under:
    - `external/boost_1_83_0/stage/lib`
- GSL:
  - [external/gsl-2.8-install](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/gsl-2.8-install)
- zlib:
  - [external/zlib-1.3.1-install](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/zlib-1.3.1-install)
- HDF5:
  - [external/hdf5-1.14.4-3-install](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/hdf5-1.14.4-3-install)
- PYTHIA8 local install already existed and was reused:
  - [external/pythia8317-src/install](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/pythia8317-src/install)

Relevant cache entries from the final X-SCAPE CMake configuration:
- `GSL_CONFIG = external/gsl-2.8-install/bin/gsl-config`
- `HDF5_CXX_LIB = external/hdf5-1.14.4-3-install/lib/libhdf5_cpp.so`
- `HEPMC_INCLUDE_DIR = /usr/include`
- `PYTHIA8_INCLUDE_DIR = external/pythia8317-src/install/include`
- `PYTHIA8_LIBRARY = external/pythia8317-src/install/lib/libpythia8.so`

From:
- [external/X-SCAPE-build/CMakeCache.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-build/CMakeCache.txt)

## How X-SCAPE was configured and built

The final result was:
- configure succeeded
- `runJetscape` built successfully

Built executable:
- [external/X-SCAPE-build/runJetscape](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-build/runJetscape)

The actual local logic was:

1. Point CMake to the local dependency stack.
2. Use the local PYTHIA8 install.
3. Allow `HEPMC_INCLUDE_DIR=/usr/include` even though a full local HepMC setup
   was not established.

Important nuance:
- `runJetscape` itself was enough for our path because we only needed
  `JetScapeWriterFinalStateHadronsAscii`
- a full HepMC-based output chain was not required for the first `e+e^-`
  truth test

## Why the first runtime failed

The first `runJetscape` smoke test failed with:
- `XML Main file not found/not properly opened! Error code : 3`

Cause:
- `runJetscape` defaults to:
  - `../config/jetscape_user.xml`
  - `../config/jetscape_main.xml`
- those defaults depend on the current working directory

See:
- [external/X-SCAPE-src/examples/runJetscape.cc](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-src/examples/runJetscape.cc)

Fix:
- always run with explicit absolute paths for:
  - the user XML
  - the main XML

That behavior is hard-coded in the wrapper now.

## Minimal `e+e^-` steering choices used locally

The local wrapper writes a user XML with:

- event count:
  - `<nEvents>`
- output:
  - `JetScapeWriterFinalStateHadronsAscii = on`
  - `JetScapeWriterAscii = off`
- random seed
- hard process:
  - `<Hard><epemGun>`
  - `eCM = 91.2`
  - `Tune:ee = 7`
  - `Next:numberShowEvent = 0`
- energy loss / shower in vacuum:
  - `lambdaQCD = 0.2`
  - `Matter`
    - `Q0 = 1.0`
    - `QS = 1.0`
    - `in_vac = 1`
    - `vir_factor = 0.25`
    - `initial_virtuality_pT = 0`
    - `recoil_on = 0`
    - `broadening_on = 0`
    - `brick_med = 0`
- hadronization:
  - `JetHadronization`
    - `name = hybrid`
    - `eCMforHadronization = 0`
    - `had_postprop = 0.0`
    - `part_prop = 0.0`
    - `pythia_decays = on`
    - `tau0Max = 10.0`
    - `reco_hadrons_in_pythia = 1`

This is intentionally minimal:
- it proves the `epemGun + hybrid` path works
- it is not yet a tuned production setup

## How the output is converted into our analysis format

The X-SCAPE final-state hadron writer produces:
- `*_final_state_hadrons.dat`

Example output directory:
- [result/20260314/xscape_hybrid](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid)

Example files:
- [xscape_epem_hybrid_zpole_truth_smoke_final_state_hadrons.dat](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_smoke_final_state_hadrons.dat)
- [xscape_epem_hybrid_zpole_truth_smoke.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_smoke.root)

The converter reads the ASCII format:
- header lines beginning with `#\tEvent`
- particle lines containing:
  - index
  - `pid`
  - `status`
  - `E`
  - `Px`
  - `Py`
  - `Pz`

The converter then computes:

1. `nChEta05`
- count charged particles with `|eta| < 0.5`

2. identified yields in `0.4 < pT < 5.0 GeV/c`
- charged pions
- charged kaons
- charged protons

3. event-level ratio
- `kPiPt0405 = nKPt0405 / nPiPt0405`

This was done to match the standalone truth schema already used by:
- PYTHIA8
- PYTHIA6
- HERWIG
- SHERPA

That alignment is what allows `make_xscape_truth_kpi_vs_dndeta.py` to reuse the
same DELPHI binning strategy as the other generators.

## Exact local run path

Current wrapper entry point:
- [tools/run_xscape_generate_truth.sh](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/run_xscape_generate_truth.sh)

Typical usage:

```bash
tools/run_xscape_generate_truth.sh 1000 \
  result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth.root
```

What it does:

1. checks whether `runJetscape` exists and builds it if needed
2. sets up `LD_LIBRARY_PATH` with local Boost/GSL/HDF5/zlib/PYTHIA8/X-SCAPE libs
3. writes a temporary XML steering file next to the requested output
4. runs:

```bash
external/X-SCAPE-build/runJetscape \
  <temporary_user_xml> \
  external/X-SCAPE-src/config/jetscape_main.xml
```

5. converts:

```bash
<output_prefix>_final_state_hadrons.dat
```

into:

```bash
<output_prefix>.root
```

## What was validated

Validated smoke-test outputs:
- [xscape_epem_hybrid_zpole_truth_smoke.root](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_smoke.root)
- [XscapeTruth_KtoPi_vs_dNdEta.pdf](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/XscapeTruth_KtoPi_vs_dNdEta.pdf)
- [XscapeTruth_KtoPi_vs_dNdEta.txt](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/result/20260314/xscape_hybrid/XscapeTruth_KtoPi_vs_dNdEta.txt)

What was checked:
- `runJetscape` initializes and runs
- final-state hadron ASCII is produced
- converter writes a ROOT tree with the expected branches
- standalone `K/pi(dN_ch/deta)` plotting script runs on that ROOT tree

## Current caveats

This path is not yet production quality.

### 1. Pythia hadronization warnings are present

The hybrid hadronization run emits occasional warnings such as:
- failed decay-channel selection for excited hadrons
- colour/junction retries
- mini-string fragmentation failures

In the current exploratory tests, the run still completes and writes output.
But this needs proper yield accounting before claiming a production-quality
X-SCAPE curve.

### 2. We have not yet propagated a validated X-SCAPE curve into the note

The note and paper currently use the standalone comparison set:
- PYTHIA8 baseline
- PYTHIA8 Ropewalk
- PYTHIA8 Dire
- PYTHIA6
- HERWIG
- SHERPA

X-SCAPE is still separate from that main package.

### 3. The current X-SCAPE observable is truth-level only

There is no detector-response or DELPHI-specific reconstruction layer here.
This is a direct truth comparison only.

### 4. The current converter treats the final-state hadron list as the analysis truth

That is fine for the standalone comparison use case, but if we later want a
more detailed hadronization-systematics study, we may need:
- a direct ROOT writer from X-SCAPE/JETSCAPE
- or a more explicit status selection in the converter

## Practical takeaway for the next agent

If you want to continue the X-SCAPE study, the correct next steps are:

1. Use the wrapper, not a hand-run `runJetscape` command.
2. Produce a moderate-statistics sample first, then inspect:
   - event count
   - inclusive `K/pi`
   - `dN_ch/deta` shape
   - run warnings
3. Only after that, consider wiring X-SCAPE into:
   - [make_kpi_vs_dndeta_comparison_plots.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/make_kpi_vs_dndeta_comparison_plots.py)
   - [make_generator_dndeta_comparison.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/make_generator_dndeta_comparison.py)
4. Do not claim a final physics comparison until the hybrid-hadronization warning rate and output completeness have been quantified.

## References to cite if this path is used scientifically

- JETSCAPE framework:
  - `arXiv:1903.07706`
- X-SCAPE README:
  - [external/X-SCAPE-src/README.md](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/X-SCAPE-src/README.md)
- Hybrid Hadronization:
  - `arXiv:2208.10428`

## Files added in this repo for this effort

- [tools/xscape_epem_hybrid_zpole.xml](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/xscape_epem_hybrid_zpole.xml)
- [tools/run_xscape_generate_truth.sh](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/run_xscape_generate_truth.sh)
- [tools/xscape_final_state_to_truth_root.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/tools/xscape_final_state_to_truth_root.py)
- [make_xscape_truth_kpi_vs_dndeta.py](/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/make_xscape_truth_kpi_vs_dndeta.py)
