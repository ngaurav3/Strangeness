#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DIR="$ROOT_DIR/external/pythia6428-src"
MAIN_SRC="$ROOT_DIR/tools/pythia6_generate_truth_kpi_dndeta.f90"
EXE="$ROOT_DIR/tools/pythia6_generate_truth.exe"
ROOT_SETUP="/raid5/root/root-v6.34.04/root/bin/thisroot.sh"

NEV="${1:-10000}"
OUTFILE="${2:-$ROOT_DIR/result/20260313/pythia6_truth/pythia6_zpole_truth_events.txt}"
OUTROOT="${3:-}"
OUTDAT="${OUTFILE%.txt}_final_state_particles.dat"

mkdir -p "$(dirname "$OUTFILE")"

gfortran -O2 -std=legacy -fallow-argument-mismatch \
  "$SRC_DIR/pythia6428.f" "$MAIN_SRC" -o "$EXE"

if [[ -n "$OUTROOT" ]]; then
  "$EXE" "$NEV" "$OUTFILE" "$OUTDAT"
  if [[ -f "$ROOT_SETUP" ]]; then
    # shellcheck disable=SC1090
    source "$ROOT_SETUP"
  fi
  python3 "$ROOT_DIR/tools/pythia6_final_state_to_truth_root.py" "$OUTDAT" "$OUTROOT"
else
  "$EXE" "$NEV" "$OUTFILE"
fi
