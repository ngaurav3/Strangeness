#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/external/X-SCAPE-build"
SOURCE_DIR="$ROOT_DIR/external/X-SCAPE-src"
MAIN_XML="$SOURCE_DIR/config/jetscape_main.xml"
ROOT_SETUP="/raid5/root/root-v6.34.04/root/bin/thisroot.sh"

NEVT="${1:-10000}"
OUT_ROOT_INPUT="${2:-$ROOT_DIR/result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth.root}"
SEED="${3:-12345}"
HADRONIZATION_NAME="${XSCAPE_HADRONIZATION_NAME:-hybrid}"
PYTHIA_DECAYS="${XSCAPE_PYTHIA_DECAYS:-on}"
TAU0_MAX="${XSCAPE_TAU0_MAX:-10.0}"
RECO_HADRONS_IN_PYTHIA="${XSCAPE_RECO_HADRONS_IN_PYTHIA:-1}"
ECM_HADRONIZATION="${XSCAPE_ECM_HADRONIZATION:-0}"
HAD_POSTPROP="${XSCAPE_HAD_POSTPROP:-0.0}"
PART_PROP="${XSCAPE_PART_PROP:-0.0}"

case "$OUT_ROOT_INPUT" in
  /*) OUT_ROOT="$OUT_ROOT_INPUT" ;;
  *) OUT_ROOT="$ROOT_DIR/$OUT_ROOT_INPUT" ;;
esac

OUT_DIR="$(dirname "$OUT_ROOT")"
OUT_PREFIX="${OUT_ROOT%.root}"
TMP_XML="$OUT_DIR/$(basename "${OUT_PREFIX}").xml"
OUT_DAT="${OUT_PREFIX}_final_state_hadrons.dat"

mkdir -p "$OUT_DIR"

if [[ ! -x "$BUILD_DIR/runJetscape" ]]; then
  make -C "$BUILD_DIR" -j8 runJetscape
fi

if [[ -f "$ROOT_SETUP" ]]; then
  # shellcheck disable=SC1090
  source "$ROOT_SETUP"
fi

export LD_LIBRARY_PATH="$ROOT_DIR/external/gsl-2.8-install/lib:$ROOT_DIR/external/hdf5-1.14.4-3-install/lib:$ROOT_DIR/external/zlib-1.3.1-install/lib:$ROOT_DIR/external/boost_1_83_0/stage/lib:$ROOT_DIR/external/pythia8317-src/install/lib:$BUILD_DIR/src/lib:$BUILD_DIR/lib:$BUILD_DIR/external_packages/gtl/lib:$BUILD_DIR/external_packages/lib:${LD_LIBRARY_PATH:-}"

cat > "$TMP_XML" <<EOF
<?xml version="1.0"?>
<jetscape>
  <nEvents> ${NEVT} </nEvents>

  <outputFilename>${OUT_PREFIX}</outputFilename>
  <JetScapeWriterAscii> off </JetScapeWriterAscii>
  <JetScapeWriterFinalStateHadronsAscii> on </JetScapeWriterFinalStateHadronsAscii>

  <Random>
    <seed> ${SEED} </seed>
  </Random>

  <Hard>
    <epemGun>
      <name>epemGun</name>
      <eCM>91.2</eCM>
      <LinesToRead>
        Tune:ee = 7
        Next:numberShowEvent = 0
      </LinesToRead>
    </epemGun>
  </Hard>

  <Eloss>
    <lambdaQCD>0.2</lambdaQCD>
    <Matter>
      <Q0>1.0</Q0>
      <QS>1.0</QS>
      <in_vac>1</in_vac>
      <vir_factor>0.25</vir_factor>
      <initial_virtuality_pT>0</initial_virtuality_pT>
      <recoil_on>0</recoil_on>
      <broadening_on>0</broadening_on>
      <brick_med>0</brick_med>
    </Matter>
  </Eloss>

  <JetHadronization>
    <name>${HADRONIZATION_NAME}</name>
    <eCMforHadronization>${ECM_HADRONIZATION}</eCMforHadronization>
    <had_postprop>${HAD_POSTPROP}</had_postprop>
    <part_prop>${PART_PROP}</part_prop>
    <pythia_decays>${PYTHIA_DECAYS}</pythia_decays>
    <tau0Max>${TAU0_MAX}</tau0Max>
    <reco_hadrons_in_pythia>${RECO_HADRONS_IN_PYTHIA}</reco_hadrons_in_pythia>
  </JetHadronization>
</jetscape>
EOF

"$BUILD_DIR/runJetscape" "$TMP_XML" "$MAIN_XML"
python3 "$ROOT_DIR/tools/xscape_final_state_to_truth_root.py" "$OUT_DAT" "$OUT_ROOT"

echo "X-SCAPE final-state hadrons: $OUT_DAT"
echo "X-SCAPE truth ROOT:         $OUT_ROOT"
echo "X-SCAPE hadronization:      $HADRONIZATION_NAME"
echo "X-SCAPE pythia_decays:      $PYTHIA_DECAYS"
echo "X-SCAPE tau0Max [mm/c]:     $TAU0_MAX"
echo "X-SCAPE reco_hadrons flag:  $RECO_HADRONS_IN_PYTHIA"
