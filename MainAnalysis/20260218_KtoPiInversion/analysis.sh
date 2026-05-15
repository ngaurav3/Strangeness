#!/usr/bin/env bash
set -euo pipefail

#------------------------------------------------------------
# MC closure workflow (for closure test only)
#------------------------------------------------------------
./ExecuteKtoPiAnalysis \
  --IsGen true \
  --UseMCTruthMatrix true \
  --UsePIDFiducial true \
  --PtMin 0.4 \
  --NtagPtMin 0.2 \
  --Input sample/Strangeness/merged_pythia8ext_v2.6.root \
  --Output output/KtoPi-MC-Gen-Closure.root

./ExecuteKtoPiAnalysis \
  --IsGen false \
  --UseMCTruthMatrix true \
  --UsePIDFiducial true \
  --PtMin 0.4 \
  --NtagPtMin 0.2 \
  --Input sample/Strangeness/merged_pythia8ext_v2.6.root \
  --Output output/KtoPi-MC-Reco-Closure.root

#------------------------------------------------------------
# Nominal workflow (data/MC comparison)
#------------------------------------------------------------
./ExecuteKtoPiAnalysis \
  --IsGen true \
  --UseMCTruthMatrix true \
  --UsePIDFiducial true \
  --PtMin 0.4 \
  --NtagPtMin 0.2 \
  --Input sample/Strangeness/merged_pythia8ext_v2.6.root \
  --Output output/KtoPi-MC-Gen-Nominal.root

./ExecuteKtoPiAnalysis \
  --IsGen false \
  --UseMCTruthMatrix false \
  --UsePIDFiducial true \
  --PtMin 0.4 \
  --NtagPtMin 0.2 \
  --Input sample/Strangeness/merged_pythia8ext_v2.6.root \
  --Output output/KtoPi-MC-Reco-Nominal.root

./ExecuteKtoPiAnalysis \
  --IsGen false \
  --UseMCTruthMatrix false \
  --UsePIDFiducial true \
  --PtMin 0.4 \
  --NtagPtMin 0.2 \
  --Input sample/Strangeness/merged_data_v2.6.root \
  --Output output/KtoPi-Data-Reco-Nominal.root
