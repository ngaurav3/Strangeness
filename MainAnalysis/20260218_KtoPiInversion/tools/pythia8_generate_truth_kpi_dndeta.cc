#include "Pythia8/Pythia.h"
#include "../include/TruthCountingPolicy.h"

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TProfile.h"
#include "TTree.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Pythia8;

namespace {
constexpr double SPEED_OF_LIGHT_CM_PER_S = 2.99792458e10;

double computeEta(double px, double py, double pz) {
  const double pt = std::sqrt(px * px + py * py);
  if (pt <= 0.) return 0.;
  return std::asinh(pz / pt);
}

double chargeFromPdg(int pdg) {
  return TruthCountingPolicy::CountedChargeFromPdg(pdg);
}

double ctauCmFromPdg(TDatabasePDG* db, int pdg) {
  static std::unordered_map<int, double> cache;
  const int apdg = std::abs(pdg);
  auto it = cache.find(apdg);
  if (it != cache.end()) return it->second;

  double ctau = 0.0;
  if (db != nullptr) {
    if (auto* particle = db->GetParticle(apdg)) {
      const double lifetime = particle->Lifetime();
      if (lifetime > 0.0 && std::isfinite(lifetime)) {
        ctau = lifetime * SPEED_OF_LIGHT_CM_PER_S;
      }
    }
  }
  cache[apdg] = ctau;
  return ctau;
}

bool hasLongLivedAncestor(const Event& event,
                          int idx,
                          TDatabasePDG* db,
                          std::vector<int>& memo,
                          std::vector<int>& state) {
  if (idx <= 0 || idx >= event.size()) return false;
  if (memo[idx] != -1) return (memo[idx] != 0);
  if (state[idx] == 1) return false;

  state[idx] = 1;
  const auto mothers = event[idx].motherList();
  for (int mother : mothers) {
    if (mother <= 0 || mother >= event.size()) continue;
    if (ctauCmFromPdg(db, event[mother].id()) > 1.0) {
      memo[idx] = 1;
      state[idx] = 2;
      return true;
    }
    if (hasLongLivedAncestor(event, mother, db, memo, state)) {
      memo[idx] = 1;
      state[idx] = 2;
      return true;
    }
  }

  memo[idx] = 0;
  state[idx] = 2;
  return false;
}
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <nEvents> <output.root> [baseline|rope|dire]\n";
    return 1;
  }

  const int nEvents = std::atoi(argv[1]);
  const std::string outName = argv[2];
  const std::string mode =
      (argc >= 4 ? std::string(argv[3]) :
       (std::getenv("PYTHIA8_MODE") ? std::getenv("PYTHIA8_MODE")
                                    : "baseline"));
  if (nEvents <= 0) {
    std::cerr << "Number of events must be positive.\n";
    return 1;
  }

  Pythia pythia;
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:eCM = 91.2");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  if (mode == "rope") {
    pythia.readString("Ropewalk:RopeHadronization = on");
    pythia.readString("Ropewalk:doShoving = off");
    pythia.readString("Ropewalk:doFlavour = on");
    pythia.readString("Ropewalk:r0 = 0.5");
    pythia.readString("Ropewalk:m0 = 0.2");
    pythia.readString("Ropewalk:beta = 0.1");
    pythia.readString("PartonVertex:setVertex = on");
  } else if (mode == "dire") {
    pythia.readString("PartonShowers:model = 3");
  } else if (mode != "baseline") {
    std::cerr << "Unknown PYTHIA8 mode: " << mode
              << " (expected baseline, rope, or dire)\n";
    return 1;
  }

  if (!pythia.init()) {
    std::cerr << "PYTHIA initialization failed.\n";
    return 1;
  }
  TDatabasePDG* db = TDatabasePDG::Instance();

  TFile outFile(outName.c_str(), "RECREATE");
  TTree tree("Events", "PYTHIA8 truth events for K/pi vs dNch/deta");

  int eventNumber = 0;
  double ecm = 91.2;
  std::string generatorMode = mode;
  int nFinal = 0;
  int nCharged = 0;
  int nChargedInclusive = 0;
  int nChEta05 = 0;
  int nChEta05Inclusive = 0;
  int nPiPt0405 = 0;
  int nPiPt0405Inclusive = 0;
  int nKPt0405 = 0;
  int nKPt0405Inclusive = 0;
  double kPiPt0405 = -1.;
  double kPiPt0405Inclusive = -1.;

  std::vector<int> pdg;
  std::vector<int> status;
  std::vector<int> isWeakDecayDaughter;
  std::vector<float> charge;
  std::vector<float> px;
  std::vector<float> py;
  std::vector<float> pz;
  std::vector<float> e;
  std::vector<float> m;
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;

  tree.Branch("eventNumber", &eventNumber);
  tree.Branch("ecm", &ecm);
  tree.Branch("generatorMode", &generatorMode);
  tree.Branch("nFinal", &nFinal);
  tree.Branch("nCharged", &nCharged);
  tree.Branch("nChargedInclusive", &nChargedInclusive);
  tree.Branch("nChEta05", &nChEta05);
  tree.Branch("nChEta05Inclusive", &nChEta05Inclusive);
  tree.Branch("nPiPt0405", &nPiPt0405);
  tree.Branch("nPiPt0405Inclusive", &nPiPt0405Inclusive);
  tree.Branch("nKPt0405", &nKPt0405);
  tree.Branch("nKPt0405Inclusive", &nKPt0405Inclusive);
  tree.Branch("kPiPt0405", &kPiPt0405);
  tree.Branch("kPiPt0405Inclusive", &kPiPt0405Inclusive);
  tree.Branch("pdg", &pdg);
  tree.Branch("status", &status);
  tree.Branch("isWeakDecayDaughter", &isWeakDecayDaughter);
  tree.Branch("charge", &charge);
  tree.Branch("px", &px);
  tree.Branch("py", &py);
  tree.Branch("pz", &pz);
  tree.Branch("e", &e);
  tree.Branch("m", &m);
  tree.Branch("pt", &pt);
  tree.Branch("eta", &eta);
  tree.Branch("phi", &phi);

  TH1D hNChEta05("hNChEta05",
                 "Charged multiplicity in |#eta|<0.5;dN_{ch}/d#eta (|#eta|<0.5);Events",
                 60, -0.5, 59.5);
  TH1D hKPiInclusive("hKPiInclusive",
                     "Inclusive K/#pi (0.4 < p_{T} < 5.0 GeV/c);K/#pi;Events",
                     100, 0.0, 1.0);
  TProfile pKPiVsDNdEta("pKPiVsDNdEta",
                        "K/#pi vs dN_{ch}/d#eta;dN_{ch}/d#eta (|#eta|<0.5);K/#pi",
                        20, -0.5, 39.5);

  int nAccepted = 0;
  for (eventNumber = 0; eventNumber < nEvents; ++eventNumber) {
    if (!pythia.next()) continue;
    ++nAccepted;

    pdg.clear();
    status.clear();
    isWeakDecayDaughter.clear();
    charge.clear();
    px.clear();
    py.clear();
    pz.clear();
    e.clear();
    m.clear();
    pt.clear();
    eta.clear();
    phi.clear();

    nFinal = 0;
    nCharged = 0;
    nChargedInclusive = 0;
    nChEta05 = 0;
    nChEta05Inclusive = 0;
    nPiPt0405 = 0;
    nPiPt0405Inclusive = 0;
    nKPt0405 = 0;
    nKPt0405Inclusive = 0;
    kPiPt0405 = -1.;
    kPiPt0405Inclusive = -1.;

    std::vector<int> memo(pythia.event.size(), -1);
    std::vector<int> state(pythia.event.size(), 0);

    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle& part = pythia.event[i];
      if (!part.isFinal()) continue;

      ++nFinal;

      const int id = part.id();
      const double q = chargeFromPdg(id);
      const double pxv = part.px();
      const double pyv = part.py();
      const double pzv = part.pz();
      const double ev = part.e();
      const double mv = part.m();
      const double ptv = std::sqrt(pxv * pxv + pyv * pyv);
      const double etav = computeEta(pxv, pyv, pzv);
      const double phiv = std::atan2(pyv, pxv);
      const bool weakDaughter = hasLongLivedAncestor(pythia.event, i, db, memo, state);

      pdg.push_back(id);
      status.push_back(part.status());
      isWeakDecayDaughter.push_back(weakDaughter ? 1 : 0);
      charge.push_back(static_cast<float>(q));
      px.push_back(static_cast<float>(pxv));
      py.push_back(static_cast<float>(pyv));
      pz.push_back(static_cast<float>(pzv));
      e.push_back(static_cast<float>(ev));
      m.push_back(static_cast<float>(mv));
      pt.push_back(static_cast<float>(ptv));
      eta.push_back(static_cast<float>(etav));
      phi.push_back(static_cast<float>(phiv));

      if (q == 0.) continue;
      ++nChargedInclusive;
      if (!weakDaughter) ++nCharged;
      if (std::abs(etav) < 0.5) {
        ++nChEta05Inclusive;
        if (!weakDaughter) ++nChEta05;
      }

      if (TruthCountingPolicy::IsCountedPionForRatio(id, pxv, pyv, pzv)) {
        ++nPiPt0405Inclusive;
        if (!weakDaughter) ++nPiPt0405;
      }
      if (TruthCountingPolicy::IsCountedKaonForRatio(id, pxv, pyv, pzv)) {
        ++nKPt0405Inclusive;
        if (!weakDaughter) ++nKPt0405;
      }
    }

    if (nPiPt0405 > 0) kPiPt0405 = double(nKPt0405) / double(nPiPt0405);
    if (nPiPt0405Inclusive > 0) {
      kPiPt0405Inclusive = double(nKPt0405Inclusive) / double(nPiPt0405Inclusive);
    }

    tree.Fill();
    hNChEta05.Fill(nChEta05);
    if (kPiPt0405 >= 0.) {
      hKPiInclusive.Fill(kPiPt0405);
      pKPiVsDNdEta.Fill(nChEta05, kPiPt0405);
    }
  }

  pythia.stat();

  outFile.cd();
  tree.Write();
  hNChEta05.Write();
  hKPiInclusive.Write();
  pKPiVsDNdEta.Write();
  TNamed("weakDecayDaughterPolicy",
         "Default nChEta05/nPiPt0405/nKPt0405 branches veto daughters of ancestors with ctau>1 cm; inclusive bookkeeping retained in *Inclusive branches")
      .Write();
  outFile.Close();

  std::cout << "Wrote truth-level sample to " << outName << "\n";
  std::cout << "Accepted events: " << nAccepted << "\n";
  return 0;
}
