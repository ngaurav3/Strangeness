#include "../include/TruthCountingPolicy.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TProfile.h"
#include "TTree.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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

bool hasLongLivedAncestor(const std::shared_ptr<HepMC3::GenParticle>& particle,
                          TDatabasePDG* db,
                          std::unordered_map<const HepMC3::GenParticle*, int>& memo,
                          std::unordered_map<const HepMC3::GenParticle*, int>& state) {
  if (!particle) return false;
  const auto* key = particle.get();
  auto it = memo.find(key);
  if (it != memo.end()) return (it->second != 0);
  if (state[key] == 1) return false;

  state[key] = 1;
  const auto vtx = particle->production_vertex();
  if (vtx) {
    for (const auto& mother : vtx->particles_in()) {
      if (!mother) continue;
      if (ctauCmFromPdg(db, mother->pid()) > 1.0) {
        memo[key] = 1;
        state[key] = 2;
        return true;
      }
      if (hasLongLivedAncestor(mother, db, memo, state)) {
        memo[key] = 1;
        state[key] = 2;
        return true;
      }
    }
  }

  memo[key] = 0;
  state[key] = 2;
  return false;
}
}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input.hepmc3> <output.root>\n";
    return 1;
  }

  const std::string inputName = argv[1];
  const std::string outputName = argv[2];

  HepMC3::ReaderAscii reader(inputName);
  if (reader.failed()) {
    std::cerr << "Cannot open HepMC3 input: " << inputName << "\n";
    return 1;
  }

  TFile outFile(outputName.c_str(), "RECREATE");
  TTree tree("Events", "Truth events for K/pi vs dNch/deta");
  auto* db = TDatabasePDG::Instance();

  int eventNumber = 0;
  double ecm = 91.2;
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

  int acceptedEvents = 0;
  HepMC3::GenEvent event;
  while (!reader.failed()) {
    reader.read_event(event);
    if (reader.failed()) break;

    ++acceptedEvents;
    eventNumber = event.event_number();
    if (eventNumber <= 0) eventNumber = acceptedEvents;

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

    std::unordered_map<const HepMC3::GenParticle*, int> memo;
    std::unordered_map<const HepMC3::GenParticle*, int> state;

    const auto particles = event.particles();
    for (const auto& p : particles) {
      if (!p) continue;
      if (p->status() != 1) continue;

      const int id = p->pid();
      const auto mom = p->momentum();
      const double pxv = mom.px();
      const double pyv = mom.py();
      const double pzv = mom.pz();
      const double ev = mom.e();
      const double mv = mom.m();
      const double ptv = std::sqrt(pxv * pxv + pyv * pyv);
      const double etav = computeEta(pxv, pyv, pzv);
      const double phiv = std::atan2(pyv, pxv);
      const double q = chargeFromPdg(id);
      const bool weakDaughter = hasLongLivedAncestor(p, db, memo, state);

      ++nFinal;
      pdg.push_back(id);
      status.push_back(1);
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

  reader.close();
  outFile.cd();
  tree.Write();
  hNChEta05.Write();
  hKPiInclusive.Write();
  pKPiVsDNdEta.Write();
  TNamed("weakDecayDaughterPolicy",
         "Default nChEta05/nPiPt0405/nKPt0405 branches veto daughters of ancestors with ctau>1 cm; inclusive bookkeeping retained in *Inclusive branches")
      .Write();
  outFile.Close();

  std::cout << "Wrote truth sample to " << outputName << "\n";
  std::cout << "Accepted events: " << acceptedEvents << "\n";
  return 0;
}
