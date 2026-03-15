#include <cmath>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "ProgressBar.h"
#include "StrangenessMessenger.h"
#include "TruthCountingPolicy.h"

namespace
{
bool IsChargedPDG(long long pdg)
{
   return TruthCountingPolicy::IsCountedChargedForActivity(pdg);
}
}

int main(int argc, char *argv[])
{
   std::string input = "sample/Strangeness/merged_pythia_v2.5.root";
   std::string output = "output/DNdEtaResponse_Nominal.root";
   bool useCentralEtaNtag = false;
   bool usePIDFiducial = true;
   if (argc > 1)
      input = argv[1];
   if (argc > 2)
      output = argv[2];
   if (argc > 3)
      useCentralEtaNtag = (std::string(argv[3]) == "1");
   if (argc > 4)
      usePIDFiducial = (std::string(argv[4]) != "0");

   const double ptMinYield = 0.4;
   const double ntagPtMin = 0.2;
   const double ptMaxYield = 5.0;
   const int maxNchTag = 60;
   const double pidTrackAbsCosMin = 0.15;
   const double pidTrackAbsCosMax = 0.675;

   TFile *fin = TFile::Open(input.c_str(), "READ");
   if (fin == nullptr || fin->IsZombie())
   {
      std::cerr << "Cannot open input: " << input << std::endl;
      return 1;
   }

   StrangenessTreeMessenger M(*fin, std::string("Tree"));
   if (M.Tree == nullptr)
   {
      std::cerr << "Missing Tree in input" << std::endl;
      return 1;
   }

   const int nbins = maxNchTag + 1;
   const double xmin = -0.5;
   const double xmax = maxNchTag + 0.5;

   TH2D *hResp = new TH2D("hDNdEtaResponse",
                          "dN_{ch}/d#eta(|#eta|<0.5) vs reco N_{tag}^{ch};dN_{ch}/d#eta (true, |#eta|<0.5);N_{tag,reco}^{ch}",
                          nbins, xmin, xmax, nbins, xmin, xmax);
   TH2D *hRespK = (TH2D *)hResp->Clone("hDNdEtaResponseK");
   TH2D *hRespPi = (TH2D *)hResp->Clone("hDNdEtaResponsePi");
   TH2D *hRespP = (TH2D *)hResp->Clone("hDNdEtaResponseP");
   hRespK->Reset();
   hRespPi->Reset();
   hRespP->Reset();

   TH1D *hDNdEtaTrue = new TH1D("hDNdEtaTrue",
                                "True dN_{ch}/d#eta distribution (|#eta|<0.5);dN_{ch}/d#eta (true, |#eta|<0.5);Events",
                                nbins, xmin, xmax);
   TH1D *hKTruedNdEta = new TH1D("hKTruedNdEta",
                                 "Generator-level K yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{K}^{gen}",
                                 nbins, xmin, xmax);
   TH1D *hPiTruedNdEta = new TH1D("hPiTruedNdEta",
                                  "Generator-level #pi yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{#pi}^{gen}",
                                  nbins, xmin, xmax);
   TH1D *hPTruedNdEta = new TH1D("hPTruedNdEta",
                                 "Generator-level p yield vs true dN_{ch}/d#eta;dN_{ch}/d#eta (true, |#eta|<0.5);N_{p}^{gen}",
                                 nbins, xmin, xmax);

   hResp->Sumw2();
   hRespK->Sumw2();
   hRespPi->Sumw2();
   hRespP->Sumw2();
   hDNdEtaTrue->Sumw2();
   hKTruedNdEta->Sumw2();
   hPiTruedNdEta->Sumw2();
   hPTruedNdEta->Sumw2();

   const int n = M.Tree->GetEntries();
   std::cout << "Entries: " << n << std::endl;

   auto passPIDFiducialFromMom = [&](double px, double py, double pz) -> bool
   {
      return TruthCountingPolicy::PassPIDFiducialFromMom(
         px, py, pz, usePIDFiducial, pidTrackAbsCosMin, pidTrackAbsCosMax);
   };

   ProgressBar bar(std::cout, n);
   bar.SetStyle(1);
   const int delta = n / 100 + 1;

   for (int iE = 0; iE < n; ++iE)
   {
      if (iE % delta == 0)
      {
         bar.Update(iE);
         bar.Print();
      }

      M.GetEntry(iE);

      if (M.NReco > STRANGE_MAX_RECO || M.NGen > STRANGE_MAX_GEN)
         continue;

      // Keep the standalone response builder aligned with the main analysis:
      // use the archived event-selection bit rather than recomputing component cuts.
      if (M.PassAll != 1)
         continue;

      int nTagReco = 0;
      for (int i = 0; i < M.NReco; ++i)
      {
         if (M.RecoGoodTrack[i] != 1)
            continue;
         if (M.RecoCharge[i] == 0.0)
            continue;
         const double pxr = M.RecoPx[i];
         const double pyr = M.RecoPy[i];
         const double pzr = M.RecoPz[i];
         const double ptr = std::sqrt(pxr * pxr + pyr * pyr);
         if (ptr < ntagPtMin)
            continue;
         if (useCentralEtaNtag)
         {
            if (ptr <= 0.0)
               continue;
            const double etaReco = std::asinh(pzr / ptr);
            if (std::abs(etaReco) >= 0.5)
               continue;
         }
         ++nTagReco;
      }
      if (nTagReco > maxNchTag)
         nTagReco = maxNchTag;

      int nChEta05 = 0;
      int nKgenEvt = 0;
      int nPigenEvt = 0;
      int nPgenEvt = 0;
      for (int i = 0; i < M.NGen; ++i)
      {
         const long long pdg = M.GenID[i];
         const long long apdg = (pdg >= 0 ? pdg : -pdg);
         if (M.GenStatus[i] != 1)
            continue;
         if (!IsChargedPDG(pdg))
            continue;

         const double px = M.GenPx[i];
         const double py = M.GenPy[i];
         const double pz = M.GenPz[i];
         const double pt = std::sqrt(px * px + py * py);

         if (pt > 0.0)
         {
            const double eta = std::asinh(pz / pt);
            if (std::abs(eta) < 0.5)
               ++nChEta05;
         }

         if (pt < ptMinYield || pt >= ptMaxYield)
            continue;
         if (!passPIDFiducialFromMom(px, py, pz))
            continue;
         if (apdg == 321)
            ++nKgenEvt;
         if (apdg == 211)
            ++nPigenEvt;
         if (apdg == 2212)
            ++nPgenEvt;
      }

      if (nChEta05 > maxNchTag)
         nChEta05 = maxNchTag;
      const double dNdEtaTrue = static_cast<double>(nChEta05);

      hResp->Fill(dNdEtaTrue, nTagReco);
      hRespK->Fill(dNdEtaTrue, nTagReco, nKgenEvt);
      hRespPi->Fill(dNdEtaTrue, nTagReco, nPigenEvt);
      hRespP->Fill(dNdEtaTrue, nTagReco, nPgenEvt);
      hDNdEtaTrue->Fill(dNdEtaTrue);
      hKTruedNdEta->Fill(dNdEtaTrue, nKgenEvt);
      hPiTruedNdEta->Fill(dNdEtaTrue, nPigenEvt);
      hPTruedNdEta->Fill(dNdEtaTrue, nPgenEvt);
   }

   bar.Update(n);
   bar.Print();
   std::cout << std::endl;

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   hResp->Write();
   hRespK->Write();
   hRespPi->Write();
   hRespP->Write();
   hDNdEtaTrue->Write();
   hKTruedNdEta->Write();
   hPiTruedNdEta->Write();
   hPTruedNdEta->Write();
   fout->Close();

   fin->Close();
   std::cout << "Wrote: " << output << std::endl;
   return 0;
}
