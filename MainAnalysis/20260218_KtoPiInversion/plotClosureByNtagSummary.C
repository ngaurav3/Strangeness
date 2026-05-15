#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include <fstream>
#include <algorithm>
#include <vector>

namespace
{
   struct SpeciesDef
   {
      const char *name;
      const char *gen2d;
      const char *reco2d;
      int color;
   };

   TH1D *IntegratedClosureVsNtag(TH2D *hReco, TH2D *hGen, const char *name)
   {
      TH1D *h = (TH1D *)hReco->ProjectionX(name, 1, hReco->GetNbinsY(), "e");
      TH1D *g = (TH1D *)hGen->ProjectionX(TString::Format("%s_gen", name), 1, hGen->GetNbinsY(), "e");
      h->SetDirectory(nullptr);
      g->SetDirectory(nullptr);
      h->Divide(g);
      delete g;
      return h;
   }

   TH2D *CellClosureMap(TH2D *hReco, TH2D *hGen, const char *name)
   {
      TH2D *h = (TH2D *)hReco->Clone(name);
      h->SetDirectory(nullptr);
      h->Divide(hGen);
      return h;
   }

   TH2D *DisplayCappedMap(TH2D *hIn, const char *name, double minV = 0.0, double maxV = 2.0)
   {
      TH2D *h = (TH2D *)hIn->Clone(name);
      h->SetDirectory(nullptr);
      for (int ix = 1; ix <= h->GetNbinsX(); ix++)
      {
         for (int iy = 1; iy <= h->GetNbinsY(); iy++)
         {
            double v = h->GetBinContent(ix, iy);
            if (v <= 0 || std::isnan(v) || std::isinf(v))
            {
               h->SetBinContent(ix, iy, 0.0);
               continue;
            }
            if (v < minV) v = minV;
            if (v > maxV) v = maxV;
            h->SetBinContent(ix, iy, v);
         }
      }
      h->SetMinimum(minV);
      h->SetMaximum(maxV);
      return h;
   }

   void WriteSpeciesTable(std::ofstream &out, const char *label, TH1D *h)
   {
      out << "\n[" << label << "]\n";
      out << "Bin  NtagCenter  Closure(RecoCorr/Gen)\n";
      for (int i = 1; i <= h->GetNbinsX(); i++)
      {
         double x = h->GetXaxis()->GetBinCenter(i);
         double v = h->GetBinContent(i);
         out << i << "  " << x << "  " << v << "\n";
      }
   }
}

void plotClosureByNtagSummary()
{
   gStyle->SetOptStat(0);
   gStyle->SetPaintTextFormat(".2f");

   TFile *fGen = TFile::Open("output/KtoPi-MC-Gen-Closure.root", "READ");
   TFile *fReco = TFile::Open("output/KtoPi-MC-Reco-Closure.root", "READ");
   if (!fGen || fGen->IsZombie() || !fReco || fReco->IsZombie())
   {
      Error("plotClosureByNtagSummary", "Failed to open closure root files");
      return;
   }

   std::vector<SpeciesDef> S = {
      {"Kaon",   "hKPt",  "hKPtCorrected",  kRed + 1},
      {"Pion",   "hPiPt", "hPiPtCorrected", kBlue + 1},
      {"Proton", "hPPt",  "hPPtCorrected",  kGreen + 2}
   };

   std::vector<TH1D *> hNtag;
   hNtag.reserve(S.size());

   std::ofstream table("Closure_ByNtag_Table.txt");
   table << "Closure diagnostics from MC closure files\n";
   table << "Definition: corrected reco / generator\n";

   TString PDF = "Closure_ByNtag_Summary.pdf";

   TCanvas cIntro("cIntro", "Closure summary", 1000, 700);
   TLatex t;
   t.SetNDC();
   t.SetTextSize(0.040);
   t.DrawLatex(0.08, 0.90, "MC Closure Diagnostics in (N_{tag}^{ch}, p_{T})");
   t.SetTextSize(0.030);
   t.DrawLatex(0.08, 0.83, "Input: output/KtoPi-MC-Gen-Closure.root and output/KtoPi-MC-Reco-Closure.root");
   t.DrawLatex(0.08, 0.78, "Species shown: K, #pi, p");
   t.DrawLatex(0.08, 0.73, "Per-cell closure map: (Reco corrected)/(Gen)");
   t.DrawLatex(0.08, 0.68, "Integrated closure: p_{T}-integrated ratio vs N_{tag}^{ch}");
   t.DrawLatex(0.08, 0.63, "Color scale of maps is clipped to [0, 2] for readability");
   cIntro.Print(PDF + "(");

   for (const auto &sp : S)
   {
      TH2D *hGen2D = dynamic_cast<TH2D *>(fGen->Get(sp.gen2d));
      TH2D *hReco2D = dynamic_cast<TH2D *>(fReco->Get(sp.reco2d));
      if (!hGen2D || !hReco2D)
      {
         Error("plotClosureByNtagSummary", "Missing histograms for %s", sp.name);
         continue;
      }

      TH2D *hMap = CellClosureMap(hReco2D, hGen2D, TString::Format("hMap_%s", sp.name));
      TH2D *hMapDisplay = DisplayCappedMap(hMap, TString::Format("hMapDisplay_%s", sp.name));
      TH1D *hInt = IntegratedClosureVsNtag(hReco2D, hGen2D, TString::Format("hNtag_%s", sp.name));
      hNtag.push_back(hInt);

      WriteSpeciesTable(table, sp.name, hInt);

      TCanvas c(TString::Format("c_%s", sp.name), sp.name, 1200, 750);
      TPad pL("pL", "", 0.00, 0.00, 0.62, 1.00);
      TPad pR("pR", "", 0.62, 0.00, 1.00, 1.00);
      pL.SetRightMargin(0.12);
      pR.SetLeftMargin(0.14);
      pL.Draw();
      pR.Draw();

      pL.cd();
      hMapDisplay->SetTitle(TString::Format("%s closure map;N_{tag}^{ch};p_{T} (GeV/c)", sp.name));
      hMapDisplay->Draw("COLZ");

      pR.cd();
      hInt->SetTitle(TString::Format("%s p_{T}-integrated closure;N_{tag}^{ch};RecoCorr/Gen", sp.name));
      hInt->SetMarkerStyle(20);
      hInt->SetMarkerSize(0.9);
      hInt->SetMarkerColor(sp.color);
      hInt->SetLineColor(sp.color);
      hInt->SetMinimum(0.0);
      hInt->SetMaximum(2.2);
      hInt->Draw("E1");
      TLine line(hInt->GetXaxis()->GetXmin(), 1.0, hInt->GetXaxis()->GetXmax(), 1.0);
      line.SetLineStyle(2);
      line.Draw();

      c.Print(PDF);

      delete hMap;
      delete hMapDisplay;
   }

   TCanvas cComp("cComp", "Integrated closure comparison", 1000, 700);
   TH1D *hFrame = nullptr;
   for (TH1D *h : hNtag)
   {
      if (h == nullptr) continue;
      hFrame = (TH1D *)h->Clone("hFrame");
      break;
   }

   if (hFrame != nullptr)
   {
      hFrame->Reset();
      hFrame->SetTitle("Integrated closure comparison;N_{tag}^{ch};RecoCorr/Gen");
      hFrame->SetMinimum(0.0);
      hFrame->SetMaximum(2.2);
      hFrame->Draw();

      TLegend leg(0.66, 0.74, 0.88, 0.89);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);

      for (size_t i = 0; i < hNtag.size(); i++)
      {
         if (hNtag[i] == nullptr) continue;
         int color = S[i].color;
         hNtag[i]->SetMarkerStyle(20 + (int)i);
         hNtag[i]->SetMarkerColor(color);
         hNtag[i]->SetLineColor(color);
         hNtag[i]->Draw("E1 SAME");
         leg.AddEntry(hNtag[i], S[i].name, "lep");
      }

      TLine line(hFrame->GetXaxis()->GetXmin(), 1.0, hFrame->GetXaxis()->GetXmax(), 1.0);
      line.SetLineStyle(2);
      line.Draw();
      leg.Draw();

      cComp.Print(PDF);
      delete hFrame;
   }

   TCanvas cEnd("cEnd", "End", 1000, 700);
   cEnd.Print(PDF + ")");

   table.close();

   for (TH1D *h : hNtag)
      delete h;

   fGen->Close();
   fReco->Close();
}
