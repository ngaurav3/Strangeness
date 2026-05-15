#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include <algorithm>

namespace
{
   void ConfigureUpperPad(TPad *p)
   {
      p->SetBottomMargin(0.00);
      p->SetTopMargin(0.08);
      p->SetLeftMargin(0.14);
      p->SetRightMargin(0.05);
   }

   void ConfigureLowerPad(TPad *p)
   {
      p->SetTopMargin(0.00);
      p->SetBottomMargin(0.32);
      p->SetLeftMargin(0.14);
      p->SetRightMargin(0.05);
   }

   void ConfigureUpperAxis(TH1 *h)
   {
      h->GetYaxis()->SetTitleSize(0.058);
      h->GetYaxis()->SetLabelSize(0.046);
      h->GetYaxis()->SetTitleOffset(1.00);
      h->GetXaxis()->SetLabelSize(0.0);
      h->GetXaxis()->SetTitleSize(0.0);
      h->GetYaxis()->ChangeLabel(1, -1, 0.0, -1, -1, -1, " ");
   }

   void ConfigureRatioAxis(TH1 *h)
   {
      h->GetXaxis()->SetTitleSize(0.13);
      h->GetYaxis()->SetTitleSize(0.12);
      h->GetXaxis()->SetLabelSize(0.11);
      h->GetYaxis()->SetLabelSize(0.10);
      h->GetXaxis()->SetTitleOffset(0.95);
      h->GetYaxis()->SetTitleOffset(0.55);
      h->GetXaxis()->SetLabelOffset(0.02);
      h->GetYaxis()->SetNdivisions(505);
      h->GetYaxis()->ChangeLabel(-1, -1, 0.0, -1, -1, -1, " ");
   }

   void MakeSpeciesClosure(TFile *fGen, TFile *fReco,
      const char *gen2DName, const char *reco2DName,
      const char *title, const char *outBase,
      double yMin = 1.0e2, double yMax = 1.0e6)
   {
      TH2D *hGen2D = dynamic_cast<TH2D *>(fGen->Get(gen2DName));
      TH2D *hReco2D = dynamic_cast<TH2D *>(fReco->Get(reco2DName));

      if (!hGen2D || !hReco2D)
      {
         Error("MakeSpeciesClosure", "Missing 2D histograms: %s or %s", gen2DName, reco2DName);
         return;
      }

      TH1D *hGen = hGen2D->ProjectionY(TString::Format("%s_gen", outBase), 1, hGen2D->GetNbinsX(), "e");
      TH1D *hReco = hReco2D->ProjectionY(TString::Format("%s_reco", outBase), 1, hReco2D->GetNbinsX(), "e");
      hGen->SetDirectory(nullptr);
      hReco->SetDirectory(nullptr);

      TH1D *hRatio = (TH1D *)hReco->Clone(TString::Format("%s_ratio", outBase));
      hRatio->SetDirectory(nullptr);
      hRatio->Divide(hGen);

	      hReco->SetMarkerStyle(20);
	      hReco->SetMarkerSize(1.4);
      hReco->SetMarkerColor(kRed + 1);
      hReco->SetLineColor(kRed + 1);

      hGen->SetLineColor(kBlue + 1);
      hGen->SetLineWidth(2);
      hGen->SetLineStyle(2);

      TCanvas *c = new TCanvas(TString::Format("c_%s", outBase), title, 820, 800);
      TPad *p1 = new TPad("p1", "", 0, 0.30, 1, 1);
      TPad *p2 = new TPad("p2", "", 0, 0.00, 1, 0.30);
      ConfigureUpperPad(p1);
      ConfigureLowerPad(p2);
      p1->Draw();
      p2->Draw();

      p1->cd();
      p1->SetLogy();
      hReco->SetTitle(TString::Format("%s; p_{T} (GeV/c); Yield", title));
      hReco->GetXaxis()->SetRangeUser(0.4, 10.0);
      ConfigureUpperAxis(hReco);
      hReco->SetMinimum(yMin);
      hReco->SetMaximum(yMax);
      hReco->Draw("E1");
      hGen->Draw("HIST SAME");

      TLegend *leg = new TLegend(0.52, 0.72, 0.86, 0.86);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(hReco, "MC reco corrected", "lep");
      leg->AddEntry(hGen, "MC generator", "l");
      leg->Draw();

      p2->cd();
      hRatio->SetTitle(";p_{T} (GeV/c);Closure");
      hRatio->GetXaxis()->SetRangeUser(0.4, 10.0);
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerSize(1.2);
      hRatio->SetMinimum(0.5);
      hRatio->SetMaximum(1.5);
      ConfigureRatioAxis(hRatio);
      hRatio->Draw("E1");

      TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw();
      hRatio->Draw("E1 SAME");

      c->SaveAs(TString::Format("%s.png", outBase));
      c->SaveAs(TString::Format("%s.pdf", outBase));

      delete c;
      delete hRatio;
      delete hReco;
      delete hGen;
   }
}

void plotClosureSpectra()
{
   gStyle->SetOptStat(0);

   TFile *fGen  = TFile::Open("output/KtoPi-MC-Gen-Closure.root", "READ");
   TFile *fReco = TFile::Open("output/KtoPi-MC-Reco-Closure.root", "READ");
   if (!fGen || fGen->IsZombie() || !fReco || fReco->IsZombie())
   {
      Error("plotClosureSpectra", "Cannot open closure input files");
      return;
   }

   MakeSpeciesClosure(fGen, fReco, "hKPt",  "hKPtCorrected",  "Kaon p_{T} closure",   "KPtClosure_Overlay", 3.0e4, 3.0e5);
   MakeSpeciesClosure(fGen, fReco, "hPiPt", "hPiPtCorrected", "Pion p_{T} closure",   "PiPtClosure_Overlay", 1.0e5, 1.0e7);
   MakeSpeciesClosure(fGen, fReco, "hPPt",  "hPPtCorrected",  "Proton p_{T} closure", "PPtClosure_Overlay", 1.0e4, 1.0e5);

   fGen->Close();
   fReco->Close();
}
