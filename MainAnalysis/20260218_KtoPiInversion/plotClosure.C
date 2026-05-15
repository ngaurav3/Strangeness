#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPad.h"
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
}

void plotClosure()
{
    gStyle->SetOptStat(0);

    TFile *fGen  = TFile::Open("output/KtoPi-MC-Gen-Closure.root", "READ");
    TFile *fReco = TFile::Open("output/KtoPi-MC-Reco-Closure.root", "READ");
    if (!fGen || fGen->IsZombie() || !fReco || fReco->IsZombie())
    {
        Error("plotClosure", "Cannot open closure input files");
        return;
    }

    TH1D *hGen  = dynamic_cast<TH1D*>(fGen->Get("hKoverPi"));
    TH1D *hReco = dynamic_cast<TH1D*>(fReco->Get("hKoverPiCorrected"));
    if (!hGen || !hReco)
    {
        Error("plotClosure", "Missing hKoverPi or hKoverPiCorrected");
        return;
    }

    TH1D *hGenDraw  = (TH1D*)hGen->Clone("hGenDraw");
    TH1D *hRecoDraw = (TH1D*)hReco->Clone("hRecoDraw");
    hGenDraw->SetDirectory(0);
    hRecoDraw->SetDirectory(0);
    fGen->Close();
    fReco->Close();

    TH1D *hRatio = (TH1D*)hRecoDraw->Clone("hClosureRatio");
    hRatio->Divide(hGenDraw);
    TH1D *hRecoResidualCorr = (TH1D*)hRecoDraw->Clone("hRecoResidualCorr");
    hRecoResidualCorr->Divide(hRatio);
    TH1D *hRatioResidualCorr = (TH1D*)hRecoResidualCorr->Clone("hRatioResidualCorr");
    hRatioResidualCorr->Divide(hGenDraw);

    hRecoDraw->SetMarkerStyle(20);
    hRecoDraw->SetMarkerSize(1.5);
    hRecoDraw->SetMarkerColor(kRed + 1);
    hRecoDraw->SetLineColor(kRed + 1);
    hGenDraw->SetLineColor(kBlue + 1);
    hGenDraw->SetLineStyle(2);
    hGenDraw->SetLineWidth(2);
    hRecoResidualCorr->SetMarkerStyle(20);
    hRecoResidualCorr->SetMarkerSize(1.5);
    hRecoResidualCorr->SetMarkerColor(kBlue + 1);
    hRecoResidualCorr->SetLineColor(kBlue + 1);

    TCanvas *c = new TCanvas("cClosure", "MC closure K/pi", 820, 800);
    TPad *p1 = new TPad("p1", "", 0, 0.30, 1, 1);
    TPad *p2 = new TPad("p2", "", 0, 0.00, 1, 0.30);
    ConfigureUpperPad(p1);
    ConfigureLowerPad(p2);
    p1->Draw();
    p2->Draw();

    p1->cd();
    hRecoDraw->SetTitle("MC Closure: corrected reco vs generator;N_{ch}^{tag};K/#pi");
    hRecoDraw->SetMinimum(0.0);
    const double yMax = 1.25 * std::max(hRecoDraw->GetMaximum(), hGenDraw->GetMaximum());
    hRecoDraw->SetMaximum(yMax > 0.0 ? yMax : 1.0);
    ConfigureUpperAxis(hRecoDraw);
    hRecoDraw->Draw("E1");
    hGenDraw->Draw("HIST SAME");
    hRecoResidualCorr->Draw("E1 SAME");

    TLegend *leg = new TLegend(0.33, 0.66, 0.83, 0.87);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.048);
    leg->AddEntry(hRecoDraw, "MC reco corrected (closure mode)", "lep");
    leg->AddEntry(hRecoResidualCorr, "MC reco after residual correction", "lep");
    leg->AddEntry(hGenDraw, "MC generator", "l");
    leg->Draw();

    p2->cd();
    hRatio->SetTitle(";N_{ch}^{tag};Closure ratio");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.3);
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    ConfigureRatioAxis(hRatio);
    hRatio->Draw("E1");
    TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
    hRatio->Draw("E1 SAME");
    hRatioResidualCorr->SetMarkerStyle(20);
    hRatioResidualCorr->SetMarkerSize(1.3);
    hRatioResidualCorr->SetMarkerColor(kBlue + 1);
    hRatioResidualCorr->SetLineColor(kBlue + 1);
    hRatioResidualCorr->Draw("E1 SAME");

    c->SaveAs("KtoPiClosure_Overlay.png");
    c->SaveAs("KtoPiClosure_Overlay.pdf");
}
