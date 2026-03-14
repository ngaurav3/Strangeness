#!/usr/bin/env python3
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

HYBRID_ROOT = os.environ.get(
    "XSCAPE_HYBRID_ROOT",
    "result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_clean20k.root",
)
COLORLESS_ROOT = os.environ.get(
    "XSCAPE_COLORLESS_ROOT",
    "result/20260314/xscape_colorless/xscape_epem_colorless_zpole_truth_20k.root",
)
TEMPLATE_ROOT = os.environ.get(
    "XSCAPE_TEMPLATE_ROOT",
    "output/systematics_20260306_dndeta/nominal_unfold_dndeta.root",
)
OUT_DIR = os.environ.get(
    "XSCAPE_COMPARISON_OUT_DIR",
    "result/20260314/xscape_controls",
)


def clone_hist(path: str, name: str):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    h = f.Get(name)
    if h is None:
        raise RuntimeError(f"Missing {name} in {path}")
    hc = h.Clone(f"{name}_{os.path.basename(path).replace('.', '_')}")
    hc.SetDirectory(0)
    f.Close()
    return hc


def ratio_graph_from_tree(path: str, template, name: str, color: int, marker: int):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    tree = f.Get("Events")
    if tree is None:
        raise RuntimeError(f"Missing Events tree in {path}")

    h_num = template.Clone(f"hNum_{name}")
    h_den = template.Clone(f"hDen_{name}")
    h_ratio = template.Clone(f"hRatio_{name}")
    for h in (h_num, h_den, h_ratio):
        h.Reset()
        h.Sumw2()
        h.SetDirectory(0)

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        n_pi = int(getattr(tree, "nPiPt0405"))
        if n_pi <= 0:
            continue
        x = float(getattr(tree, "nChEta05"))
        h_num.Fill(x, int(getattr(tree, "nKPt0405")))
        h_den.Fill(x, n_pi)

    h_ratio.Divide(h_num, h_den, 1.0, 1.0, "")

    g = ROOT.TGraphErrors()
    j = 0
    for ib in range(1, h_ratio.GetNbinsX() + 1):
        y = h_ratio.GetBinContent(ib)
        ey = h_ratio.GetBinError(ib)
        if y == 0 and ey == 0:
            continue
        x = h_ratio.GetBinCenter(ib)
        ex = 0.45 * h_ratio.GetXaxis().GetBinWidth(ib)
        g.SetPoint(j, x, y)
        g.SetPointError(j, ex, ey)
        j += 1
    g.SetName(f"gRatio_{name}")
    g.SetMarkerStyle(marker)
    g.SetMarkerSize(1.3)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.SetLineWidth(2)
    return h_ratio, g


def dndeta_hist_from_tree(path: str, name: str, color: int, style: int):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    tree = f.Get("Events")
    if tree is None:
        raise RuntimeError(f"Missing Events tree in {path}")

    h = ROOT.TH1D(name, "", 31, -0.5, 30.5)
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        h.Fill(float(getattr(tree, "nChEta05")))
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetDirectory(0)
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(3)
    h.SetMarkerSize(0)
    f.Close()
    return h


def style_legend(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    template = clone_hist(TEMPLATE_ROOT, "hRatioDataBayes_dNdEta")

    h_hybrid_ratio, g_hybrid = ratio_graph_from_tree(
        HYBRID_ROOT, template, "hybrid", ROOT.kAzure + 7, 20
    )
    h_colorless_ratio, g_colorless = ratio_graph_from_tree(
        COLORLESS_ROOT, template, "colorless", ROOT.kRed + 1, 24
    )
    h_hybrid_dn = dndeta_hist_from_tree(
        HYBRID_ROOT, "hXscapeHybridDNdEta", ROOT.kAzure + 7, 1
    )
    h_colorless_dn = dndeta_hist_from_tree(
        COLORLESS_ROOT, "hXscapeColorlessDNdEta", ROOT.kRed + 1, 2
    )

    c = ROOT.TCanvas("cXscapeCompare", "", 1200, 600)
    c.Divide(2, 1)

    c.cd(1)
    frame1 = template.Clone("frameXscapeRatio")
    frame1.Reset()
    frame1.SetDirectory(0)
    frame1.SetTitle("")
    frame1.SetMinimum(0.05)
    frame1.SetMaximum(0.22)
    frame1.GetXaxis().SetTitle("dN_{ch}/d#eta (|#eta|<0.5)")
    frame1.GetYaxis().SetTitle("K/#pi")
    frame1.GetXaxis().SetLimits(-1.0, 31.0)
    frame1.GetXaxis().SetTitleSize(0.050)
    frame1.GetYaxis().SetTitleSize(0.050)
    frame1.GetXaxis().SetLabelSize(0.042)
    frame1.GetYaxis().SetLabelSize(0.042)
    frame1.GetYaxis().SetTitleOffset(1.10)
    frame1.Draw("AXIS")
    g_hybrid.Draw("P SAME")
    g_colorless.Draw("P SAME")
    leg1 = ROOT.TLegend(0.20, 0.75, 0.58, 0.88)
    style_legend(leg1)
    leg1.SetTextSize(0.034)
    leg1.AddEntry(g_hybrid, "X-SCAPE hybrid hadronization", "p")
    leg1.AddEntry(g_colorless, "X-SCAPE colorless hadronization", "p")
    leg1.Draw()

    c.cd(2)
    frame2 = h_hybrid_dn.Clone("frameXscapeDNdEta")
    frame2.Reset()
    frame2.SetDirectory(0)
    frame2.SetTitle("")
    frame2.SetMinimum(0.0)
    frame2.SetMaximum(1.35 * max(h_hybrid_dn.GetMaximum(), h_colorless_dn.GetMaximum()))
    frame2.GetXaxis().SetTitle("dN_{ch}/d#eta (|#eta|<0.5)")
    frame2.GetYaxis().SetTitle("Normalized events")
    frame2.GetXaxis().SetTitleSize(0.050)
    frame2.GetYaxis().SetTitleSize(0.050)
    frame2.GetXaxis().SetLabelSize(0.042)
    frame2.GetYaxis().SetLabelSize(0.042)
    frame2.GetYaxis().SetTitleOffset(1.10)
    frame2.Draw("HIST")
    h_hybrid_dn.Draw("HIST SAME")
    h_colorless_dn.Draw("HIST SAME")
    leg2 = ROOT.TLegend(0.20, 0.75, 0.58, 0.88)
    style_legend(leg2)
    leg2.SetTextSize(0.034)
    leg2.AddEntry(h_hybrid_dn, "X-SCAPE hybrid hadronization", "l")
    leg2.AddEntry(h_colorless_dn, "X-SCAPE colorless hadronization", "l")
    leg2.Draw()

    pdf = os.path.join(OUT_DIR, "XscapeHadronization_Comparison.pdf")
    png = os.path.join(OUT_DIR, "XscapeHadronization_Comparison.png")
    root_out = os.path.join(OUT_DIR, "XscapeHadronization_Comparison.root")
    txt_out = os.path.join(OUT_DIR, "XscapeHadronization_Comparison.txt")

    c.SaveAs(pdf)
    c.SaveAs(png)

    fout = ROOT.TFile(root_out, "RECREATE")
    h_hybrid_ratio.Write()
    h_colorless_ratio.Write()
    g_hybrid.Write()
    g_colorless.Write()
    h_hybrid_dn.Write()
    h_colorless_dn.Write()
    fout.Close()

    with open(txt_out, "w", encoding="ascii") as fout_txt:
        fout_txt.write("bin_center,hybrid_kpi,colorless_kpi,hybrid_over_colorless\n")
        for ib in range(1, h_hybrid_ratio.GetNbinsX() + 1):
            hy = h_hybrid_ratio.GetBinContent(ib)
            co = h_colorless_ratio.GetBinContent(ib)
            ratio = hy / co if co > 0 else -1.0
            fout_txt.write(
                f"{h_hybrid_ratio.GetBinCenter(ib):.6f},{hy:.6f},{co:.6f},{ratio:.6f}\n"
            )

    print("Wrote:", pdf)
    print("Wrote:", root_out)
    print("Wrote:", txt_out)


if __name__ == "__main__":
    main()
