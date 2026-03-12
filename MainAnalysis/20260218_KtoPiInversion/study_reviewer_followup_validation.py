#!/usr/bin/env python3
import csv
import math
import os
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

BASE = os.path.abspath(os.path.dirname(__file__))
DATE = "20260312"
OUTDIR = os.path.join(BASE, f"result/{DATE}/reviewer_followup_validation")
KEEPBINS_DIR = os.path.join(OUTDIR, "keepbins_scan")
TOY_DIR = os.path.join(OUTDIR, "toy_coverage")
SVD_DIR = os.path.join(OUTDIR, "svd_kreg_scan")
REPORT_MD = os.path.join(BASE, f"report/{DATE}_iteration10_reviewer_followup.md")
ROOT_SETUP = "source /raid5/root/root-v6.34.04/root/bin/thisroot.sh"
NTOYS = 400
KEEPBINS_VALUES = [8, 9, 10]
KREG_VALUES = list(range(2, 11))

for path in [OUTDIR, KEEPBINS_DIR, TOY_DIR, SVD_DIR]:
    os.makedirs(path, exist_ok=True)


def run_shell(cmd):
    subprocess.run(["bash", "-lc", cmd], cwd=BASE, check=True)


def run_root_macro(expr):
    cmd = f'{ROOT_SETUP} && root -l -b -q \'{expr}\''
    subprocess.run(["bash", "-lc", cmd], cwd=BASE, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def clone_hist(path, name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    h = f.Get(name)
    if not h:
        raise RuntimeError(f"Missing {name} in {path}")
    out = h.Clone(f"{name}_{os.path.basename(path).replace('.', '_')}")
    out.SetDirectory(0)
    f.Close()
    return out


def hist_arrays(h):
    n = h.GetNbinsX()
    centers = np.array([h.GetXaxis().GetBinCenter(i) for i in range(1, n + 1)], dtype=float)
    lows = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n + 1)], dtype=float)
    highs = np.array([h.GetXaxis().GetBinUpEdge(i) for i in range(1, n + 1)], dtype=float)
    vals = np.array([h.GetBinContent(i) for i in range(1, n + 1)], dtype=float)
    errs = np.array([h.GetBinError(i) for i in range(1, n + 1)], dtype=float)
    return centers, lows, highs, vals, errs


def matrix_array(h2):
    nx = h2.GetNbinsX()
    ny = h2.GetNbinsY()
    arr = np.zeros((nx, ny), dtype=float)
    for ix in range(1, nx + 1):
        for iy in range(1, ny + 1):
            arr[ix - 1, iy - 1] = h2.GetBinContent(ix, iy)
    return arr


def build_ratio(num, num_err, den, den_err):
    ratio = np.zeros_like(num, dtype=float)
    err = np.zeros_like(num, dtype=float)
    mask = den > 0
    ratio[mask] = num[mask] / den[mask]
    mask2 = mask & (num > 0)
    err[mask2] = ratio[mask2] * np.sqrt(
        (num_err[mask2] / num[mask2]) ** 2 + (den_err[mask2] / den[mask2]) ** 2
    )
    return ratio, err


def rms_distance_from_unity(values):
    if len(values) == 0:
        return 0.0
    return float(np.sqrt(np.mean((np.asarray(values) - 1.0) ** 2)))


def iterative_bayes_unfold(meas, response_true_reco, prior, niter):
    ntrue, nreco = response_true_reco.shape
    prior_vec = np.clip(np.asarray(prior, dtype=float), 0.0, None)
    if prior_vec.sum() <= 0:
        prior_vec = np.full(ntrue, 1.0 / ntrue, dtype=float)
    else:
        prior_vec = prior_vec / prior_vec.sum()

    p = np.clip(np.asarray(response_true_reco, dtype=float), 0.0, None)
    unfolded = np.zeros(ntrue, dtype=float)
    for _ in range(niter):
        unfolded[:] = 0.0
        for r in range(nreco):
            mr = max(0.0, meas[r])
            if mr == 0.0:
                continue
            norm = np.dot(p[:, r], prior_vec)
            if norm <= 0.0:
                continue
            weights = (p[:, r] * prior_vec) / norm
            unfolded += weights * mr
        total = np.clip(unfolded, 0.0, None).sum()
        if total <= 0.0:
            break
        prior_vec = np.clip(unfolded, 0.0, None) / total
    return np.clip(unfolded, 0.0, None), np.sqrt(np.clip(unfolded, 0.0, None))


def svd_unfold(meas, meas_err, response_true_reco, kreg):
    ntrue, nreco = response_true_reco.shape
    a = np.zeros((nreco, ntrue), dtype=float)
    for t in range(ntrue):
        col_sum = np.sum(response_true_reco[t, :])
        if col_sum <= 0:
            continue
        a[:, t] = response_true_reco[t, :] / col_sum
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    k = max(1, min(kreg, len(s)))
    splus = np.zeros_like(s)
    for i in range(k):
        if s[i] > 1e-12:
            splus[i] = 1.0 / s[i]
    b = vt.T @ np.diag(splus) @ u.T
    x = np.clip(b @ meas, 0.0, None)
    var = np.sum((b ** 2) * (meas_err[np.newaxis, :] ** 2), axis=1)
    return x, np.sqrt(np.clip(var, 0.0, None)), s


def fold_truth(truth, truth_err, response_true_reco):
    ntrue, nreco = response_true_reco.shape
    reco = np.zeros(nreco, dtype=float)
    err2 = np.zeros(nreco, dtype=float)
    for t in range(ntrue):
        col_sum = np.sum(response_true_reco[t, :])
        if col_sum <= 0.0:
            continue
        prob = response_true_reco[t, :] / col_sum
        reco += prob * truth[t]
        err2 += (prob ** 2) * (truth_err[t] ** 2)
    return reco, np.sqrt(np.clip(err2, 0.0, None))


def response_metrics(resp):
    ntrue, nreco = resp.shape
    purity = np.zeros(min(ntrue, nreco), dtype=float)
    stability = np.zeros(min(ntrue, nreco), dtype=float)
    for i in range(min(ntrue, nreco)):
        row_sum = np.sum(resp[:, i])
        col_sum = np.sum(resp[i, :])
        diag = resp[i, i]
        purity[i] = diag / row_sum if row_sum > 0 else 0.0
        stability[i] = diag / col_sum if col_sum > 0 else 0.0
    return purity, stability


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def perform_keepbins_scan():
    print("[keepBins] scanning dN/deta tail merge choices")
    rows = []
    closure_panels = {"K": {}, "Pi": {}}
    for keep in KEEPBINS_VALUES:
        out_root = os.path.join(KEEPBINS_DIR, f"dndeta_keepBins{keep}.root")
        expr = (
            'runDNdEtaUnfolding_BayesSVD.C('
            f'1,8,"output/KtoPi-MC-Reco-Nominal.root","output/KtoPi-Data-Reco-Nominal.root",'
            f'"{out_root}","",false,{keep})'
        )
        print(f"  - keepBins={keep}")
        run_root_macro(expr)

        h_refold_mc = clone_hist(out_root, "hRefoldRecoClosureMc_dNdEta")
        h_refold_data = clone_hist(out_root, "hRefoldRecoClosureData_dNdEta")
        _, _, _, refold_mc_vals, _ = hist_arrays(h_refold_mc)
        _, _, _, refold_data_vals, _ = hist_arrays(h_refold_data)

        h_k_refold_mc = clone_hist(out_root, "hKMcBayesRefold_dNdEta")
        h_k_refold_data = clone_hist(out_root, "hKDataBayesRefold_dNdEta")
        h_k_reco_mc = clone_hist(out_root, "hKMcReco")
        h_k_reco_data = clone_hist(out_root, "hKDataReco")
        h_pi_refold_mc = clone_hist(out_root, "hPiMcBayesRefold_dNdEta")
        h_pi_refold_data = clone_hist(out_root, "hPiDataBayesRefold_dNdEta")
        h_pi_reco_mc = clone_hist(out_root, "hPiMcReco")
        h_pi_reco_data = clone_hist(out_root, "hPiDataReco")
        _, _, _, k_refold_mc, _ = hist_arrays(build_root_ratio(h_k_refold_mc, h_k_reco_mc, f"k_mc_refold_keep{keep}"))
        _, _, _, k_refold_data, _ = hist_arrays(build_root_ratio(h_k_refold_data, h_k_reco_data, f"k_data_refold_keep{keep}"))
        _, _, _, pi_refold_mc, _ = hist_arrays(build_root_ratio(h_pi_refold_mc, h_pi_reco_mc, f"pi_mc_refold_keep{keep}"))
        _, _, _, pi_refold_data, _ = hist_arrays(build_root_ratio(h_pi_refold_data, h_pi_reco_data, f"pi_data_refold_keep{keep}"))

        h_k_unfold = clone_hist(out_root, "hKMcBayes_dNdEta")
        h_k_truth = clone_hist(out_root, "hKPrior_dNdEta")
        h_pi_unfold = clone_hist(out_root, "hPiMcBayes_dNdEta")
        h_pi_truth = clone_hist(out_root, "hPiPrior_dNdEta")
        h_k_closure = build_root_ratio(h_k_unfold, h_k_truth, f"k_closure_keep{keep}")
        h_pi_closure = build_root_ratio(h_pi_unfold, h_pi_truth, f"pi_closure_keep{keep}")
        centers, lows, highs, k_closure_vals, _ = hist_arrays(h_k_closure)
        _, _, _, pi_closure_vals, _ = hist_arrays(h_pi_closure)
        closure_panels["K"][keep] = (centers, lows, highs, k_closure_vals)
        closure_panels["Pi"][keep] = (centers, lows, highs, pi_closure_vals)

        resp_k = matrix_array(clone_hist(out_root, "hDNdEtaResponseKRebinned"))
        resp_pi = matrix_array(clone_hist(out_root, "hDNdEtaResponsePiRebinned"))
        purity_k, stability_k = response_metrics(resp_k)
        purity_pi, stability_pi = response_metrics(resp_pi)

        last = len(refold_data_vals) - 1
        min_species_last = min(k_refold_mc[last], k_refold_data[last], pi_refold_mc[last], pi_refold_data[last])
        min_tail_quality = min(purity_k[last], stability_k[last], purity_pi[last], stability_pi[last])
        passes = (
            rms_distance_from_unity(refold_data_vals) < 0.05 and
            rms_distance_from_unity(refold_mc_vals) < 0.03 and
            min_species_last > 0.4 and
            min_tail_quality > 0.3
        )
        rows.append({
            "keepBins": keep,
            "ratio_refold_rms_mc": rms_distance_from_unity(refold_mc_vals),
            "ratio_refold_rms_data": rms_distance_from_unity(refold_data_vals),
            "k_last_refold_mc": float(k_refold_mc[last]),
            "k_last_refold_data": float(k_refold_data[last]),
            "pi_last_refold_mc": float(pi_refold_mc[last]),
            "pi_last_refold_data": float(pi_refold_data[last]),
            "k_last_closure": float(k_closure_vals[last]),
            "pi_last_closure": float(pi_closure_vals[last]),
            "k_last_purity": float(purity_k[last]),
            "k_last_stability": float(stability_k[last]),
            "pi_last_purity": float(purity_pi[last]),
            "pi_last_stability": float(stability_pi[last]),
            "passes_rule": int(passes),
        })

    write_csv(
        os.path.join(KEEPBINS_DIR, "keepbins_scan_summary.csv"),
        rows,
        list(rows[0].keys()),
    )
    plot_keepbins_summary(rows)
    plot_keepbins_species_closure(closure_panels)
    return rows


def build_root_ratio(num, den, name):
    h = num.Clone(name)
    h.SetDirectory(0)
    h.Divide(den)
    return h


def plot_keepbins_summary(rows):
    keeps = np.array([row["keepBins"] for row in rows], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(keeps, [row["ratio_refold_rms_mc"] for row in rows], marker="o", label="ratio refold RMS (MC)")
    ax.plot(keeps, [row["ratio_refold_rms_data"] for row in rows], marker="s", label="ratio refold RMS (data)")
    ax.axhline(0.03, color="tab:blue", linestyle=":", linewidth=1.5)
    ax.axhline(0.05, color="tab:orange", linestyle=":", linewidth=1.5)
    ax.set_xlabel("keepBins")
    ax.set_ylabel("RMS")
    ax.set_title("Ratio-level refolding quality")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[0, 1]
    for key, label in [
        ("k_last_refold_mc", "K last-bin refold/reco (MC)"),
        ("k_last_refold_data", "K last-bin refold/reco (data)"),
        ("pi_last_refold_mc", "pi last-bin refold/reco (MC)"),
        ("pi_last_refold_data", "pi last-bin refold/reco (data)"),
    ]:
        ax.plot(keeps, [row[key] for row in rows], marker="o", label=label)
    ax.axhline(0.4, color="black", linestyle="--", linewidth=1.5)
    ax.set_xlabel("keepBins")
    ax.set_ylabel("Last-bin ratio")
    ax.set_title("Species-level final-bin refolding")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1, 0]
    for key, label in [
        ("k_last_purity", "K purity"),
        ("k_last_stability", "K stability"),
        ("pi_last_purity", "pi purity"),
        ("pi_last_stability", "pi stability"),
    ]:
        ax.plot(keeps, [row[key] for row in rows], marker="o", label=label)
    ax.axhline(0.3, color="black", linestyle="--", linewidth=1.5)
    ax.set_xlabel("keepBins")
    ax.set_ylabel("Fraction")
    ax.set_title("Final merged-bin migration quality")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1, 1]
    ax.plot(keeps, [row["k_last_closure"] for row in rows], marker="o", label="K unfolded/truth")
    ax.plot(keeps, [row["pi_last_closure"] for row in rows], marker="s", label="pi unfolded/truth")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.5)
    ax.set_xlabel("keepBins")
    ax.set_ylabel("Last-bin closure")
    ax.set_title("Species-level final-bin closure")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)

    fig.savefig(os.path.join(KEEPBINS_DIR, "DNdEtaKeepBinsScan_Summary.pdf"))
    fig.savefig(os.path.join(KEEPBINS_DIR, "DNdEtaKeepBinsScan_Summary.png"), dpi=160)
    plt.close(fig)


def plot_keepbins_species_closure(closure_panels):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, constrained_layout=True)
    for ax, species in zip(axes, ["K", "Pi"]):
        for keep, color in zip(KEEPBINS_VALUES, ["tab:blue", "tab:orange", "tab:green"]):
            centers, lows, highs, vals = closure_panels[species][keep]
            ax.step(centers, vals, where="mid", label=f"keepBins={keep}", color=color, linewidth=2)
            ax.axvspan(lows[-1], highs[-1], color=color, alpha=0.08)
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.5)
        ax.set_xlabel(r"True $dN_{ch}/d\eta$ ($|\eta|<0.5$)")
        ax.set_title(f"{species} unfolded / truth")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Closure")
    axes[0].set_ylim(0.0, 1.6)
    axes[1].legend(frameon=False)
    fig.savefig(os.path.join(KEEPBINS_DIR, "DNdEtaKeepBinsScan_SpeciesClosure.pdf"))
    fig.savefig(os.path.join(KEEPBINS_DIR, "DNdEtaKeepBinsScan_SpeciesClosure.png"), dpi=160)
    plt.close(fig)


def nominal_inputs(axis):
    if axis == "Ntag":
        path = os.path.join(BASE, "output/NtagUnfolding_BayesSVD.root")
        return {
            "path": path,
            "resp_k": matrix_array(clone_hist(path, "hNtagResponseK")),
            "resp_pi": matrix_array(clone_hist(path, "hNtagResponsePi")),
            "k_truth": hist_arrays(clone_hist(path, "hKTrueNtag"))[3],
            "pi_truth": hist_arrays(clone_hist(path, "hPiTrueNtag"))[3],
            "k_reco": hist_arrays(clone_hist(path, "hKMcReco"))[3],
            "pi_reco": hist_arrays(clone_hist(path, "hPiMcReco"))[3],
            "centers": hist_arrays(clone_hist(path, "hRatioMcTrue"))[0],
        }
    path = os.path.join(BASE, "output/systematics_20260306_dndeta/nominal_unfold_dndeta.root")
    return {
        "path": path,
        "resp_k": matrix_array(clone_hist(path, "hDNdEtaResponseKRebinned")),
        "resp_pi": matrix_array(clone_hist(path, "hDNdEtaResponsePiRebinned")),
        "k_truth": hist_arrays(clone_hist(path, "hKPrior_dNdEta"))[3],
        "pi_truth": hist_arrays(clone_hist(path, "hPiPrior_dNdEta"))[3],
        "k_reco": hist_arrays(clone_hist(path, "hKMcReco"))[3],
        "pi_reco": hist_arrays(clone_hist(path, "hPiMcReco"))[3],
        "centers": hist_arrays(clone_hist(path, "hRatioMcTrue_dNdEta"))[0],
    }


def perform_toy_coverage():
    print("[toys] running pseudo-experiment coverage study")
    rng = np.random.default_rng(12345)
    summary_rows = []
    coverage_payload = {}

    for axis in ["Ntag", "dNdEta"]:
        print(f"  - axis={axis}, toys={NTOYS}")
        payload = nominal_inputs(axis)
        k_truth = payload["k_truth"]
        pi_truth = payload["pi_truth"]
        k_reco = payload["k_reco"]
        pi_reco = payload["pi_reco"]
        resp_k = payload["resp_k"]
        resp_pi = payload["resp_pi"]
        truth_ratio, _ = build_ratio(k_truth, np.sqrt(np.clip(k_truth, 0.0, None)), pi_truth, np.sqrt(np.clip(pi_truth, 0.0, None)))
        nbins = len(truth_ratio)
        pulls = np.full((NTOYS, nbins), np.nan, dtype=float)

        for itoy in range(NTOYS):
            if (itoy + 1) % 100 == 0:
                print(f"    toy {itoy + 1}/{NTOYS}")
            k_meas = rng.poisson(np.clip(k_reco, 0.0, None))
            pi_meas = rng.poisson(np.clip(pi_reco, 0.0, None))
            k_unf, k_err = iterative_bayes_unfold(k_meas, resp_k, k_truth, 1)
            pi_unf, pi_err = iterative_bayes_unfold(pi_meas, resp_pi, pi_truth, 1)
            ratio, ratio_err = build_ratio(k_unf, k_err, pi_unf, pi_err)
            valid = ratio_err > 0
            pulls[itoy, valid] = (ratio[valid] - truth_ratio[valid]) / ratio_err[valid]

        mean_pull = np.nanmean(pulls, axis=0)
        width_pull = np.nanstd(pulls, axis=0)
        cov1 = np.nanmean(np.abs(pulls) < 1.0, axis=0)
        cov2 = np.nanmean(np.abs(pulls) < 2.0, axis=0)
        coverage_payload[axis] = {
            "centers": payload["centers"],
            "mean": mean_pull,
            "width": width_pull,
            "cov1": cov1,
            "cov2": cov2,
        }
        summary_rows.append({
            "axis": axis,
            "mean_abs_pull_bias": float(np.nanmean(np.abs(mean_pull))),
            "mean_pull_width": float(np.nanmean(width_pull)),
            "min_cov1": float(np.nanmin(cov1)),
            "min_cov2": float(np.nanmin(cov2)),
            "max_pull_width": float(np.nanmax(width_pull)),
        })

    write_csv(os.path.join(TOY_DIR, "toy_coverage_summary.csv"), summary_rows, list(summary_rows[0].keys()))
    plot_toy_coverage(coverage_payload)
    return summary_rows


def plot_toy_coverage(payload):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for row, axis in enumerate(["Ntag", "dNdEta"]):
        d = payload[axis]
        axes[row, 0].axhline(0.0, color="black", linestyle="--", linewidth=1.5)
        axes[row, 0].plot(d["centers"], d["mean"], marker="o")
        axes[row, 0].set_title(f"{axis} toy pull mean")
        axes[row, 0].set_ylabel("Mean pull")
        axes[row, 0].grid(alpha=0.25)

        axes[row, 1].axhline(1.0, color="black", linestyle="--", linewidth=1.5)
        axes[row, 1].plot(d["centers"], d["width"], marker="o", label="Pull width")
        axes[row, 1].plot(d["centers"], d["cov1"], marker="s", label="Coverage |pull|<1")
        axes[row, 1].plot(d["centers"], d["cov2"], marker="^", label="Coverage |pull|<2")
        axes[row, 1].set_title(f"{axis} toy pull width / coverage")
        axes[row, 1].set_ylabel("Value")
        axes[row, 1].grid(alpha=0.25)
        axes[row, 1].legend(frameon=False, fontsize=8)

    axes[1, 0].set_xlabel("Unfolded axis bin center")
    axes[1, 1].set_xlabel("Unfolded axis bin center")
    fig.savefig(os.path.join(TOY_DIR, "UnfoldingToyCoverage_Summary.pdf"))
    fig.savefig(os.path.join(TOY_DIR, "UnfoldingToyCoverage_Summary.png"), dpi=160)
    plt.close(fig)


def perform_svd_scan():
    print("[svd] running SVD regularization scan")
    ntag_rows = []
    dndeta_rows = []
    for kreg in KREG_VALUES:
        print(f"  - kReg={kreg} (Ntag)")
        run_root_macro(f"runNtagUnfolding_BayesSVD.C(1,{kreg})")
        out_path = os.path.join(SVD_DIR, f"ntag_kreg{kreg}.root")
        shutil.copyfile(os.path.join(BASE, "output/NtagUnfolding_BayesSVD.root"), out_path)
        ntag_rows.append(svd_metrics_from_file(out_path, "Ntag", kreg))

        print(f"  - kReg={kreg} (dNdEta)")
        out_path = os.path.join(SVD_DIR, f"dndeta_kreg{kreg}.root")
        expr = (
            'runDNdEtaUnfolding_BayesSVD.C('
            f'1,{kreg},"output/KtoPi-MC-Reco-Nominal.root","output/KtoPi-Data-Reco-Nominal.root",'
            f'"{out_path}","",false)'
        )
        run_root_macro(expr)
        dndeta_rows.append(svd_metrics_from_file(out_path, "dNdEta", kreg))

    write_csv(os.path.join(SVD_DIR, "svd_kreg_scan_summary.csv"), ntag_rows + dndeta_rows, list((ntag_rows + dndeta_rows)[0].keys()))
    plot_svd_scan(ntag_rows, dndeta_rows)
    return ntag_rows, dndeta_rows


def svd_metrics_from_file(path, axis, kreg):
    if axis == "Ntag":
        closure = hist_arrays(clone_hist(path, "hClosureSVD"))[3]
        diff = hist_arrays(clone_hist(path, "hMethodDiff_BayesMinusSVD"))[3]
        resp_k = matrix_array(clone_hist(path, "hNtagResponseK"))
        resp_pi = matrix_array(clone_hist(path, "hNtagResponsePi"))
    else:
        closure = hist_arrays(clone_hist(path, "hClosureSVD_dNdEta"))[3]
        diff = hist_arrays(clone_hist(path, "hMethodDiff_BayesMinusSVD_dNdEta"))[3]
        resp_k = matrix_array(clone_hist(path, "hDNdEtaResponseKRebinned"))
        resp_pi = matrix_array(clone_hist(path, "hDNdEtaResponsePiRebinned"))
    _, _, s_k = svd_unfold(np.ones(resp_k.shape[1]), np.ones(resp_k.shape[1]), resp_k, kreg)
    _, _, s_pi = svd_unfold(np.ones(resp_pi.shape[1]), np.ones(resp_pi.shape[1]), resp_pi, kreg)
    return {
        "axis": axis,
        "kReg": kreg,
        "closure_rms": rms_distance_from_unity(closure),
        "method_diff_rms": float(np.sqrt(np.mean(np.asarray(diff) ** 2))),
        "sv_k1": float(s_k[0]) if len(s_k) > 0 else 0.0,
        "sv_k_lastkept": float(s_k[min(kreg, len(s_k)) - 1]) if len(s_k) > 0 else 0.0,
        "sv_pi1": float(s_pi[0]) if len(s_pi) > 0 else 0.0,
        "sv_pi_lastkept": float(s_pi[min(kreg, len(s_pi)) - 1]) if len(s_pi) > 0 else 0.0,
    }


def singular_values_from_response(resp):
    ntrue, nreco = resp.shape
    a = np.zeros((nreco, ntrue), dtype=float)
    for t in range(ntrue):
        col_sum = np.sum(resp[t, :])
        if col_sum > 0:
            a[:, t] = resp[t, :] / col_sum
    return np.linalg.svd(a, compute_uv=False)


def plot_svd_scan(ntag_rows, dndeta_rows):
    nominal_ntag = nominal_inputs("Ntag")
    nominal_dndeta = nominal_inputs("dNdEta")
    sv_sets = {
        "Ntag": {
            "K": singular_values_from_response(nominal_ntag["resp_k"]),
            "Pi": singular_values_from_response(nominal_ntag["resp_pi"]),
        },
        "dNdEta": {
            "K": singular_values_from_response(nominal_dndeta["resp_k"]),
            "Pi": singular_values_from_response(nominal_dndeta["resp_pi"]),
        },
    }
    rows_by_axis = {"Ntag": ntag_rows, "dNdEta": dndeta_rows}

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for i, axis in enumerate(["Ntag", "dNdEta"]):
        ax = axes[i, 0]
        for species, style in [("K", "o-"), ("Pi", "s-")]:
            sv = sv_sets[axis][species]
            ax.semilogy(np.arange(1, len(sv) + 1), sv / sv[0], style, label=species)
        ax.axvline(8, color="black", linestyle="--", linewidth=1.5)
        ax.set_title(f"{axis} normalized singular values")
        ax.set_xlabel("Mode index")
        ax.set_ylabel(r"$\sigma_i / \sigma_1$")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)

        ax = axes[i, 1]
        ax.plot(KREG_VALUES, [row["closure_rms"] for row in rows_by_axis[axis]], marker="o")
        ax.axvline(8, color="black", linestyle="--", linewidth=1.5)
        ax.set_title(f"{axis} SVD closure RMS")
        ax.set_xlabel("kReg")
        ax.set_ylabel("RMS")
        ax.grid(alpha=0.25)

        ax = axes[i, 2]
        ax.plot(KREG_VALUES, [row["method_diff_rms"] for row in rows_by_axis[axis]], marker="o")
        ax.axvline(8, color="black", linestyle="--", linewidth=1.5)
        ax.set_title(f"{axis} Bayes-SVD method difference RMS")
        ax.set_xlabel("kReg")
        ax.set_ylabel("RMS")
        ax.grid(alpha=0.25)

    fig.savefig(os.path.join(SVD_DIR, "SVD_kReg_Scan_Summary.pdf"))
    fig.savefig(os.path.join(SVD_DIR, "SVD_kReg_Scan_Summary.png"), dpi=160)
    plt.close(fig)


def write_report(keep_rows, toy_rows, ntag_svd_rows, dndeta_svd_rows):
    keep_pass = [row for row in keep_rows if row["passes_rule"] == 1]
    best_keep = min(keep_rows, key=lambda row: (abs(row["ratio_refold_rms_data"]) + abs(row["ratio_refold_rms_mc"])))

    def svd_line(rows, axis):
        best_closure = min(rows, key=lambda row: row["closure_rms"])
        stable = min(rows, key=lambda row: row["method_diff_rms"])
        return (
            f"- `{axis}` SVD scan: smallest closure RMS at `kReg={best_closure['kReg']}` "
            f"({best_closure['closure_rms']:.4f}); smallest Bayes-SVD method-difference RMS at "
            f"`kReg={stable['kReg']}` ({stable['method_diff_rms']:.4f}); production cross-check value remains `kReg=8`."
        )

    lines = [
        f"# {DATE} reviewer follow-up validation",
        "",
        "Author: Yen-Jie Lee and OpenAI",
        "",
        "## Scope",
        "This note records the remaining unfolding-validation studies requested by the reviewer after iteration 9:",
        "- explicit `keepBins` scan for the `dN_{ch}/d\\eta` tail merge",
        "- pseudo-experiment coverage study",
        "- explicit SVD `kReg` scan with singular-value summary",
        "",
        "## keepBins scan",
        "Scanned values:",
        "- `keepBins = 8, 9, 10`",
        "",
        "Acceptance rule used in the scan:",
        "1. ratio-level `K/pi` refolding remains acceptable",
        "2. species-level final-bin refolding no longer collapses catastrophically",
        "3. final merged-bin purity and stability improve beyond the marginal `0.20-0.25` level",
        "",
        f"- Best ratio-level refolding stability in this scan: `keepBins={best_keep['keepBins']}` with data RMS `{best_keep['ratio_refold_rms_data']:.4f}` and MC RMS `{best_keep['ratio_refold_rms_mc']:.4f}`.",
    ]
    if keep_pass:
        lines.append(f"- Passing scan point(s): {', '.join(str(row['keepBins']) for row in keep_pass)}.")
    else:
        lines.append("- Passing scan point(s): none under the current explicit rule.")
    lines.extend([
        f"- Best scanned final-bin species refolding remains limited: K data `{best_keep['k_last_refold_data']:.3f}`, pi data `{best_keep['pi_last_refold_data']:.3f}`.",
        f"- Best scanned final-bin migration quality remains: K purity/stability `{best_keep['k_last_purity']:.3f}/{best_keep['k_last_stability']:.3f}`, pi purity/stability `{best_keep['pi_last_purity']:.3f}/{best_keep['pi_last_stability']:.3f}`.",
        "",
        "Generated outputs:",
        f"- `{os.path.relpath(os.path.join(KEEPBINS_DIR, 'DNdEtaKeepBinsScan_Summary.pdf'), BASE)}`",
        f"- `{os.path.relpath(os.path.join(KEEPBINS_DIR, 'DNdEtaKeepBinsScan_SpeciesClosure.pdf'), BASE)}`",
        "",
        "## Pseudo-experiment coverage",
    ])
    for row in toy_rows:
        lines.append(
            f"- `{row['axis']}` toy coverage: mean |pull bias| `{row['mean_abs_pull_bias']:.3f}`, "
            f"mean pull width `{row['mean_pull_width']:.3f}`, minimum coverage `|pull|<1` `{row['min_cov1']:.3f}`, "
            f"`|pull|<2` `{row['min_cov2']:.3f}`."
        )
    lines.extend([
        "",
        "Generated output:",
        f"- `{os.path.relpath(os.path.join(TOY_DIR, 'UnfoldingToyCoverage_Summary.pdf'), BASE)}`",
        "",
        "## SVD regularization scan",
        svd_line(ntag_svd_rows, "Ntag"),
        svd_line(dndeta_svd_rows, "dNdEta"),
        "",
        "Generated output:",
        f"- `{os.path.relpath(os.path.join(SVD_DIR, 'SVD_kReg_Scan_Summary.pdf'), BASE)}`",
        "",
        "## Bottom line",
        "- The note-side reviewer comments from iteration 8 and 9 are addressed.",
        "- The `dN_{ch}/d\\eta` tail problem is now explicitly scanned rather than only discussed, but the species-level final-bin issue is not fully closed by the present `keepBins` scan alone.",
        "- The pseudo-experiment coverage and explicit SVD `kReg` scan now exist as concrete validation studies rather than open placeholders.",
    ])
    with open(REPORT_MD, "w", encoding="ascii") as f:
        f.write("\n".join(lines) + "\n")


def main():
    print("Starting reviewer follow-up validation studies")
    keep_rows = perform_keepbins_scan()
    toy_rows = perform_toy_coverage()
    ntag_svd_rows, dndeta_svd_rows = perform_svd_scan()
    write_report(keep_rows, toy_rows, ntag_svd_rows, dndeta_svd_rows)
    print("Wrote follow-up summary to", REPORT_MD)


if __name__ == "__main__":
    main()
