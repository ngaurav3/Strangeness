#!/usr/bin/env python3
import math
import os
import sys
from array import array

import ROOT

ROOT.gROOT.SetBatch(True)

from truth_counting_policy import (
    is_counted_charged_for_activity,
    is_counted_kaon_for_ratio,
    is_counted_pion_for_ratio,
    is_counted_proton_for_ratio,
    policy_charge_from_pdg,
)


def compute_eta(px: float, py: float, pz: float) -> float:
    p = math.sqrt(px * px + py * py + pz * pz)
    if p <= abs(pz):
        return 1e9 if pz >= 0 else -1e9
    return 0.5 * math.log((p + pz) / (p - pz))

def main() -> int:
    if len(sys.argv) != 3:
        print(
            "Usage: xscape_final_state_to_truth_root.py input_final_state_hadrons.dat output.root",
            file=sys.stderr,
        )
        return 1

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    db = ROOT.TDatabasePDG.Instance()

    fout = ROOT.TFile(output_path, "RECREATE")
    tree = ROOT.TTree("Events", "Converted X-SCAPE final-state hadrons")

    n_ch_eta05 = array("i", [0])
    n_ch_eta05_inclusive = array("i", [0])
    n_pi_pt0405 = array("i", [0])
    n_pi_pt0405_inclusive = array("i", [0])
    n_k_pt0405 = array("i", [0])
    n_k_pt0405_inclusive = array("i", [0])
    n_p_pt0405 = array("i", [0])
    k_pi_pt0405 = array("d", [0.0])
    k_pi_pt0405_inclusive = array("d", [0.0])

    tree.Branch("nChEta05", n_ch_eta05, "nChEta05/I")
    tree.Branch("nChEta05Inclusive", n_ch_eta05_inclusive, "nChEta05Inclusive/I")
    tree.Branch("nPiPt0405", n_pi_pt0405, "nPiPt0405/I")
    tree.Branch("nPiPt0405Inclusive", n_pi_pt0405_inclusive, "nPiPt0405Inclusive/I")
    tree.Branch("nKPt0405", n_k_pt0405, "nKPt0405/I")
    tree.Branch("nKPt0405Inclusive", n_k_pt0405_inclusive, "nKPt0405Inclusive/I")
    tree.Branch("nPPt0405", n_p_pt0405, "nPPt0405/I")
    tree.Branch("kPiPt0405", k_pi_pt0405, "kPiPt0405/D")
    tree.Branch("kPiPt0405Inclusive", k_pi_pt0405_inclusive, "kPiPt0405Inclusive/D")

    pdg = ROOT.std.vector("int")()
    is_weak_decay_daughter = ROOT.std.vector("int")()
    charge = ROOT.std.vector("float")()
    px = ROOT.std.vector("float")()
    py = ROOT.std.vector("float")()
    pz = ROOT.std.vector("float")()
    energy = ROOT.std.vector("float")()
    mass = ROOT.std.vector("float")()
    pt = ROOT.std.vector("float")()
    eta = ROOT.std.vector("float")()
    phi = ROOT.std.vector("float")()

    tree.Branch("pdg", pdg)
    tree.Branch("isWeakDecayDaughter", is_weak_decay_daughter)
    tree.Branch("charge", charge)
    tree.Branch("px", px)
    tree.Branch("py", py)
    tree.Branch("pz", pz)
    tree.Branch("e", energy)
    tree.Branch("m", mass)
    tree.Branch("pt", pt)
    tree.Branch("eta", eta)
    tree.Branch("phi", phi)

    h_nch = ROOT.TH1D("hNChEta05", "", 80, 0.0, 80.0)
    h_kpi = ROOT.TH1D("hKPiInclusive", "", 1, 0.0, 1.0)

    total_k = 0
    total_pi = 0
    event_count = 0

    def flush_event():
        nonlocal total_k, total_pi, event_count
        if event_count == 0 and n_ch_eta05[0] == 0 and n_pi_pt0405[0] == 0 and n_k_pt0405[0] == 0:
            return
        k_pi_pt0405[0] = (
            float(n_k_pt0405[0]) / float(n_pi_pt0405[0]) if n_pi_pt0405[0] > 0 else -1.0
        )
        k_pi_pt0405_inclusive[0] = (
            float(n_k_pt0405_inclusive[0]) / float(n_pi_pt0405_inclusive[0])
            if n_pi_pt0405_inclusive[0] > 0
            else -1.0
        )
        tree.Fill()
        h_nch.Fill(float(n_ch_eta05[0]))
        total_k += n_k_pt0405[0]
        total_pi += n_pi_pt0405[0]
        event_count += 1

    def reset_event():
        n_ch_eta05[0] = 0
        n_ch_eta05_inclusive[0] = 0
        n_pi_pt0405[0] = 0
        n_pi_pt0405_inclusive[0] = 0
        n_k_pt0405[0] = 0
        n_k_pt0405_inclusive[0] = 0
        n_p_pt0405[0] = 0
        k_pi_pt0405[0] = -1.0
        k_pi_pt0405_inclusive[0] = -1.0
        pdg.clear()
        is_weak_decay_daughter.clear()
        charge.clear()
        px.clear()
        py.clear()
        pz.clear()
        energy.clear()
        mass.clear()
        pt.clear()
        eta.clear()
        phi.clear()

    reset_event()
    seen_header = False

    with open(input_path, encoding="ascii") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#\tEvent"):
                if seen_header:
                    flush_event()
                    reset_event()
                seen_header = True
                continue
            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            pid = int(parts[1])
            ev = float(parts[3])
            pxv = float(parts[4])
            pyv = float(parts[5])
            pzv = float(parts[6])
            ptv = math.sqrt(pxv * pxv + pyv * pyv)
            etav = compute_eta(pxv, pyv, pzv)
            phiv = math.atan2(pyv, pxv)
            q = policy_charge_from_pdg(pid)
            particle = db.GetParticle(pid)
            mv = particle.Mass() if particle else max(ev * ev - pxv * pxv - pyv * pyv - pzv * pzv, 0.0) ** 0.5

            pdg.push_back(pid)
            is_weak_decay_daughter.push_back(0)
            charge.push_back(float(q))
            px.push_back(float(pxv))
            py.push_back(float(pyv))
            pz.push_back(float(pzv))
            energy.push_back(float(ev))
            mass.push_back(float(mv))
            pt.push_back(float(ptv))
            eta.push_back(float(etav))
            phi.push_back(float(phiv))

            if is_counted_charged_for_activity(pid) and abs(etav) < 0.5:
                n_ch_eta05[0] += 1
                n_ch_eta05_inclusive[0] += 1

            if is_counted_pion_for_ratio(pid, pxv, pyv, pzv):
                n_pi_pt0405[0] += 1
                n_pi_pt0405_inclusive[0] += 1
            elif is_counted_kaon_for_ratio(pid, pxv, pyv, pzv):
                n_k_pt0405[0] += 1
                n_k_pt0405_inclusive[0] += 1
            elif is_counted_proton_for_ratio(pid, pxv, pyv, pzv):
                n_p_pt0405[0] += 1

    if seen_header:
        flush_event()

    if total_pi > 0:
        h_kpi.SetBinContent(1, float(total_k) / float(total_pi))

    tree.Write()
    h_nch.Write()
    h_kpi.Write()
    fout.Close()

    summary_path = os.path.splitext(output_path)[0] + ".summary.txt"
    with open(summary_path, "w", encoding="ascii") as f:
        f.write(f"input = {input_path}\n")
        f.write(f"output = {output_path}\n")
        f.write(f"events = {event_count}\n")
        f.write(f"total_K_ch_pt0405 = {total_k}\n")
        f.write(f"total_pi_ch_pt0405 = {total_pi}\n")
        if total_pi > 0:
            f.write(f"inclusive_Kpi_pt0405 = {float(total_k) / float(total_pi):.8f}\n")

    print(f"Wrote ROOT: {output_path}")
    print(f"Wrote summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
