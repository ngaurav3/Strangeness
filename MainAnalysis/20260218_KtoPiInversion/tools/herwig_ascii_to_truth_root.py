#!/usr/bin/env python3
import math
import os
import sys
from array import array

import ROOT

ROOT.gROOT.SetBatch(True)

SPEED_OF_LIGHT_CM_PER_S = 2.99792458e10

from truth_counting_policy import (
    is_counted_charged_for_activity,
    is_counted_kaon_for_ratio,
    is_counted_pion_for_ratio,
    policy_charge_from_pdg,
)


def compute_eta(px, py, pz):
    pt = math.hypot(px, py)
    if pt <= 0.0:
        return 0.0
    return math.asinh(pz / pt)


def build_ctau_lookup():
    db = ROOT.TDatabasePDG.Instance()
    cache = {}

    def ctau_cm(pdg):
        apdg = abs(int(pdg))
        if apdg in cache:
            return cache[apdg]
        ctau = 0.0
        particle = db.GetParticle(apdg) if db else None
        if particle:
            lifetime = particle.Lifetime()
            if lifetime > 0.0 and math.isfinite(lifetime):
                ctau = lifetime * SPEED_OF_LIGHT_CM_PER_S
        cache[apdg] = ctau
        return ctau

    return ctau_cm


def parse_incoming_particles(token):
    token = token.strip()
    if not token.startswith("[") or not token.endswith("]"):
        return []
    body = token[1:-1].strip()
    if not body:
        return []
    return [int(x) for x in body.split(",") if x]


def finalize_event(event_number_value, particles, vertex_incoming, tree, buffers, h_nch, h_ratio, profile, ctau_cm):
    if event_number_value is None:
        return 0

    pdg = buffers["pdg"]
    status = buffers["status"]
    is_weak = buffers["isWeakDecayDaughter"]
    charge = buffers["charge"]
    px = buffers["px"]
    py = buffers["py"]
    pz = buffers["pz"]
    energy = buffers["e"]
    mass = buffers["m"]
    pt = buffers["pt"]
    eta = buffers["eta"]
    phi = buffers["phi"]

    pdg.clear()
    status.clear()
    is_weak.clear()
    charge.clear()
    px.clear()
    py.clear()
    pz.clear()
    energy.clear()
    mass.clear()
    pt.clear()
    eta.clear()
    phi.clear()

    memo = {}
    visiting = set()

    def has_long_lived_ancestor(barcode):
        if barcode in memo:
            return memo[barcode]
        if barcode in visiting:
            return False
        visiting.add(barcode)
        particle = particles.get(barcode)
        result = False
        if particle:
            prod_vtx = particle["prod_vtx"]
            if prod_vtx in vertex_incoming:
                for mother_barcode in vertex_incoming[prod_vtx]:
                    mother = particles.get(mother_barcode)
                    if mother is None:
                        continue
                    if ctau_cm(mother["pid"]) > 1.0 or has_long_lived_ancestor(mother_barcode):
                        result = True
                        break
        visiting.discard(barcode)
        memo[barcode] = result
        return result

    buffers["eventNumber"][0] = int(event_number_value)
    buffers["ecm"][0] = 91.2
    buffers["nFinal"][0] = 0
    buffers["nCharged"][0] = 0
    buffers["nChargedInclusive"][0] = 0
    buffers["nChEta05"][0] = 0
    buffers["nChEta05Inclusive"][0] = 0
    buffers["nPiPt0405"][0] = 0
    buffers["nPiPt0405Inclusive"][0] = 0
    buffers["nKPt0405"][0] = 0
    buffers["nKPt0405Inclusive"][0] = 0
    buffers["kPiPt0405"][0] = -1.0
    buffers["kPiPt0405Inclusive"][0] = -1.0

    for barcode in sorted(particles):
        particle = particles[barcode]
        if particle["status"] != 1:
            continue

        weak_daughter = has_long_lived_ancestor(barcode)
        q = policy_charge_from_pdg(particle["pid"])
        pxv = particle["px"]
        pyv = particle["py"]
        pzv = particle["pz"]
        ev = particle["e"]
        mv = particle["m"]
        ptv = math.hypot(pxv, pyv)
        etav = compute_eta(pxv, pyv, pzv)
        phiv = math.atan2(pyv, pxv)

        buffers["nFinal"][0] += 1
        pdg.push_back(int(particle["pid"]))
        status.push_back(int(particle["status"]))
        is_weak.push_back(1 if weak_daughter else 0)
        charge.push_back(float(q))
        px.push_back(float(pxv))
        py.push_back(float(pyv))
        pz.push_back(float(pzv))
        energy.push_back(float(ev))
        mass.push_back(float(mv))
        pt.push_back(float(ptv))
        eta.push_back(float(etav))
        phi.push_back(float(phiv))

        if q == 0.0:
            continue
        buffers["nChargedInclusive"][0] += 1
        if not weak_daughter:
            buffers["nCharged"][0] += 1
        if abs(etav) < 0.5:
            buffers["nChEta05Inclusive"][0] += 1
            if not weak_daughter:
                buffers["nChEta05"][0] += 1

        if is_counted_pion_for_ratio(particle["pid"], pxv, pyv, pzv):
            buffers["nPiPt0405Inclusive"][0] += 1
            if not weak_daughter:
                buffers["nPiPt0405"][0] += 1
        if is_counted_kaon_for_ratio(particle["pid"], pxv, pyv, pzv):
            buffers["nKPt0405Inclusive"][0] += 1
            if not weak_daughter:
                buffers["nKPt0405"][0] += 1

    if buffers["nPiPt0405"][0] > 0:
        buffers["kPiPt0405"][0] = float(buffers["nKPt0405"][0]) / float(buffers["nPiPt0405"][0])
    if buffers["nPiPt0405Inclusive"][0] > 0:
        buffers["kPiPt0405Inclusive"][0] = float(buffers["nKPt0405Inclusive"][0]) / float(buffers["nPiPt0405Inclusive"][0])

    tree.Fill()
    h_nch.Fill(buffers["nChEta05"][0])
    if buffers["kPiPt0405"][0] >= 0.0:
        h_ratio.Fill(buffers["kPiPt0405"][0])
        profile.Fill(buffers["nChEta05"][0], buffers["kPiPt0405"][0])
    return 1


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} input.hepmc output.root", file=sys.stderr)
        return 1

    input_name = sys.argv[1]
    output_name = sys.argv[2]
    if not os.path.exists(input_name):
        print(f"Missing input: {input_name}", file=sys.stderr)
        return 1

    out_file = ROOT.TFile(output_name, "RECREATE")
    tree = ROOT.TTree("Events", "Herwig truth events for K/pi vs dNch/deta")

    buffers = {
        "eventNumber": array("i", [0]),
        "ecm": array("d", [91.2]),
        "nFinal": array("i", [0]),
        "nCharged": array("i", [0]),
        "nChargedInclusive": array("i", [0]),
        "nChEta05": array("i", [0]),
        "nChEta05Inclusive": array("i", [0]),
        "nPiPt0405": array("i", [0]),
        "nPiPt0405Inclusive": array("i", [0]),
        "nKPt0405": array("i", [0]),
        "nKPt0405Inclusive": array("i", [0]),
        "kPiPt0405": array("d", [-1.0]),
        "kPiPt0405Inclusive": array("d", [-1.0]),
        "pdg": ROOT.std.vector("int")(),
        "status": ROOT.std.vector("int")(),
        "isWeakDecayDaughter": ROOT.std.vector("int")(),
        "charge": ROOT.std.vector("float")(),
        "px": ROOT.std.vector("float")(),
        "py": ROOT.std.vector("float")(),
        "pz": ROOT.std.vector("float")(),
        "e": ROOT.std.vector("float")(),
        "m": ROOT.std.vector("float")(),
        "pt": ROOT.std.vector("float")(),
        "eta": ROOT.std.vector("float")(),
        "phi": ROOT.std.vector("float")(),
    }

    tree.Branch("eventNumber", buffers["eventNumber"], "eventNumber/I")
    tree.Branch("ecm", buffers["ecm"], "ecm/D")
    tree.Branch("nFinal", buffers["nFinal"], "nFinal/I")
    tree.Branch("nCharged", buffers["nCharged"], "nCharged/I")
    tree.Branch("nChargedInclusive", buffers["nChargedInclusive"], "nChargedInclusive/I")
    tree.Branch("nChEta05", buffers["nChEta05"], "nChEta05/I")
    tree.Branch("nChEta05Inclusive", buffers["nChEta05Inclusive"], "nChEta05Inclusive/I")
    tree.Branch("nPiPt0405", buffers["nPiPt0405"], "nPiPt0405/I")
    tree.Branch("nPiPt0405Inclusive", buffers["nPiPt0405Inclusive"], "nPiPt0405Inclusive/I")
    tree.Branch("nKPt0405", buffers["nKPt0405"], "nKPt0405/I")
    tree.Branch("nKPt0405Inclusive", buffers["nKPt0405Inclusive"], "nKPt0405Inclusive/I")
    tree.Branch("kPiPt0405", buffers["kPiPt0405"], "kPiPt0405/D")
    tree.Branch("kPiPt0405Inclusive", buffers["kPiPt0405Inclusive"], "kPiPt0405Inclusive/D")
    tree.Branch("pdg", buffers["pdg"])
    tree.Branch("status", buffers["status"])
    tree.Branch("isWeakDecayDaughter", buffers["isWeakDecayDaughter"])
    tree.Branch("charge", buffers["charge"])
    tree.Branch("px", buffers["px"])
    tree.Branch("py", buffers["py"])
    tree.Branch("pz", buffers["pz"])
    tree.Branch("e", buffers["e"])
    tree.Branch("m", buffers["m"])
    tree.Branch("pt", buffers["pt"])
    tree.Branch("eta", buffers["eta"])
    tree.Branch("phi", buffers["phi"])

    h_nch = ROOT.TH1D(
        "hNChEta05",
        "Charged multiplicity in |#eta|<0.5;dN_{ch}/d#eta (|#eta|<0.5);Events",
        60,
        -0.5,
        59.5,
    )
    h_ratio = ROOT.TH1D(
        "hKPiInclusive",
        "Inclusive K/#pi (0.4 < p_{T} < 5.0 GeV/c);K/#pi;Events",
        100,
        0.0,
        1.0,
    )
    profile = ROOT.TProfile(
        "pKPiVsDNdEta",
        "K/#pi vs dN_{ch}/d#eta;dN_{ch}/d#eta (|#eta|<0.5);K/#pi",
        20,
        -0.5,
        39.5,
    )

    ctau_cm = build_ctau_lookup()

    accepted_events = 0
    current_event_number = None
    current_particles = {}
    current_vertex_incoming = {}

    with open(input_name, encoding="ascii", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("HepMC::"):
                continue
            prefix = line[0]
            if prefix == "E":
                accepted_events += finalize_event(
                    current_event_number,
                    current_particles,
                    current_vertex_incoming,
                    tree,
                    buffers,
                    h_nch,
                    h_ratio,
                    profile,
                    ctau_cm,
                )
                parts = line.split()
                current_event_number = int(parts[1])
                current_particles = {}
                current_vertex_incoming = {}
            elif prefix == "P":
                parts = line.split()
                if len(parts) < 10:
                    continue
                barcode = int(parts[1])
                current_particles[barcode] = {
                    "prod_vtx": int(parts[2]),
                    "pid": int(parts[3]),
                    "px": float(parts[4]),
                    "py": float(parts[5]),
                    "pz": float(parts[6]),
                    "e": float(parts[7]),
                    "m": float(parts[8]),
                    "status": int(parts[9]),
                }
            elif prefix == "V":
                parts = line.split()
                if len(parts) < 4:
                    continue
                vertex_id = int(parts[1])
                current_vertex_incoming[vertex_id] = parse_incoming_particles(parts[3])

    accepted_events += finalize_event(
        current_event_number,
        current_particles,
        current_vertex_incoming,
        tree,
        buffers,
        h_nch,
        h_ratio,
        profile,
        ctau_cm,
    )

    out_file.cd()
    tree.Write()
    h_nch.Write()
    h_ratio.Write()
    profile.Write()
    ROOT.TNamed(
        "weakDecayDaughterPolicy",
        "Default nChEta05/nPiPt0405/nKPt0405 branches veto daughters of ancestors with ctau>1 cm; inclusive bookkeeping retained in *Inclusive branches",
    ).Write()
    out_file.Close()

    print(f"Wrote Herwig truth sample to {os.path.abspath(output_name)}")
    print(f"Accepted events: {accepted_events}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
