import math

COUNTED_CHARGED_PDGS = {
    11,
    13,
    15,
    211,
    321,
    2212,
    3112,
    3222,
    3312,
    3334,
    411,
    431,
    521,
    541,
    24,
}

PDG_CHARGE_MAP = {
    11: -1.0,
    -11: +1.0,
    13: -1.0,
    -13: +1.0,
    15: -1.0,
    -15: +1.0,
    211: +1.0,
    -211: -1.0,
    321: +1.0,
    -321: -1.0,
    2212: +1.0,
    -2212: -1.0,
    3112: -1.0,
    -3112: +1.0,
    3222: +1.0,
    -3222: -1.0,
    3312: -1.0,
    -3312: +1.0,
    3334: -1.0,
    -3334: +1.0,
    411: +1.0,
    -411: -1.0,
    431: +1.0,
    -431: -1.0,
    521: +1.0,
    -521: -1.0,
    541: +1.0,
    -541: -1.0,
    24: +1.0,
    -24: -1.0,
}


def is_counted_charged_for_activity(pdg: int) -> bool:
    return abs(int(pdg)) in COUNTED_CHARGED_PDGS


def policy_charge_from_pdg(pdg: int) -> float:
    return PDG_CHARGE_MAP.get(int(pdg), 0.0)


def pass_pid_fiducial_from_mom(
    px: float,
    py: float,
    pz: float,
    use_pid_fiducial: bool = True,
    abs_cos_min: float = 0.15,
    abs_cos_max: float = 0.675,
) -> bool:
    if not use_pid_fiducial:
        return True
    p2 = px * px + py * py + pz * pz
    if p2 <= 0.0:
        return False
    abs_cos = abs(pz / math.sqrt(p2))
    return abs_cos > abs_cos_min and abs_cos < abs_cos_max


def pass_pt_window(px: float, py: float, pt_min: float = 0.4, pt_max: float = 5.0) -> bool:
    pt = math.hypot(px, py)
    return pt >= pt_min and pt < pt_max


def is_counted_pion_for_ratio(
    pdg: int,
    px: float,
    py: float,
    pz: float,
    use_pid_fiducial: bool = True,
    abs_cos_min: float = 0.15,
    abs_cos_max: float = 0.675,
    pt_min: float = 0.4,
    pt_max: float = 5.0,
) -> bool:
    return (
        abs(int(pdg)) == 211
        and pass_pt_window(px, py, pt_min, pt_max)
        and pass_pid_fiducial_from_mom(px, py, pz, use_pid_fiducial, abs_cos_min, abs_cos_max)
    )


def is_counted_kaon_for_ratio(
    pdg: int,
    px: float,
    py: float,
    pz: float,
    use_pid_fiducial: bool = True,
    abs_cos_min: float = 0.15,
    abs_cos_max: float = 0.675,
    pt_min: float = 0.4,
    pt_max: float = 5.0,
) -> bool:
    return (
        abs(int(pdg)) == 321
        and pass_pt_window(px, py, pt_min, pt_max)
        and pass_pid_fiducial_from_mom(px, py, pz, use_pid_fiducial, abs_cos_min, abs_cos_max)
    )


def is_counted_proton_for_ratio(
    pdg: int,
    px: float,
    py: float,
    pz: float,
    use_pid_fiducial: bool = True,
    abs_cos_min: float = 0.15,
    abs_cos_max: float = 0.675,
    pt_min: float = 0.4,
    pt_max: float = 5.0,
) -> bool:
    return (
        abs(int(pdg)) == 2212
        and pass_pt_window(px, py, pt_min, pt_max)
        and pass_pid_fiducial_from_mom(px, py, pz, use_pid_fiducial, abs_cos_min, abs_cos_max)
    )
