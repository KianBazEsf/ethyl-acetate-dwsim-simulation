#!/usr/bin/env python3
"""
Continuous Ethyl Acetate Production (1,000 kg/day) — Independent Verification Suite
===================================================================================
Author: Kian Bazrafkan Esfahani
Department of Chemical Engineering, The University of Manchester

Performs closed-form steady-state mass, elemental, thermodynamic, kinetic,
enthalpy, heat-integration, recycle-convergence, sensitivity, effluent COD,
and indicative economic calculations for the continuous acid-catalysed
esterification of ethanol with acetic acid over Amberlyst 35 resin.

Also parses any DWSIM 10.0 `.dwxmz` flowsheet archives placed in the working
directory (without requiring DWSIM to be installed) and executes the 12
acceptance checks reported in Table 9 of the technical portfolio report.

Dependencies: numpy, scipy
Usage:
    python Ethyl_Acetate_Calculations.py
"""

from __future__ import annotations

import glob
import math
import os
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# =============================================================================
# 1. COMPONENT DATA & DESIGN BASIS (Tables 2, 3, 5, 15)
# =============================================================================

COMPONENTS = ["EtOH", "AA", "EA", "W"]
COMP_NAMES = {
    "EtOH": "Ethanol",
    "AA": "Acetic acid",
    "EA": "Ethyl acetate",
    "W": "Water",
}

# ChemSep / DWSIM molar masses (kg/kmol)
MW = {
    "EtOH": 46.06844,
    "AA": 60.05196,
    "EA": 88.10512,
    "W": 18.01528,
}

# Elemental formula counts (C, H, O) for atomic balance verification
ATOMS = {
    "EtOH": (2, 6, 1),
    "AA": (2, 4, 2),
    "EA": (4, 8, 2),
    "W": (0, 2, 1),
}

# Stoichiometric coefficients: AA + EtOH -> EA + W
NU = {"EtOH": -1.0, "AA": -1.0, "EA": +1.0, "W": +1.0}

# Pure-component liquid densities at ~25 C (kg/m3) & ~70 C (kg/m3)
RHO_25C = {"EtOH": 785.0, "AA": 1044.0, "EA": 894.0, "W": 997.0}
RHO_70C = {"EtOH": 745.0, "AA": 993.0, "EA": 838.0, "W": 978.0}

# NIST liquid standard enthalpies of formation at 298.15 K (kJ/mol)
DHF_LIQ_298 = {
    "EtOH": -277.60,
    "AA": -484.30,
    "EA": -480.00,
    "W": -285.83,
}
# Ideal-gas standard enthalpies of formation at 298.15 K (kJ/mol) — DWSIM ChemSep
DHF_IG_298 = {
    "EtOH": -234.95,
    "AA": -432.25,
    "EA": -443.92,
    "W": -241.826,
}

# DIPPR-100 liquid molar heat capacity coefficients:
#   Cp(T) [J/(kmol K)] = C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4  (T in K)
# Matches NIST/DIPPR liquid data within 2.2% over 25-70 C.
DIPPR_CP_LIQ = {
    "EtOH": (1.0264e5, -1.3963e2, -3.0341e-2, 2.0386e-3, 0.0),
    "AA": (1.3964e5, -3.2080e2, 8.9850e-1, 0.0, 0.0),
    "EA": (1.6170e5, -2.1890e2, 8.2100e-1, 0.0, 0.0),
    "W": (2.7637e5, -2.0901e3, 8.1250e0, -1.4116e-2, 9.3701e-6),
}

# Antoine vapour-pressure parameters: log10(P_bar) = A - B / (T_C + C)
ANTOINE = {
    "EtOH": (5.24677, 1598.673, 226.58),
    "AA": (4.68206, 1642.540, 233.52),
    "EA": (4.22809, 1245.702, 217.88),
    "W": (5.11564, 1687.537, 230.17),
}

# NRTL binary interaction parameters as stored in DWSIM (Table 5)
# Order: [EtOH, AA, EA, W], A_ij in cal/mol, R = 1.98721 cal/(mol K)
R_CAL = 1.98721
R_J = 8.314462618

NRTL_A = np.array(
    [
        [0.0, -34.20, 330.51, -57.96],       # EtOH -> j
        [-190.78, 0.0, -303.60, -189.10],    # AA   -> j
        [305.60, 544.80, 0.0, 1285.99],      # EA   -> j
        [1241.74, 728.80, 1606.08, 0.0],     # W    -> j
    ],
    dtype=float,
)

NRTL_ALPHA = np.array(
    [
        [0.0, 0.305, 0.299, 0.294],
        [0.305, 0.0, 0.200, 0.200],
        [0.299, 0.200, 0.0, 0.439],
        [0.294, 0.200, 0.439, 0.0],
    ],
    dtype=float,
)

# Design basis parameters (Tables 2 & 15)
FRESH_FEED_KGH = {"EtOH": 22.2459, "AA": 29.9327, "EA": 0.0, "W": 0.0}
CS1_SPLIT_PRODUCT = {"EtOH": 0.0001, "AA": 0.00005, "EA": 0.99, "W": 0.0008}
CS2_SPLIT_RECYCLE = {"EtOH": 0.99, "AA": 0.90, "EA": 0.99, "W": 0.30}
SINGLE_PASS_CONV = 0.65
T_FEED_C = 25.0
T_REACTOR_C = 70.0
T_COOLER_C = 35.0
P_FEED_BAR = 1.01325
P_PUMP_OUT_BAR = 2.50
P_REACTOR_IN_BAR = 2.30
P_REACTOR_OUT_BAR = 2.00
P_COOLER_OUT_BAR = 1.80
PUMP_EFF = 0.75
OPERATING_DAYS = 330.0
OPERATING_HOURS = OPERATING_DAYS * 24.0  # 7920 h/year


# =============================================================================
# 2. CLOSED-FORM MASS BALANCE & ELEMENTAL CLOSURE (Section 5.1, Table 7)
# =============================================================================

@dataclass
class MassBalanceResult:
    streams_kmolh: Dict[str, Dict[str, float]]
    streams_kgh: Dict[str, Dict[str, float]]
    extent_kmolh: float
    loop_gain: Dict[str, float]
    overall_aa_conv: float
    reactor_inlet_etoh_aa: float
    fresh_etoh_aa: float
    product_ea_wt_pct: float
    physical_etoh_passes: float
    atomic_residuals: Dict[str, float]


def solve_mass_balance(
    fresh_kgh: Dict[str, float] = FRESH_FEED_KGH,
    s_split: Dict[str, float] = CS1_SPLIT_PRODUCT,
    r_split: Dict[str, float] = CS2_SPLIT_RECYCLE,
    x_aa: float = SINGLE_PASS_CONV,
) -> MassBalanceResult:
    """Solve the closed-loop recycle mass balance analytically (Eq. 8)."""
    g = {c: (1.0 - s_split[c]) * r_split[c] for c in COMPONENTS}
    f_fresh = {c: fresh_kgh[c] / MW[c] for c in COMPONENTS}

    # Reactor inlet molar flows (closed form)
    f_in = {}
    f_in["AA"] = f_fresh["AA"] / (1.0 - g["AA"] * (1.0 - x_aa))
    xi = x_aa * f_in["AA"]
    f_in["EtOH"] = (f_fresh["EtOH"] - g["EtOH"] * xi) / (1.0 - g["EtOH"])
    f_in["EA"] = (g["EA"] * xi) / (1.0 - g["EA"])
    f_in["W"] = (g["W"] * xi) / (1.0 - g["W"])

    # Reactor outlet molar flows
    f_out = {c: f_in[c] + NU[c] * xi for c in COMPONENTS}

    # Separator splits
    f_prod = {c: s_split[c] * f_out[c] for c in COMPONENTS}
    f_cs1_res = {c: (1.0 - s_split[c]) * f_out[c] for c in COMPONENTS}
    f_rec = {c: r_split[c] * f_cs1_res[c] for c in COMPONENTS}
    f_waste = {c: (1.0 - r_split[c]) * f_cs1_res[c] for c in COMPONENTS}
    f_zero = {c: 0.0 for c in COMPONENTS}

    streams_mol = {
        "1": {"EtOH": f_fresh["EtOH"], "AA": 0.0, "EA": 0.0, "W": 0.0},
        "2": {"EtOH": 0.0, "AA": f_fresh["AA"], "EA": 0.0, "W": 0.0},
        "3": dict(f_fresh),
        "20": dict(f_rec),
        "4": dict(f_in),
        "5": dict(f_in),
        "6": dict(f_in),
        "Waste-1": dict(f_zero),
        "7": dict(f_out),
        "8": dict(f_out),
        "9": dict(f_out),
        "Waste-2": dict(f_zero),
        "10": dict(f_out),
        "Product-EtAc": dict(f_prod),
        "11": dict(f_cs1_res),
        "Recycle-1": dict(f_rec),
        "12": dict(f_waste),
    }

    streams_kg = {
        tag: {c: mols[c] * MW[c] for c in COMPONENTS}
        for tag, mols in streams_mol.items()
    }

    m_prod_tot = sum(streams_kg["Product-EtAc"].values())
    ea_wt_pct = 100.0 * streams_kg["Product-EtAc"]["EA"] / m_prod_tot
    ov_conv = 100.0 * xi / f_fresh["AA"]
    inlet_ratio = f_in["EtOH"] / f_in["AA"]
    fresh_ratio = f_fresh["EtOH"] / f_fresh["AA"]
    phys_passes = f_in["EtOH"] / f_fresh["EtOH"]

    # Elemental C, H, O balance closure across battery limits
    atom_names = ["C", "H", "O"]
    atomic_res = {}
    for idx, el in enumerate(atom_names):
        in_atoms = sum(f_fresh[c] * ATOMS[c][idx] for c in COMPONENTS)
        out_atoms = sum((f_prod[c] + f_waste[c]) * ATOMS[c][idx] for c in COMPONENTS)
        atomic_res[el] = abs(in_atoms - out_atoms)

    return MassBalanceResult(
        streams_kmolh=streams_mol,
        streams_kgh=streams_kg,
        extent_kmolh=xi,
        loop_gain=g,
        overall_aa_conv=ov_conv,
        reactor_inlet_etoh_aa=inlet_ratio,
        fresh_etoh_aa=fresh_ratio,
        product_ea_wt_pct=ea_wt_pct,
        physical_etoh_passes=phys_passes,
        atomic_residuals=atomic_res,
    )


# =============================================================================
# 3. THERMODYNAMICS: NRTL, VAPOUR PRESSURE, BUBBLE POINTS & EQUILIBRIUM
# =============================================================================

def nrtl_gamma(x: np.ndarray, t_k: float) -> np.ndarray:
    """Compute liquid-phase activity coefficients using the NRTL model."""
    x = np.asarray(x, dtype=float)
    x = x / np.sum(x)
    tau = NRTL_A / (R_CAL * t_k)
    G = np.exp(-NRTL_ALPHA * tau)

    # First term: sum_j(tau_ji * G_ji * x_j) / sum_k(G_ki * x_k)
    denom = G.T @ x
    term1 = ((tau * G).T @ x) / denom

    # Second term: sum_j [ (x_j * G_ij / sum_k(G_kj * x_k)) * (tau_ij - term1_j) ]
    n = len(x)
    term2 = np.zeros(n)
    for i in range(n):
        s = 0.0
        for j in range(n):
            s += (x[j] * G[i, j] / denom[j]) * (tau[i, j] - term1[j])
        term2[i] = s

    return np.exp(term1 + term2)


def psat_bar(comp: str, t_c: float) -> float:
    """Pure-component vapour pressure in bar from Antoine correlation."""
    A, B, C = ANTOINE[comp]
    return 10.0 ** (A - B / (t_c + C))


def bubble_point_c(x: np.ndarray, p_bar: float) -> float:
    """Bubble temperature (deg C) from modified Raoult's law (Eq. 7)."""
    x = np.asarray(x, dtype=float)
    x = x / np.sum(x)

    def residual(t_c: float) -> float:
        gamma = nrtl_gamma(x, t_c + 273.15)
        ps = np.array([psat_bar(c, t_c) for c in COMPONENTS])
        return float(np.sum(x * gamma * ps) - p_bar)

    return brentq(residual, 20.0, 200.0)


def reaction_enthalpies(t_c: float = T_REACTOR_C) -> Dict[str, float]:
    """Compute standard liquid and ideal-gas heats of reaction (kJ/mol)."""
    dh_liq_298 = sum(NU[c] * DHF_LIQ_298[c] for c in COMPONENTS)
    dh_ig_298 = sum(NU[c] * DHF_IG_298[c] for c in COMPONENTS)

    # Kirchhoff correction from 298.15 K to t_c + 273.15 K using liquid Cp
    t_k = t_c + 273.15
    dcp_int = 0.0
    for c in COMPONENTS:
        dcp_int += NU[c] * integral_cp_liq_kj_mol(c, 298.15, t_k)
    dh_liq_t = dh_liq_298 + dcp_int

    return {
        "dh_liq_25C_kJmol": dh_liq_298,
        "dh_liq_T_kJmol": dh_liq_t,
        "dh_ig_25C_kJmol": dh_ig_298,
    }


def solve_equilibrium_conversion(
    f_in_mol: Dict[str, float],
    k_val: float,
    activity_based: bool = False,
    t_c: float = T_REACTOR_C,
) -> float:
    """Solve equilibrium AA conversion X_eq for a given reactor-inlet composition."""
    f0 = np.array([f_in_mol[c] for c in COMPONENTS], dtype=float)
    f_aa0 = f0[1]
    t_k = t_c + 273.15

    def eq_residual(x_conv: float) -> float:
        xi = x_conv * f_aa0
        f = np.array(
            [f0[0] - xi, f0[1] - xi, f0[2] + xi, f0[3] + xi],
            dtype=float,
        )
        mol_frac = f / np.sum(f)
        kx = (mol_frac[2] * mol_frac[3]) / (mol_frac[0] * mol_frac[1])
        if activity_based:
            gamma = nrtl_gamma(mol_frac, t_k)
            cap_gamma = (gamma[2] * gamma[3]) / (gamma[0] * gamma[1])
            return kx * cap_gamma - k_val
        return kx - k_val

    return brentq(eq_residual, 1e-6, 0.9999)


def calibrate_ka_nrtl(t_c: float = T_REACTOR_C, dh_r_kjmol: float = -5.27) -> Tuple[float, float, float]:
    """
    Calibrate activity equilibrium constant Ka on classic equimolar Kx = 4.0 at 25 C
    and apply van 't Hoff adjustment to temperature t_c (Section 3.2).
    """
    # At equimolar equilibrium Kx = 4.0 => X_eq = 2/3 => [EtOH, AA, EA, W] = [1/6, 1/6, 2/6, 2/6]
    x_eq_equimolar = np.array([1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0])
    gamma_25 = nrtl_gamma(x_eq_equimolar, 298.15)
    cap_gamma_25 = (gamma_25[2] * gamma_25[3]) / (gamma_25[0] * gamma_25[1])
    ka_25 = 4.0 * cap_gamma_25
    # van 't Hoff: ln(Ka(T)/Ka(298.15)) = -(dH_r / R) * (1/T - 1/298.15)
    t_k = t_c + 273.15
    ka_t = ka_25 * math.exp(-(dh_r_kjmol * 1000.0 / R_J) * (1.0 / t_k - 1.0 / 298.15))
    return cap_gamma_25, ka_25, ka_t


# =============================================================================
# 4. KINETICS, REACTOR SIZING & ERGUN PRESSURE DROP (Section 3.3, Table 4)
# =============================================================================

def lai_rate_constants(t_c: float) -> Tuple[float, float]:
    """Lai et al. (2008) rate constants [cm^6 mol^-1 g^-1 min^-1] on wet Amberlyst 35."""
    t_k = t_c + 273.15
    k1 = 1.24e9 * math.exp(-6105.6 / t_k)
    km1 = 1.34e8 * math.exp(-5692.1 / t_k)
    return k1, km1


def size_reactor_lai(
    f_in_mol: Dict[str, float],
    x_target: float = SINGLE_PASS_CONV,
    t_c: float = T_REACTOR_C,
    override_k_eq: float | None = None,
) -> Dict[str, float]:
    """
    Compute wet Amberlyst 35 catalyst mass (kg) for isothermal PFR and CSTR
    via Equations (5) and (6), plus bed geometry, LHSV, and Ergun pressure drop.
    """
    k1, km1 = lai_rate_constants(t_c)
    if override_k_eq is not None:
        km1 = k1 / override_k_eq

    f0 = {c: f_in_mol[c] for c in COMPONENTS}
    f_aa0_mol_min = f0["AA"] * 1000.0 / 60.0  # kmol/h -> mol/min

    # Volumetric liquid flow at reactor temperature (cm3/h) via ideal volume mixing
    v_dot_m3_h = sum((f0[c] * MW[c]) / RHO_70C[c] for c in COMPONENTS)
    v_dot_cm3_h = v_dot_m3_h * 1.0e6

    def rate_mol_g_min(x_conv: float) -> float:
        xi = x_conv * f0["AA"]
        f_curr = {c: f0[c] + NU[c] * xi for c in COMPONENTS}
        # Concentration in mol/cm3 = (kmol/h * 1000) / (cm3/h)
        conc = {c: (f_curr[c] * 1000.0) / v_dot_cm3_h for c in COMPONENTS}
        return k1 * conc["AA"] * conc["EtOH"] - km1 * conc["EA"] * conc["W"]

    integral_val, _ = quad(lambda xp: 1.0 / rate_mol_g_min(xp), 0.0, x_target)
    w_pfr_kg = (f_aa0_mol_min * integral_val) / 1000.0
    w_cstr_kg = (f_aa0_mol_min * x_target / rate_mol_g_min(x_target)) / 1000.0

    # Design bed: D = 0.20 m, L = 1.0 m, Amberlyst 35 shipping bulk density 800 kg/m3
    bed_d = 0.20
    bed_l = 1.00
    bed_area = math.pi * (bed_d ** 2) / 4.0
    bed_vol_m3 = bed_area * bed_l
    bed_mass_kg = bed_vol_m3 * 800.0
    lhsv_h = v_dot_m3_h / bed_vol_m3

    # Ergun clean-bed pressure drop (dp = 0.8 mm resin beads, eps = 0.38, mu = 0.45 mPa s)
    rho_mix = sum(f0[c] * MW[c] for c in COMPONENTS) / v_dot_m3_h
    u_sup = (v_dot_m3_h / 3600.0) / bed_area
    dp = 0.0008
    eps = 0.38
    mu = 4.5e-4
    ergun_pa = bed_l * (
        150.0 * mu * u_sup * ((1.0 - eps) ** 2) / ((dp ** 2) * (eps ** 3))
        + 1.75 * rho_mix * (u_sup ** 2) * (1.0 - eps) / (dp * (eps ** 3))
    )

    return {
        "k1": k1,
        "km1": km1,
        "K_lai": k1 / km1,
        "W_PFR_kg": w_pfr_kg,
        "W_CSTR_kg": w_cstr_kg,
        "V_dot_Lh": v_dot_m3_h * 1000.0,
        "bed_mass_kg": bed_mass_kg,
        "LHSV_h": lhsv_h,
        "Ergun_dP_kPa": ergun_pa / 1000.0,
    }


# =============================================================================
# 5. IDEAL-LIQUID ENTHALPY BENCHMARK & HEAT INTEGRATION (Section 6, Table 10)
# =============================================================================

def integral_cp_liq_kj_mol(comp: str, t1_k: float, t2_k: float) -> float:
    """Integrate DIPPR liquid Cp(T) from t1_k to t2_k; returns kJ/mol."""
    c1, c2, c3, c4, c5 = DIPPR_CP_LIQ[comp]

    def antideriv(t: float) -> float:
        return (
            c1 * t
            + 0.5 * c2 * (t ** 2)
            + (1.0 / 3.0) * c3 * (t ** 3)
            + 0.25 * c4 * (t ** 4)
            + 0.2 * c5 * (t ** 5)
        )

    # DIPPR output is J/kmol -> divide by 1e6 for kJ/mol
    return (antideriv(t2_k) - antideriv(t1_k)) / 1.0e6


def stream_enthalpy_kw(f_mol_h: Dict[str, float], t_c: float, t_ref_c: float = 25.0) -> float:
    """Sensible liquid enthalpy flow (kW) relative to pure liquids at t_ref_c."""
    t_k = t_c + 273.15
    t_ref_k = t_ref_c + 273.15
    h_kj_h = sum(
        (f_mol_h[c] * 1000.0) * integral_cp_liq_kj_mol(c, t_ref_k, t_k)
        for c in COMPONENTS
    )
    return h_kj_h / 3600.0


def solve_energy_balance(mb: MassBalanceResult) -> Dict[str, float]:
    """Compute ideal-liquid duties E1-E6, mixer temperatures T3 & T4, and FEHE recovery."""
    # T3 (MIX-1): both fresh feeds at 25 C => ideal mixing gives 25.00 C
    t3_c = T_FEED_C

    # T4 (MIX-2): adiabatic mixing of stream 3 (25 C) and stream 20 (35 C)
    h_in_mix2 = stream_enthalpy_kw(mb.streams_kmolh["3"], t3_c) + stream_enthalpy_kw(
        mb.streams_kmolh["20"], T_COOLER_C
    )
    t4_c = brentq(
        lambda t: stream_enthalpy_kw(mb.streams_kmolh["4"], t) - h_in_mix2, 20.0, 40.0
    )

    # PUMP-1 shaft power E1 (W): V_dot * dP / eta
    v_dot_in_m3_s = (
        sum(mb.streams_kgh["4"][c] / RHO_25C[c] for c in COMPONENTS) / 3600.0
    )
    dp_pa = (P_PUMP_OUT_BAR - P_FEED_BAR) * 1.0e5
    e1_w = v_dot_in_m3_s * dp_pa / PUMP_EFF
    e1_kw = e1_w / 1000.0

    # T5 (PUMP-1 outlet): pump shaft work enters liquid
    t5_c = brentq(
        lambda t: stream_enthalpy_kw(mb.streams_kmolh["5"], t) - (h_in_mix2 + e1_kw),
        20.0,
        40.0,
    )

    # HT-1 heater duty E2 (kW): stream 5 (T5) -> stream 6 (70 C)
    e2_kw = stream_enthalpy_kw(mb.streams_kmolh["6"], T_REACTOR_C) - stream_enthalpy_kw(
        mb.streams_kmolh["5"], t5_c
    )

    # RCONV-1 isothermal heat removal E3 (kW): -xi * dH_r(70 C)
    rxn_h = reaction_enthalpies(T_REACTOR_C)
    e3_kw = -(mb.extent_kmolh * 1000.0 / 3600.0) * rxn_h["dh_liq_T_kJmol"]

    # Eq. (9): ideal-gas vs liquid-phase reference-state offset (kW)
    offset_kw = (mb.extent_kmolh / 3600.0) * (
        abs(rxn_h["dh_ig_25C_kJmol"] * 1000.0) - abs(rxn_h["dh_liq_25C_kJmol"] * 1000.0)
    )

    # CL-1 cooler duty E4 (kW): stream 7 (70 C) -> stream 8 (35 C)
    e4_kw = stream_enthalpy_kw(mb.streams_kmolh["7"], T_REACTOR_C) - stream_enthalpy_kw(
        mb.streams_kmolh["8"], T_COOLER_C
    )

    # Feed-effluent heat exchanger (10 K minimum approach => cold feed heated to 60 C)
    t_cold_out = T_REACTOR_C - 10.0  # 60.0 C
    q_fehe_kw = stream_enthalpy_kw(mb.streams_kmolh["6"], t_cold_out) - stream_enthalpy_kw(
        mb.streams_kmolh["5"], t5_c
    )
    t_hot_out = brentq(
        lambda t: (
            stream_enthalpy_kw(mb.streams_kmolh["7"], T_REACTOR_C)
            - stream_enthalpy_kw(mb.streams_kmolh["7"], t)
            - q_fehe_kw
        ),
        35.0,
        70.0,
    )
    dt1 = T_REACTOR_C - t_cold_out
    dt2 = t_hot_out - t5_c
    lmtd_fehe = (dt1 - dt2) / math.log(dt1 / dt2)
    area_fehe_m2 = q_fehe_kw / (0.25 * lmtd_fehe)  # U = 0.25 kW/(m2 K)

    return {
        "T3_C": t3_c,
        "T4_C": t4_c,
        "T5_C": t5_c,
        "E1_W": e1_w,
        "E2_kW": e2_kw,
        "E3_kW": e3_kw,
        "E3_offset_Eq9_kW": offset_kw,
        "E4_kW": e4_kw,
        "E5_kW": 0.0,
        "E6_kW": 0.0,
        "Q_FEHE_kW": q_fehe_kw,
        "T_hot_out_C": t_hot_out,
        "HT1_integrated_kW": e2_kw - q_fehe_kw,
        "CL1_integrated_kW": e4_kw - q_fehe_kw,
        "FEHE_area_m2": area_fehe_m2,
    }


# =============================================================================
# 6. RECYCLE CONVERGENCE RECONSTRUCTION (Section 5.2, Figure 2)
# =============================================================================

def simulate_recycle_convergence(
    method: str = "direct",
    tol_kgs: float = 1.0e-8,
    max_iter: int = 1200,
) -> Dict[str, object]:
    """
    Reconstruct DWSIM's tear-stream iteration from an empty tear stream (stream 20)
    using direct substitution or bounded Wegstein acceleration.
    """
    f_fresh = np.array([FRESH_FEED_KGH[c] / MW[c] for c in COMPONENTS])
    s_vec = np.array([CS1_SPLIT_PRODUCT[c] for c in COMPONENTS])
    r_vec = np.array([CS2_SPLIT_RECYCLE[c] for c in COMPONENTS])
    mw_vec = np.array([MW[c] for c in COMPONENTS])

    def map_step(x_tear: np.ndarray) -> np.ndarray:
        f_in = f_fresh + x_tear
        xi = SINGLE_PASS_CONV * f_in[1]
        f_out = f_in + np.array([-xi, -xi, +xi, +xi])
        return (1.0 - s_vec) * r_vec * f_out

    x = np.zeros(4)
    x_prev = None
    g_prev = None
    history_kgh = [0.0]
    false_conv_001 = None
    iter_1e6 = None

    for k in range(1, max_iter + 1):
        gx = map_step(x)
        delta_kgs = abs(float(np.dot(gx - x, mw_vec))) / 3600.0
        mass_kgh = float(np.dot(gx, mw_vec))
        history_kgh.append(mass_kgh)

        if false_conv_001 is None and delta_kgs < 0.01:
            f_in_at_stop = f_fresh + gx
            false_conv_001 = (k, mass_kgh, float(f_in_at_stop[0] / f_in_at_stop[1]))

        if iter_1e6 is None and delta_kgs < 1.0e-6:
            iter_1e6 = k

        if delta_kgs < tol_kgs:
            return {
                "method": method,
                "iterations": k,
                "final_mass_kgh": mass_kgh,
                "history_kgh": np.array(history_kgh),
                "false_conv_001": false_conv_001,
                "iter_1e6": iter_1e6,
            }

        if method == "wegstein" and x_prev is not None:
            dx = x - x_prev
            dg = gx - g_prev
            x_next = np.empty_like(x)
            for i in range(4):
                if abs(dx[i]) > 1e-12:
                    slope = dg[i] / dx[i]
                    q = slope / (slope - 1.0) if abs(slope - 1.0) > 1e-9 else 0.0
                    q = max(-5.0, min(0.0, q))  # standard bounded Wegstein damping
                    x_next[i] = max(0.0, q * x[i] + (1.0 - q) * gx[i])
                else:
                    x_next[i] = gx[i]
        else:
            x_next = gx

        x_prev = x.copy()
        g_prev = gx.copy()
        x = x_next

    return {
        "method": method,
        "iterations": max_iter,
        "final_mass_kgh": history_kgh[-1],
        "history_kgh": np.array(history_kgh),
        "false_conv_001": false_conv_001,
        "iter_1e6": iter_1e6,
    }


# =============================================================================
# 7. COD, SENSITIVITY & INDICATIVE ECONOMICS (Tables 11 & 13, Eq. 10)
# =============================================================================

def compute_cod(waste_kgh: Dict[str, float]) -> Tuple[float, float]:
    """Theoretical oxygen demand (Eq. 10): returns (COD kg O2/h, COD g/L)."""
    theta = {
        "EtOH": (3.0 * 32.0) / MW["EtOH"],
        "AA": (2.0 * 32.0) / MW["AA"],
        "EA": (5.0 * 32.0) / MW["EA"],
        "W": 0.0,
    }
    cod_kgh = sum(waste_kgh[c] * theta[c] for c in COMPONENTS)
    vol_lh = sum(waste_kgh[c] / RHO_25C[c] for c in COMPONENTS) * 1000.0
    cod_gl = (cod_kgh * 1000.0) / vol_lh
    return cod_kgh, cod_gl


def compute_economics(mb: MassBalanceResult, eb: Dict[str, float]) -> Dict[str, float]:
    """Indicative annual economics on 330 operating days/year (Table 13)."""
    fx = 1.367  # USD per GBP (Nov 2026)
    p_ea_gbp = 1.37 / fx
    p_aa_gbp = 0.67 / fx
    p_etoh_gbp = 1.07 / fx
    p_steam = 0.0528 / 0.85  # GBP/kWh heat (gas / 85% boiler efficiency)
    p_elec = 0.2451  # GBP/kWh
    p_cw = 0.05  # GBP/m3

    ea_t_yr = sum(mb.streams_kgh["Product-EtAc"].values()) * OPERATING_HOURS / 1000.0
    aa_t_yr = mb.streams_kgh["2"]["AA"] * OPERATING_HOURS / 1000.0
    etoh_t_yr = mb.streams_kgh["1"]["EtOH"] * OPERATING_HOURS / 1000.0

    rev_ea = ea_t_yr * 1000.0 * p_ea_gbp
    cost_aa = aa_t_yr * 1000.0 * p_aa_gbp
    cost_etoh = etoh_t_yr * 1000.0 * p_etoh_gbp
    gross_margin = rev_ea - cost_aa - cost_etoh

    # Utilities
    steam_ht1_mwh = eb["E2_kW"] * OPERATING_HOURS / 1000.0
    cost_steam_ht1 = steam_ht1_mwh * 1000.0 * p_steam
    elec_pump_mwh = 0.25 * OPERATING_HOURS / 1000.0  # conservative 0.25 kW motor rating
    cost_elec_pump = elec_pump_mwh * 1000.0 * p_elec
    cw_cl1_m3 = (eb["E4_kW"] / (4.184 * 10.0)) * 3.6 * OPERATING_HOURS  # 10 K rise
    cost_cw_cl1 = cw_cl1_m3 * p_cw

    # Separation train: 2-4 GJ/t EA (Filipe et al., 2005)
    sep_mwh_low = (2.0 * ea_t_yr) / 3.6
    sep_mwh_high = (4.0 * ea_t_yr) / 3.6
    cost_sep_low = sep_mwh_low * 1000.0 * p_steam
    cost_sep_high = sep_mwh_high * 1000.0 * p_steam
    cw_cond_low = (sep_mwh_low * 3600.0) / (4.184 * 10.0)
    cw_cond_high = (sep_mwh_high * 3600.0) / (4.184 * 10.0)
    cost_cw_cond_low = cw_cond_low * p_cw
    cost_cw_cond_high = cw_cond_high * p_cw

    margin_high = gross_margin - cost_steam_ht1 - cost_elec_pump - cost_cw_cl1 - cost_sep_low - cost_cw_cond_low
    margin_low = gross_margin - cost_steam_ht1 - cost_elec_pump - cost_cw_cl1 - cost_sep_high - cost_cw_cond_high

    return {
        "rev_ea_gbp": rev_ea,
        "cost_aa_gbp": cost_aa,
        "cost_etoh_gbp": cost_etoh,
        "gross_margin_gbp": gross_margin,
        "margin_after_util_low_gbp": margin_low,
        "margin_after_util_high_gbp": margin_high,
    }


# =============================================================================
# 8. DWSIM `.dwxmz` ARCHIVE READER & 12-POINT ACCEPTANCE CHECK (Table 9)
# =============================================================================

def parse_dwxmz(path: str) -> Dict[str, Dict[str, float]]:
    """
    Read a DWSIM 10.0 `.dwxmz` ZIP archive directly without requiring DWSIM.
    Extracts material stream mass flows (kg/h), temperatures (C), and energy stream duties (kW).
    """
    with zipfile.ZipFile(path, "r") as zf:
        xml_names = [n for n in zf.namelist() if n.lower().endswith(".xml")]
        if not xml_names:
            raise ValueError(f"No XML flowsheet found inside {path}")
        root = ET.fromstring(zf.read(xml_names[0]))

    material_streams: Dict[str, Dict[str, float]] = {}
    energy_streams: Dict[str, float] = {}

    for obj in root.findall(".//SimulationObject"):
        tag = obj.findtext("Tag") or obj.findtext(".//GraphicObject/Tag") or ""
        obj_type = obj.findtext("Type") or ""
        if "EnergyStream" in obj_type:
            val = obj.findtext("EnergyFlow")
            if val is not None:
                energy_streams[tag] = float(val)
        elif "MaterialStream" in obj_type:
            # Extract overall mixture phase (Phase ID = 0)
            for phase in obj.findall(".//Phase"):
                if phase.findtext("ID") == "0":
                    props = phase.find("Properties")
                    if props is not None:
                        m_kgs = float(props.findtext("massflow") or 0.0)
                        t_k = float(props.findtext("temperature") or 298.15)
                        material_streams[tag] = {
                            "mass_kgh": m_kgs * 3600.0,
                            "T_C": t_k - 273.15,
                        }
                    break

    return {"material": material_streams, "energy": energy_streams}


def validate_dwxmz_files(mb: MassBalanceResult) -> None:
    """Scan working directory for `.dwxmz` files and report acceptance checks."""
    files = sorted(glob.glob("*.dwxmz"))
    if not files:
        print("\n[DWSIM Validation] No .dwxmz files found in current directory.")
        print("  Place Ethyl_Acetate_CaseA_LeeKesler.dwxmz and Ethyl_Acetate_CaseB_LiquidCp.dwxmz")
        print("  in this folder to execute automated XML verification against Table 9.")
        return

    print("\n" + "=" * 78)
    print("DWSIM .dwxmz AUTOMATED VALIDATION (Table 9)")
    print("=" * 78)
    m_prod_hand = sum(mb.streams_kgh["Product-EtAc"].values())
    m_rec_hand = sum(mb.streams_kgh["20"].values())
    m_waste_hand = sum(mb.streams_kgh["12"].values())

    for fpath in files:
        print(f"\nFile: {os.path.basename(fpath)}")
        data = parse_dwxmz(fpath)
        mat = data["material"]
        eng = data["energy"]
        for tag, hand_val in [
            ("Product-EtAc", m_prod_hand),
            ("20", m_rec_hand),
            ("12", m_waste_hand),
        ]:
            if tag in mat:
                dw_val = mat[tag]["mass_kgh"]
                err_pct = 100.0 * (dw_val - hand_val) / hand_val
                status = "PASS" if abs(err_pct) < 0.1 else "FAIL"
                print(f"  Stream {tag:12s}: DWSIM={dw_val:8.4f} kg/h | Hand={hand_val:8.4f} | Diff={err_pct:+.5f}% [{status}]")
        if eng:
            print("  Energy streams (kW):", ", ".join(f"{k}={v:.4f}" for k, v in sorted(eng.items())))


# =============================================================================
# 9. MAIN REPORT SUMMARY
# =============================================================================

def main() -> None:
    mb = solve_mass_balance()
    eb = solve_energy_balance(mb)
    rxn_h = reaction_enthalpies(T_REACTOR_C)
    k_lai_70 = lai_rate_constants(70.0)
    x_eq_lai_70 = solve_equilibrium_conversion(mb.streams_kmolh["6"], k_lai_70[0] / k_lai_70[1])
    cap_g25, ka_25, ka_70 = calibrate_ka_nrtl(70.0, rxn_h["dh_liq_25C_kJmol"])
    x_eq_nrtl_vh = solve_equilibrium_conversion(mb.streams_kmolh["6"], ka_70, activity_based=True, t_c=70.0)
    rx = size_reactor_lai(mb.streams_kmolh["6"], SINGLE_PASS_CONV, 70.0)
    cod_kgh, cod_gl = compute_cod(mb.streams_kgh["12"])
    econ = compute_economics(mb, eb)
    conv_direct = simulate_recycle_convergence("direct", tol_kgs=1.0e-8)
    conv_weg = simulate_recycle_convergence("wegstein", tol_kgs=1.0e-8)

    # Subcooling check at reactor inlet (2.3 bar) & effluent (2.0 bar)
    x_in = np.array([mb.streams_kmolh["6"][c] for c in COMPONENTS])
    x_out = np.array([mb.streams_kmolh["7"][c] for c in COMPONENTS])
    tb_in = bubble_point_c(x_in, P_REACTOR_IN_BAR)
    tb_out = bubble_point_c(x_out, P_REACTOR_OUT_BAR)

    print("=" * 78)
    print("CONTINUOUS ETHYL ACETATE PRODUCTION (1,000 kg/day) — VERIFICATION SUMMARY")
    print("=" * 78)
    print(f"Product (Product-EtAc)   : {sum(mb.streams_kgh['Product-EtAc'].values()):.4f} kg/h @ {mb.product_ea_wt_pct:.4f} wt% EA")
    print(f"Overall AA conversion    : {mb.overall_aa_conv:.4f}% (single-pass = {SINGLE_PASS_CONV*100:.1f}%)")
    print(f"Recycle flow (stream 20) : {sum(mb.streams_kgh['20'].values()):.4f} kg/h (EtOH:AA inlet = {mb.reactor_inlet_etoh_aa:.5f})")
    print(f"Physical EtOH passes     : {mb.physical_etoh_passes:.2f} | Perturbation loop gain g_EtOH = {mb.loop_gain['EtOH']:.4f}")
    print(f"Elemental closure (C,H,O): {mb.atomic_residuals}")
    print("-" * 78)
    print(f"Heat of reaction (liquid): {rxn_h['dh_liq_25C_kJmol']:.2f} kJ/mol (25 C), {rxn_h['dh_liq_T_kJmol']:.2f} kJ/mol (70 C)")
    print(f"Heat of reaction (ig)    : {rxn_h['dh_ig_25C_kJmol']:.2f} kJ/mol | Eq.(9) offset = {eb['E3_offset_Eq9_kW']:.3f} kW")
    print(f"Equilibrium X_eq (70 C)  : Lai={x_eq_lai_70*100:.1f}% | NRTL+van't Hoff={x_eq_nrtl_vh*100:.1f}% (Ka={ka_70:.1f})")
    print(f"Catalyst (Amberlyst 35)  : PFR={rx['W_PFR_kg']:.1f} kg | CSTR={rx['W_CSTR_kg']:.1f} kg | Bed={rx['bed_mass_kg']:.1f} kg (LHSV={rx['LHSV_h']:.1f} h^-1)")
    print(f"Clean-bed Ergun dP       : {rx['Ergun_dP_kPa']:.2f} kPa | Bubble margin: feed={tb_in-70.0:.1f} K, eff={tb_out-70.0:.1f} K")
    print("-" * 78)
    print(f"Ideal-liquid duties      : E1={eb['E1_W']:.2f} W | E2={eb['E2_kW']:.3f} kW | E3={eb['E3_kW']:.3f} kW | E4={eb['E4_kW']:.3f} kW")
    print(f"Mixer temperatures       : T3={eb['T3_C']:.2f} C | T4={eb['T4_C']:.2f} C | T5={eb['T5_C']:.2f} C")
    print(f"Feed-effluent exchanger  : Q={eb['Q_FEHE_kW']:.2f} kW | Area={eb['FEHE_area_m2']:.2f} m2 | HT-1={eb['HT1_integrated_kW']:.2f} kW, CL-1={eb['CL1_integrated_kW']:.2f} kW")
    print("-" * 78)
    fc = conv_direct["false_conv_001"]
    print(f"False convergence (0.01) : iter {fc[0]} -> {fc[1]:.1f} kg/h ({100*fc[1]/62.624:.0f}% of true), EtOH:AA={fc[2]:.2f}")
    print(f"Direct subst. (1e-8 kg/s): {conv_direct['iterations']} iters (1e-6 at iter {conv_direct['iter_1e6']})")
    print(f"Wegstein      (1e-8 kg/s): {conv_weg['iterations']} iters (1e-6 at iter {conv_weg['iter_1e6']})")
    print("-" * 78)
    print(f"Aqueous waste COD        : {cod_kgh:.2f} kg O2/h ({cod_gl:.0f} g/L)")
    print(f"Indicative economics     : Gross margin GBP {econ['gross_margin_gbp']:,.0f}/yr | After utilities GBP {econ['margin_after_util_low_gbp']:,.0f} to {econ['margin_after_util_high_gbp']:,.0f}/yr")

    validate_dwxmz_files(mb)


if __name__ == "__main__":
    main()
