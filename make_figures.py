#!/usr/bin/env python3
"""
Generate Pure Black-and-White Publication Figures (Figures 2, 3, and 4)
=======================================================================
Author: Kian Bazrafkan Esfahani
Department of Chemical Engineering, The University of Manchester

Imports the verified calculation routines from `Ethyl_Acetate_Calculations.py`
and saves 300-DPI monochrome figures in `./figures/`:
  - Figure2_Recycle_Convergence.png
  - Figure3_Enthalpy_Comparison.png
  - Figure4_Reactor_Tradeoffs.png

(Figure 1, the vector Process Flow Diagram, is embedded natively as pure B&W
TikZ in `Ethyl_Acetate_Report_BW.tex`.)

Dependencies: numpy, scipy, matplotlib
Usage:
    python make_figures.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import Ethyl_Acetate_Calculations as calc

OUT_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 9,
        "axes.edgecolor": "black",
        "axes.linewidth": 0.9,
        "text.color": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.dpi": 300,
    }
)


def plot_figure2_recycle_convergence() -> None:
    """Figure 2: Recycle flow vs iteration (direct substitution vs Wegstein)."""
    direct = calc.simulate_recycle_convergence("direct", tol_kgs=1.0e-8, max_iter=1000)
    weg = calc.simulate_recycle_convergence("wegstein", tol_kgs=1.0e-8, max_iter=1000)

    n_plot = 500
    h_dir = direct["history_kgh"]
    h_weg = weg["history_kgh"]
    x_dir = np.arange(min(len(h_dir), n_plot + 1))
    y_dir = h_dir[: len(x_dir)]
    x_weg = np.arange(n_plot + 1)
    y_weg = np.array([h_weg[i] if i < len(h_weg) else h_weg[-1] for i in x_weg])

    fig, ax = plt.subplots(figsize=(7.6, 3.2))
    ax.plot(x_dir, y_dir, color="black", lw=1.8, label="Direct substitution (DWSIM default)")
    ax.plot(x_weg, y_weg, color="black", lw=1.8, ls="--", label="Wegstein acceleration (DWSIM option)")
    ax.axhline(62.624, color="black", lw=0.9, ls=":")
    ax.text(495, 64.2, "true steady state: 62.6 kg/h, EtOH:AA = 2.00", ha="right", fontsize=8)

    # Mark false convergence at iteration 2 (default 0.01 kg/s tolerance)
    fc_iter, fc_mass, fc_ratio = direct["false_conv_001"]
    ax.scatter([fc_iter], [fc_mass], color="black", s=28, zorder=5)
    ax.annotate(
        f"False convergence at default tolerance (0.01 kg/s)\n"
        f"iter {fc_iter}: recycle {fc_mass:.1f} kg/h (41% of true), EtOH:AA = {fc_ratio:.2f}\n"
        f"equilibrium estimate 58%, below imposed 65%",
        xy=(fc_iter, fc_mass),
        xytext=(55, 7),
        fontsize=7.5,
        bbox=dict(boxstyle="square,pad=0.35", fc="white", ec="black", lw=0.8),
        arrowprops=dict(arrowstyle="->", color="black", lw=0.9),
    )

    # Mark 1e-6 kg/s checkpoints
    it_dir = direct["iter_1e6"]
    it_weg = weg["iter_1e6"]
    if it_weg and it_weg <= n_plot:
        ax.scatter([it_weg], [y_weg[it_weg]], color="black", s=24, zorder=5)
        ax.annotate(
            f"10$^{{-6}}$ kg/s at iter {it_weg}",
            xy=(it_weg, y_weg[it_weg]),
            xytext=(it_weg + 45, 42),
            fontsize=7.5,
            bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="black", lw=0.7),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
        )
    if it_dir and it_dir <= n_plot:
        ax.scatter([it_dir], [y_dir[it_dir]], color="black", s=24, zorder=5)
        ax.annotate(
            f"10$^{{-6}}$ kg/s at iter {it_dir}",
            xy=(it_dir, y_dir[it_dir]),
            xytext=(330, 42),
            fontsize=7.5,
            bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="black", lw=0.7),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
        )

    ax.set_xlim(0, n_plot)
    ax.set_ylim(0, 72)
    ax.set_xlabel("Recycle iteration")
    ax.set_ylabel("Recycle flow, stream 20 (kg/h)")
    ax.grid(True, ls=":", color="0.75", lw=0.6)
    ax.legend(loc="upper left", frameon=True, edgecolor="black", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "Figure2_Recycle_Convergence.png"))
    plt.close(fig)


def plot_figure3_enthalpy_comparison() -> None:
    """Figure 3: Pure B&W bar chart comparing hand calc, Case A (Lee-Kesler), and Case B (Liquid Cp)."""
    labels = ["E2\nheater HT-1", "E3\nreactor RCONV-1", "E4\ncooler CL-1", "E5 + E6\nseparators CS-1, CS-2"]
    hand = [3.161, 0.680, 2.842, 0.000]
    case_a = [3.675, 0.750, 3.387, 0.470]
    case_b_base = [3.165, 0.681, 2.842, 0.000]
    case_b_offset = [0.000, 1.747, 0.000, 0.000]

    x = np.arange(len(labels))
    w = 0.24

    fig, ax = plt.subplots(figsize=(7.6, 3.5))
    ax.bar(x - w, hand, w, color="0.20", edgecolor="black", lw=0.8, label="Hand calculation (NIST/DIPPR ideal-liquid)")
    ax.bar(x, case_a, w, color="0.65", edgecolor="black", lw=0.8, label="DWSIM Case A: Lee–Kesler (default)")
    ax.bar(x + w, case_b_base, w, color="white", edgecolor="black", lw=0.8, label="DWSIM Case B: Liquid $C_p$ + $\\Delta H_{\\mathrm{vap}}$")
    ax.bar(
        x + w,
        case_b_offset,
        w,
        bottom=case_b_base,
        color="white",
        edgecolor="black",
        hatch="///",
        lw=0.8,
        label="Case B reactor: ideal-gas $\\Delta_{\\mathrm{r}}H$ offset (1.75 kW)",
    )

    for i in range(len(labels)):
        ax.text(x[i] - w, hand[i] + 0.07, f"{hand[i]:.2f}", ha="center", va="bottom", fontsize=7.5)
        ax.text(x[i], case_a[i] + 0.07, f"{case_a[i]:.2f}", ha="center", va="bottom", fontsize=7.5)
        tot_b = case_b_base[i] + case_b_offset[i]
        ax.text(x[i] + w, tot_b + 0.07, f"{tot_b:.2f}", ha="center", va="bottom", fontsize=7.5)

    ax.annotate(
        "1.75 kW offset:\n2.43 − 1.75 = 0.68 kW",
        xy=(1 + w, 2.55),
        xytext=(1.0, 3.35),
        ha="center",
        fontsize=7.5,
        bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="black", lw=0.7),
        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 4.35)
    ax.set_ylabel("Duty (kW)")
    ax.grid(axis="y", ls=":", color="0.75", lw=0.6)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.24), ncol=2, frameon=True, edgecolor="black", fontsize=7.8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "Figure3_Enthalpy_Comparison.png"), bbox_inches="tight")
    plt.close(fig)


def plot_figure4_reactor_tradeoffs() -> None:
    """Figure 4: 2x2 panel of reactor design trade-offs at 70 C (Lai et al. kinetics)."""
    k1, km1 = calc.lai_rate_constants(70.0)
    k_lai = k1 / km1

    # Sweep 1: EtOH:AA inlet ratio at X = 0.65 (fixed EA production = 41.67 kg/h)
    ratios = np.linspace(1.5, 3.0, 31)
    xeq_ratio = []
    rec_ratio = []
    mb_base = calc.solve_mass_balance()
    xi_fixed = mb_base.extent_kmolh
    f_aa_in_fixed = xi_fixed / 0.65
    g = mb_base.loop_gain
    f_ea_in = (g["EA"] * xi_fixed) / (1.0 - g["EA"])
    f_w_in = (g["W"] * xi_fixed) / (1.0 - g["W"])

    for r_val in ratios:
        f_etoh_in = r_val * f_aa_in_fixed
        f_in = {"EtOH": f_etoh_in, "AA": f_aa_in_fixed, "EA": f_ea_in, "W": f_w_in}
        xeq = calc.solve_equilibrium_conversion(f_in, k_lai)
        f_out = {c: f_in[c] + calc.NU[c] * xi_fixed for c in calc.COMPONENTS}
        m_rec = sum(g[c] * f_out[c] * calc.MW[c] for c in calc.COMPONENTS)
        xeq_ratio.append(100.0 * xeq)
        rec_ratio.append(m_rec)

    # Sweep 2: single-pass AA conversion at EtOH:AA = 2.0
    convs = np.linspace(0.40, 0.74, 35)
    w_pfr_list = []
    w_cstr_list = []
    rec_conv = []
    for x_val in convs:
        f_aa_in = xi_fixed / x_val
        f_in = {"EtOH": 2.0 * f_aa_in, "AA": f_aa_in, "EA": f_ea_in, "W": f_w_in}
        rx = calc.size_reactor_lai(f_in, x_val, 70.0)
        f_out = {c: f_in[c] + calc.NU[c] * xi_fixed for c in calc.COMPONENTS}
        m_rec = sum(g[c] * f_out[c] * calc.MW[c] for c in calc.COMPONENTS)
        w_pfr_list.append(rx["W_PFR_kg"])
        w_cstr_list.append(rx["W_CSTR_kg"])
        rec_conv.append(m_rec)

    fig, axs = plt.subplots(2, 2, figsize=(7.8, 5.2))

    # (a) Xeq vs ratio
    axs[0, 0].plot(ratios, xeq_ratio, color="black", lw=1.6)
    axs[0, 0].axhline(65.0, color="black", lw=0.9, ls="--")
    axs[0, 0].scatter([2.0], [76.5], color="black", s=26, zorder=5)
    axs[0, 0].set_title("(a) Equilibrium conversion vs inlet ratio", fontsize=9, fontweight="bold")
    axs[0, 0].set_xlabel("EtOH:AA molar ratio at reactor inlet")
    axs[0, 0].set_ylabel("$X_{\\mathrm{eq}}$ at 70 °C (%)")
    axs[0, 0].set_ylim(60, 88)
    axs[0, 0].grid(True, ls=":", color="0.75", lw=0.6)

    # (b) Recycle vs ratio
    axs[0, 1].plot(ratios, rec_ratio, color="black", lw=1.6)
    axs[0, 1].scatter([2.0], [62.62], color="black", s=26, zorder=5)
    axs[0, 1].set_title("(b) Recycle flow vs inlet ratio", fontsize=9, fontweight="bold")
    axs[0, 1].set_xlabel("EtOH:AA molar ratio at reactor inlet")
    axs[0, 1].set_ylabel("Recycle, stream 20 (kg/h)")
    axs[0, 1].set_ylim(0, 105)
    axs[0, 1].grid(True, ls=":", color="0.75", lw=0.6)

    # (c) Catalyst vs conversion
    axs[1, 0].plot(convs * 100.0, w_pfr_list, color="black", lw=1.6, label="plug-flow bed")
    axs[1, 0].plot(convs * 100.0, w_cstr_list, color="black", lw=1.4, ls="--", label="stirred tank")
    axs[1, 0].axvline(76.5, color="black", lw=0.9, ls=":")
    axs[1, 0].scatter([65.0, 65.0], [15.4, 48.9], color="black", s=26, zorder=5)
    axs[1, 0].text(74.8, 56, "equilibrium\n76.5%", ha="right", fontsize=7.5)
    axs[1, 0].set_title("(c) Catalyst vs single-pass conversion", fontsize=9, fontweight="bold")
    axs[1, 0].set_xlabel("Single-pass AA conversion (%)")
    axs[1, 0].set_ylabel("Amberlyst 35, wet (kg)")
    axs[1, 0].set_xlim(39, 78)
    axs[1, 0].set_ylim(0, 105)
    axs[1, 0].legend(loc="upper left", frameon=True, edgecolor="black", fontsize=7.5)
    axs[1, 0].grid(True, ls=":", color="0.75", lw=0.6)

    # (d) Recycle vs conversion
    axs[1, 1].plot(convs * 100.0, rec_conv, color="black", lw=1.6)
    axs[1, 1].scatter([65.0], [62.62], color="black", s=26, zorder=5)
    axs[1, 1].set_title("(d) Recycle flow vs single-pass conversion", fontsize=9, fontweight="bold")
    axs[1, 1].set_xlabel("Single-pass AA conversion (%)")
    axs[1, 1].set_ylabel("Recycle, stream 20 (kg/h)")
    axs[1, 1].set_xlim(39, 78)
    axs[1, 1].set_ylim(0, 130)
    axs[1, 1].grid(True, ls=":", color="0.75", lw=0.6)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "Figure4_Reactor_Tradeoffs.png"))
    plt.close(fig)


def main() -> None:
    plot_figure2_recycle_convergence()
    plot_figure3_enthalpy_comparison()
    plot_figure4_reactor_tradeoffs()
    print(f"Saved monochrome publication figures to ./{OUT_DIR}/")


if __name__ == "__main__":
    main()
