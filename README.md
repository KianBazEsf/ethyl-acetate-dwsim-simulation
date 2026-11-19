# Continuous Ethyl Acetate Production (1,000 kg/day)

**DWSIM 10.0 Process Simulation (NRTL), Independent Python Verification & Thermodynamic Sensitivity Study**  
*Kian Bazrafkan Esfahani — Department of Chemical Engineering, The University of Manchester*

---

## Project Overview

This repository contains the complete technical portfolio report, DWSIM stream exports, Excel validation workbook, and independent Python verification suite for a conceptual continuous esterification plant producing **1,000 kg/day (`41.67 kg/h`, `330 t/year`)** of ethyl acetate from ethanol and acetic acid over **Amberlyst 35** catalyst resin.

Rather than stopping at a converged simulator flowsheet, this study audits DWSIM's default solver and thermodynamic assumptions against closed-form analytical mass, energy, equilibrium, and kinetic calculations.

### Key Verified Outcomes

| Metric | Verified Result |
| :--- | :--- |
| **Production Basis** | `1,000 kg/day` EA (`41.67 kg/h`; `330 t/year` over 330 operating days) |
| **Product Stream (`Product-EtAc`)** | `41.68 kg/h` at `99.96 wt%` EA (meets sales composition limits by CS-1 target splits) |
| **Fresh Feeds** | Ethanol `22.25 kg/h`, Acetic Acid `29.93 kg/h` (`EtOH:AA = 0.97 mol/mol`) |
| **Recycle Flow (`Stream 20`)** | `62.62 kg/h` (`EtOH:AA = 2.00 mol/mol` at reactor inlet) |
| **Acetic Acid Conversion** | `65.0%` single-pass (imposed), `94.9%` overall |
| **Equilibrium Screening Estimate** | `73.5–81.2%` at `70 °C` and `EtOH:AA = 2.0` (model-dependent estimates) |
| **Reactor Sizing (Amberlyst 35)** | `15.4 kg` wet catalyst in PFR (`48.9 kg` in CSTR); design bed `0.20 m × 1.0 m` (`25.1 kg`, `LHSV = 4.4 h⁻¹`) |
| **Ideal-Liquid Duties (`E2 / E3 / E4`)** | `3.16 / 0.68 / 2.84 kW`; feed–effluent exchange recovers `2.32 kW` (`73%` heater / `82%` cooler savings) |
| **DWSIM Cross-Check** | Two solved runs (`Case A: Lee–Kesler` & `Case B: Liquid Cp + ΔHvap`); every component flow agrees within `0.009%`; `12/12` checks pass |
| **Aqueous Waste (`Stream 12`)** | `10.50 kg/h`, theoretical COD `2.58 kg O₂/h` (`~243 g/L`) |
| **Indicative Margin** | `£50–62k/year` after raw materials and utilities (learning-scale basis; excludes CAPEX, labour, and effluent treatment) |

---

## Three Simulation Pitfalls Exposed

1. **False Recycle Convergence & Snowball Sensitivity:**  
   Physically, `32.5%` of reactor-inlet ethanol reacts per pass so total ethanol averages **~3.0 physical reactor passes** (`1.455 / 0.483 = 3.01`). For excess-ethanol perturbations, however, acetic acid limits reaction extent so perturbations circulate with numerical loop gain `g_EtOH = 0.990` (amplification factor `≈ 99`). A Python reconstruction of DWSIM's default direct-substitution tolerance (`0.01 kg/s`) predicts false convergence after 2 iterations at `41%` of true recycle (`25.4 kg/h`, `EtOH:AA = 0.99`). Tightening tolerance to `1 × 10⁻⁸ kg/s` achieves `<0.00007%` external mass-balance residual.
2. **Enthalpy-Method Sensitivity (`Case A` vs. `Case B`):**  
   DWSIM's default **Lee–Kesler** option predicts heater (`E2`) and cooler (`E4`) duties `16–19%` above the ideal-liquid NIST/DIPPR benchmark alongside artificial mixer warming (`+4.9 K`) and non-zero separator enthalpy imbalances (`0.285` and `0.186 kW`). Switching to **Liquid `Cp` + Enthalpy of Vaporization** aligns `E2` and `E4` within `0.11%` and `0.03%` of the ideal-liquid benchmark (experimental excess-enthalpy validation still required) while raising reactor duty `E3` to `2.428 kW`—numerically reconciled by a `1.747 kW` ideal-gas vs. liquid reaction-enthalpy reference offset.
3. **NRTL Database Audit:**  
   Two of the six binary pairs are UNIFAC-estimated; the ternary azeotrope ester mole fraction differs from literature by `0.087` (`0.666` vs. `0.579`), a spurious acetic acid–water azeotrope appears, and the VLE-fitted ethyl acetate–water miscibility gap is too narrow for decanter design.

---

## Repository Contents

| File | Description |
| :--- | :--- |
| `Ethyl_Acetate_Report_BW.pdf` | Final 20-page pure black-and-white technical portfolio report |
| `Ethyl_Acetate_Report_BW.tex` | Full LaTeX source code (including native TikZ PFD & PGFPlots figures) |
| `Ethyl_Acetate_Calculations.py` | Independent NumPy/SciPy verification suite & automated `.dwxmz` XML validator |
| `make_figures.py` | Matplotlib script regenerating monochrome publication Figures 2–4 |
| `Ethyl_Acetate_Two_Model_Validation_Comparison.xlsx` | Formula-driven comparison of `Case A` (Lee–Kesler) vs. `Case B` (Liquid `Cp`) |
| `dwsim_exports/` | Raw DWSIM 10.0 stream & energy spreadsheet exports for both solved runs |
| `requirements.txt` | Python package dependencies |

### Adding Your Solved DWSIM Flowsheet Files Before Pushing to GitHub
Copy your two solved DWSIM 10.0 `.dwxmz` files from your Desktop into the repository root with the filenames referenced in Appendix A:
- `Final version for report(1).dwxmz` → **`Ethyl_Acetate_CaseA_LeeKesler.dwxmz`**
- `Final version for report liquidcp test.dwxmz` → **`Ethyl_Acetate_CaseB_LiquidCp.dwxmz`**

---

## How to Run the Verification Scripts

```bash
pip install -r requirements.txt
python Ethyl_Acetate_Calculations.py
python make_figures.py
```

When `Ethyl_Acetate_CaseA_LeeKesler.dwxmz` and `Ethyl_Acetate_CaseB_LiquidCp.dwxmz` are placed in the root directory, `Ethyl_Acetate_Calculations.py` automatically unpacks and validates them without requiring DWSIM to be installed.

---

## Limitations & Planned Next Stage

This flowsheet uses black-box component separators (`CS-1`, `CS-2`) and a fixed `65%` conversion reactor. Next engineering stage:
- Replace `CS-1` and `CS-2` with a rigorous **VLLE decanter + azeotropic distillation train** using regressed LLE NRTL parameters (`23–46 kW` true separation load).
- Replace the conversion reactor with a **kinetic PFR** using Lai et al. (2008) / Calvar et al. (2007) rate laws and add the feed–effluent heat exchanger directly into DWSIM.
