#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR = Path(".")

CHARGE_CSV = PROJECT_DIR / "01_charge_off" / "01_charge_off_dvdl_summary.csv"
VDW_CSV    = PROJECT_DIR / "02_vdw_off" / "02_vdw_off_dvdl_summary.csv"

# Optional charged-ion correction
CHARGE_CORRECTION_KCAL_MOL = 0.0

# ============================================================
# FUNCTIONS
# ============================================================

def load_ti_csv(csv_file):

    df = pd.read_csv(csv_file)

    required = {"lambda", "mean_dvdl", "sem_dvdl"}

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"{csv_file} missing columns: {missing}")

    df = df.dropna(subset=["lambda", "mean_dvdl", "sem_dvdl"])

    df = df.sort_values("lambda")

    return df


def integrate_leg(df):

    lambdas = df["lambda"].to_numpy(dtype=float)
    dvdl = df["mean_dvdl"].to_numpy(dtype=float)

    dg = np.trapz(dvdl, lambdas)

    return dg


def estimate_leg_uncertainty(df):
    """
    Approximate propagated uncertainty from SEM values.

    This assumes independent lambda windows.
    """

    lambdas = df["lambda"].to_numpy(dtype=float)
    sems = df["sem_dvdl"].to_numpy(dtype=float)

    variances = []

    for i in range(len(lambdas) - 1):

        dlam = lambdas[i + 1] - lambdas[i]

        # trapezoid uncertainty estimate
        local_var = (
            (0.5 * dlam * sems[i])**2 +
            (0.5 * dlam * sems[i + 1])**2
        )

        variances.append(local_var)

    total_std = np.sqrt(np.sum(variances))

    return total_std


def summarize_leg(name, csv_file):

    df = load_ti_csv(csv_file)

    dg = integrate_leg(df)

    std = estimate_leg_uncertainty(df)

    print(name)
    print("-" * len(name))
    print(f"CSV file: {csv_file}")
    print(f"Lambda range: {df['lambda'].min():.3f} -> {df['lambda'].max():.3f}")
    print(f"Number of windows: {len(df)}")
    print(f"ΔG decoupling = {dg: .6f} kcal/mol")
    print(f"STD estimate  = {std: .6f} kcal/mol")
    print(f"±2 STD        = {2.0 * std: .6f} kcal/mol")
    print()

    return dg, std, df


# ============================================================
# MAIN
# ============================================================

dg_charge, std_charge, df_charge = summarize_leg(
    "Charge-removal leg",
    CHARGE_CSV
)

dg_vdw, std_vdw, df_vdw = summarize_leg(
    "vdW-removal leg",
    VDW_CSV
)

# ------------------------------------------------------------
# Total decoupling free energy
# ------------------------------------------------------------

dg_decouple_raw = dg_charge + dg_vdw

# Propagate uncertainties assuming independent legs
std_total = np.sqrt(std_charge**2 + std_vdw**2)

# Hydration free energy is reverse of decoupling
dg_hyd_raw = -dg_decouple_raw

dg_hyd_corrected = dg_hyd_raw + CHARGE_CORRECTION_KCAL_MOL

# ============================================================
# OUTPUT
# ============================================================

print("=" * 70)
print("TI Hydration Free Energy Summary")
print("=" * 70)

print(f"ΔG_charge decoupling   = {dg_charge: .6f} kcal/mol")
print(f"ΔG_vdW decoupling      = {dg_vdw: .6f} kcal/mol")
print()

print(f"Raw ΔG_decouple        = {dg_decouple_raw: .6f} kcal/mol")
print(f"Raw ΔG_hydration       = {dg_hyd_raw: .6f} kcal/mol")
print()

print(f"Charge correction      = {CHARGE_CORRECTION_KCAL_MOL: .6f} kcal/mol")
print(f"Corrected hydration ΔG = {dg_hyd_corrected: .6f} kcal/mol")
print()

print(f"STD estimate           = {std_total: .6f} kcal/mol")
print(f"±2 STD                 = {2.0 * std_total: .6f} kcal/mol")

print("=" * 70)
