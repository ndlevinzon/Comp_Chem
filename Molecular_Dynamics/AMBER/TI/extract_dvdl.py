#!/usr/bin/env python3

from pathlib import Path
import re
import csv
import numpy as np

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR = Path(".")

LEGS = [
    "01_charge_off",
    "02_vdw_off",
]

OUTPUT_SUFFIX = "dvdl_summary.csv"

# Which output file to parse
# Usually v1 contains the alchemical transformation
OUTFILE_NAME = "prod_v1.out"

# Discard first fraction of DV/DL values
DISCARD_FRACTION = 0.20

# ============================================================
# REGEX
# ============================================================

lambda_pattern = re.compile(r"lambda_(\d+\.\d+)")
dvdl_pattern = re.compile(r"DV/DL\s*=\s*([-+0-9.Ee]+)")

# ============================================================
# FUNCTIONS
# ============================================================

def extract_lambda(dirname):
    """
    Extract lambda value from directory name:
    lambda_0.050 -> 0.050
    """
    m = lambda_pattern.search(dirname.name)

    if not m:
        raise ValueError(f"Could not parse lambda from {dirname}")

    return float(m.group(1))


def extract_dvdl_values(outfile):
    """
    Parse all DV/DL values from AMBER output file.
    """
    values = []

    with open(outfile, "r", errors="ignore") as f:

        for line in f:

            m = dvdl_pattern.search(line)

            if m:
                values.append(float(m.group(1)))

    return np.array(values, dtype=float)


def compute_stats(values, discard_fraction=0.20):

    if len(values) == 0:
        return np.nan, np.nan, 0

    discard_n = int(len(values) * discard_fraction)

    trimmed = values[discard_n:]

    if len(trimmed) == 0:
        return np.nan, np.nan, 0

    mean = np.mean(trimmed)

    if len(trimmed) > 1:
        sem = np.std(trimmed, ddof=1) / np.sqrt(len(trimmed))
    else:
        sem = np.nan

    return mean, sem, len(trimmed)


def process_leg(leg_dir):

    rows = []

    lambda_dirs = sorted(
        leg_dir.glob("lambda_*"),
        key=lambda x: extract_lambda(x)
    )

    if len(lambda_dirs) == 0:
        print(f"[WARNING] No lambda directories found in {leg_dir}")
        return rows

    for wdir in lambda_dirs:

        try:
            lam = extract_lambda(wdir)

        except Exception as e:

            print(f"[WARNING] Could not parse lambda in {wdir}: {e}")
            continue

        outfile = wdir / OUTFILE_NAME

        # ----------------------------------------------------
        # Missing output file
        # ----------------------------------------------------

        if not outfile.exists():

            print(f"[SKIP] Missing output: {outfile}")
            continue

        # ----------------------------------------------------
        # Empty output file
        # ----------------------------------------------------

        try:

            if outfile.stat().st_size == 0:

                print(f"[SKIP] Empty output file: {outfile}")
                continue

        except Exception as e:

            print(f"[SKIP] Could not stat file {outfile}: {e}")
            continue

        # ----------------------------------------------------
        # Parse DV/DL values
        # ----------------------------------------------------

        try:

            values = extract_dvdl_values(outfile)

        except Exception as e:

            print(f"[SKIP] Failed parsing DV/DL from {outfile}: {e}")
            continue

        # ----------------------------------------------------
        # No DV/DL values found
        # ----------------------------------------------------

        if len(values) == 0:

            print(f"[SKIP] No DV/DL values found in {outfile}")
            continue

        # ----------------------------------------------------
        # Compute statistics
        # ----------------------------------------------------

        try:

            mean, sem, n = compute_stats(
                values,
                discard_fraction=DISCARD_FRACTION
            )

        except Exception as e:

            print(f"[SKIP] Failed computing stats for {outfile}: {e}")
            continue

        print(
            f"{leg_dir.name:15s} "
            f"lambda={lam:6.3f} "
            f"N={n:6d} "
            f"mean={mean:12.6f} "
            f"SEM={sem:12.6f}"
        )

        rows.append({
            "lambda": lam,
            "n_values": n,
            "mean_dvdl": mean,
            "sem_dvdl": sem,
        })

    return rows


def write_csv(outfile, rows):

    fieldnames = [
        "lambda",
        "n_values",
        "mean_dvdl",
        "sem_dvdl",
    ]

    with open(outfile, "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=fieldnames)

        writer.writeheader()

        if len(rows) == 0:

            print(f"[WARNING] No valid windows found for {outfile}")

        for row in rows:
            writer.writerow(row)


# ============================================================
# MAIN
# ============================================================

print("========================================")
print("AMBER TI DV/DL extraction")
print("========================================")

for leg in LEGS:

    leg_dir = PROJECT_DIR / leg

    if not leg_dir.exists():

        print(f"[WARNING] Missing leg directory: {leg_dir}")
        continue

    print()
    print("----------------------------------------")
    print(f"Processing leg: {leg}")
    print("----------------------------------------")

    rows = process_leg(leg_dir)

    outfile = leg_dir / f"{leg}_{OUTPUT_SUFFIX}"

    write_csv(outfile, rows)

    print()
    print(f"Wrote: {outfile}")

print()
print("========================================")
print("Done")
print("========================================")
