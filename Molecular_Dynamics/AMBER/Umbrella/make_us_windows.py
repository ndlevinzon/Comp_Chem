#!/usr/bin/env python3
import argparse
from pathlib import Path

# -----------------------
# USER SETTINGS (edit me)
# -----------------------
# Atom indices for chi dihedral (AMBER atom numbering, 1-based)
# Purine:  O4' C1' N9 C4
# Pyrimidine: O4' C1' N1 C2
IAT = (5772, 5773, 5747, 5761)

START_DEG = 0
END_DEG_EXCLUSIVE = 360   # recommended: 0..350 (avoid duplicate 360)
STEP_DEG = 10

# Production umbrella width + strength
HALF_WIDTH_DEG_PROD = 10.0
RK_PROD = 100.0  # kcal/mol/rad^2

# Seeding (stronger, narrower) to pull chi to target
HALF_WIDTH_DEG_SEED = 2.0
RK_SEED = 500.0  # kcal/mol/rad^2

# Template filenames expected in input-dir
PROD_TEMPLATE_NAME = "us_prod.mdin.template"
SEED_MIN_TEMPLATE_NAME = "us_seed_min.mdin.template"
SEED_MD_TEMPLATE_NAME  = "us_seed_md.mdin.template"
# -----------------------


def tag_from_deg(deg: int) -> str:
    return f"{deg:03d}" if deg >= 0 else f"m{abs(deg):03d}"


def write_disang(target_deg: float, half_width_deg: float, rk: float, iat, outpath: Path):
    r1 = target_deg - half_width_deg
    r2 = target_deg
    r3 = target_deg
    r4 = target_deg + half_width_deg

    text = f"""# Dihedral umbrella restraint at {target_deg:.1f} degrees
&rst
  iat={iat[0]},{iat[1]},{iat[2]},{iat[3]},
  r1={r1:.3f}, r2={r2:.3f}, r3={r3:.3f}, r4={r4:.3f},
  rk2={rk:.3f}, rk3={rk:.3f},
/
"""
    outpath.write_text(text)


def fill_template(template_path: Path, target_deg: int, target_tag: str, outpath: Path):
    txt = template_path.read_text()
    txt = txt.replace("{TARGET_DEG}", str(target_deg))
    txt = txt.replace("{TARGET_TAG}", target_tag)
    outpath.write_text(txt)


def main():
    ap = argparse.ArgumentParser(
        description="Generate chi umbrella window directories with seed_min/seed_md/prod mdins + DISANG files."
    )
    ap.add_argument("--input-dir", default=".", help="Directory containing mdin templates (default: .)")
    ap.add_argument("--output-dir", default="windows", help="Output root (default: ./windows)")
    args = ap.parse_args()

    input_dir = Path(args.input_dir).resolve()
    out_root = Path(args.output_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    prod_template = input_dir / PROD_TEMPLATE_NAME
    seedmin_template = input_dir / SEED_MIN_TEMPLATE_NAME
    seedmd_template = input_dir / SEED_MD_TEMPLATE_NAME

    if not prod_template.exists():
        raise FileNotFoundError(f"Missing required template: {prod_template}")
    if not seedmin_template.exists():
        raise FileNotFoundError(f"Missing required template: {seedmin_template}")
    if not seedmd_template.exists():
        raise FileNotFoundError(f"Missing required template: {seedmd_template}")

    degrees = list(range(START_DEG, END_DEG_EXCLUSIVE, STEP_DEG))

    for deg in degrees:
        tag = tag_from_deg(deg)
        wdir = out_root / f"chi_{tag}"
        wdir.mkdir(parents=True, exist_ok=True)

        # DISANG files: strong for seeding, normal for production
        write_disang(deg, HALF_WIDTH_DEG_SEED, RK_SEED, IAT, wdir / f"disang_seed_{tag}.RST")
        write_disang(deg, HALF_WIDTH_DEG_PROD, RK_PROD, IAT, wdir / f"disang_prod_{tag}.RST")

        # mdins from templates
        fill_template(seedmin_template, deg, tag, wdir / f"seed_min_{tag}.mdin")
        fill_template(seedmd_template,  deg, tag, wdir / f"seed_md_{tag}.mdin")
        fill_template(prod_template,    deg, tag, wdir / f"prod_{tag}.mdin")

    print(f"[DONE] Wrote {len(degrees)} windows in: {out_root}")


if __name__ == "__main__":
    main()
