#!/usr/bin/env python3
import os
import re
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

# --- Protonation enumeration (Dimorphite-DL) ---
# pip install dimorphite_dl
try:
    # Most common import style
    from dimorphite_dl import DimorphiteDL
except Exception:
    DimorphiteDL = None

# --- PDBQT writing (Meeko) ---
from meeko import MoleculePreparation


@dataclass
class PrepConfig:
    ph: float = 7.4
    ph_tolerance: float = 1.0
    max_protomers: int = 12
    max_tautomers: int = 8
    num_confs: int = 30
    max_confs_written: int = 10
    prune_rms: float = 0.5
    minimize: bool = True
    seed: int = 2025


# -----------------------------
# Utilities
# -----------------------------
def safe_name(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s[:200] if len(s) > 200 else s


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    if smiles is None or str(smiles).strip() == "":
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    # sanitize explicitly to fail early if needed
    Chem.SanitizeMol(m)
    return m


def enumerate_protonation_states(smiles: str, cfg: PrepConfig) -> List[str]:
    """
    Returns a list of protonated SMILES candidates.
    Uses Dimorphite-DL if available; otherwise returns the original SMILES.
    """
    if DimorphiteDL is None:
        return [smiles]

    # DimorphiteDL can be used in different ways depending on version.
    # This pattern works for many installations:
    dim = DimorphiteDL(
        min_ph=cfg.ph - cfg.ph_tolerance,
        max_ph=cfg.ph + cfg.ph_tolerance,
        pka_precision=1.0,
        max_variants=cfg.max_protomers,
        label_states=False,
    )
    out = dim.protonate(smiles)

    # Ensure list + uniqueness
    if isinstance(out, str):
        out = [out]
    out = [s for s in out if s and isinstance(s, str)]
    # Canonicalize to dedupe
    uniq = []
    seen = set()
    for s in out:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        can = Chem.MolToSmiles(m, isomericSmiles=True)
        if can not in seen:
            seen.add(can)
            uniq.append(can)
        if len(uniq) >= cfg.max_protomers:
            break
    return uniq if uniq else [smiles]


def enumerate_tautomers_rdkit(m: Chem.Mol, cfg: PrepConfig) -> List[Chem.Mol]:
    """
    RDKit tautomer enumeration is available via rdMolStandardize in many builds.
    If not available, fall back to returning [m].
    """
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        te = rdMolStandardize.TautomerEnumerator()
        tauts = te.Enumerate(m)
        # Deduplicate by canonical SMILES
        uniq = []
        seen = set()
        for tm in tauts:
            can = Chem.MolToSmiles(tm, isomericSmiles=True)
            if can not in seen:
                seen.add(can)
                uniq.append(tm)
            if len(uniq) >= cfg.max_tautomers:
                break
        return uniq if uniq else [m]
    except Exception:
        return [m]


def embed_and_minimize(m: Chem.Mol, cfg: PrepConfig) -> Tuple[Chem.Mol, List[int], List[float]]:
    """
    Add H, generate conformers (ETKDGv3), optionally MMFF minimize.
    Returns: (mol_with_confs, conf_ids, conf_energies)
    """
    mol = Chem.AddHs(m)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(cfg.seed)
    params.pruneRmsThresh = float(cfg.prune_rms)

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(cfg.num_confs), params=params))
    if not conf_ids:
        raise RuntimeError("3D embedding failed (no conformers generated).")

    energies = []
    if cfg.minimize:
        # MMFF preferred; fallback to UFF
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if mp is not None:
            for cid in conf_ids:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
                ff.Minimize(maxIts=200)
                energies.append(float(ff.CalcEnergy()))
        else:
            for cid in conf_ids:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                ff.Minimize(maxIts=200)
                energies.append(float(ff.CalcEnergy()))
    else:
        energies = [0.0 for _ in conf_ids]

    return mol, conf_ids, energies


def select_top_conformers(conf_ids: List[int], energies: List[float], max_keep: int) -> List[int]:
    order = np.argsort(np.asarray(energies, dtype=float))
    keep = [conf_ids[i] for i in order[:max_keep]]
    return keep


def write_pdbqt_meeko(molH: Chem.Mol, conf_id: int, out_pdbqt: str) -> None:
    # --- build a single-conformer mol (most compatible across Meeko versions) ---
    m = Chem.Mol(molH)
    conf = molH.GetConformer(conf_id)
    m.RemoveAllConformers()
    m.AddConformer(Chem.Conformer(conf), assignId=True)

    prep = MoleculePreparation()

    setups = prep.prepare(m)
    if isinstance(setups, list):
        if len(setups) == 0:
            raise RuntimeError("Meeko prepare() returned empty list.")
        setup = setups[0]
    else:
        setup = setups

    # --- Export PDBQT depending on Meeko version ---
    # Newer Meeko: use PDBQTWriterLegacy (setup has no write_pdbqt_string)
    try:
        from meeko import PDBQTWriterLegacy
        pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setup)
        if not ok:
            raise RuntimeError(f"PDBQTWriterLegacy failed: {err}")
    except ImportError:
        # Older Meeko (rare): setup might have the method
        if hasattr(setup, "write_pdbqt_string"):
            pdbqt_string = setup.write_pdbqt_string()
        else:
            raise RuntimeError(
                "Cannot export PDBQT: Meeko has no PDBQTWriterLegacy and setup has no write_pdbqt_string(). "
                "Check your meeko version."
            )

    with open(out_pdbqt, "w") as f:
        f.write(pdbqt_string)




# -----------------------------
# Main pipeline
# -----------------------------
def process_one_ligand(lig_id: str, smiles: str, outdir: str, cfg: PrepConfig) -> List[dict]:
    """
    Returns a list of result dicts (one per written pdbqt).
    """
    lig_id_safe = safe_name(lig_id)
    lig_dir = os.path.join(outdir, lig_id_safe)
    os.makedirs(lig_dir, exist_ok=True)

    results = []
    prot_smiles_list = enumerate_protonation_states(smiles, cfg)

    state_counter = 0
    for prot_smiles in prot_smiles_list:
        try:
            base_m = mol_from_smiles(prot_smiles)
            if base_m is None:
                continue
        except Exception:
            continue

        # Optional tautomer enumeration (often helpful for heterocycles)
        tautomers = enumerate_tautomers_rdkit(base_m, cfg)

        for t_i, tm in enumerate(tautomers):
            state_counter += 1
            if state_counter > (cfg.max_protomers * cfg.max_tautomers):
                break

            # Generate conformers
            try:
                mol3d, conf_ids, energies = embed_and_minimize(tm, cfg)
            except Exception as e:
                results.append({
                    "Identifier": lig_id,
                    "input_smiles": smiles,
                    "prot_smiles": prot_smiles,
                    "tautomer_index": t_i,
                    "status": f"FAILED_3D: {e}",
                    "pdbqt": ""
                })
                continue

            keep_conf_ids = select_top_conformers(conf_ids, energies, cfg.max_confs_written)

            for k, cid in enumerate(keep_conf_ids):
                out_name = f"{lig_id_safe}__p{state_counter:02d}__t{t_i:02d}__c{k:02d}.pdbqt"
                out_path = os.path.join(lig_dir, out_name)

                try:
                    write_pdbqt_meeko(mol3d, cid, out_path)
                    results.append({
                        "Identifier": lig_id,
                        "input_smiles": smiles,
                        "prot_smiles": prot_smiles,
                        "tautomer_index": t_i,
                        "conf_rank": k,
                        "status": "OK",
                        "pdbqt": out_path
                    })
                except Exception as e:
                    results.append({
                        "Identifier": lig_id,
                        "input_smiles": smiles,
                        "prot_smiles": prot_smiles,
                        "tautomer_index": t_i,
                        "conf_rank": k,
                        "status": f"FAILED_PDBQT: {e}",
                        "pdbqt": ""
                    })

    return results


def main():
    ap = argparse.ArgumentParser(description="CSV (Identifier,Canonical.SMILES) -> protonated/tautomerized 3D conformers -> PDBQT (AutoDock-GPU ready)")
    ap.add_argument("--csv", required=True, help="Input CSV with headers Identifier,Canonical.SMILES")
    ap.add_argument("--outdir", default="pdbqt_out", help="Output directory")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--ph-tol", type=float, default=1.0)
    ap.add_argument("--max-protomers", type=int, default=12)
    ap.add_argument("--max-tautomers", type=int, default=8)
    ap.add_argument("--num-confs", type=int, default=30)
    ap.add_argument("--max-confs-written", type=int, default=10)
    ap.add_argument("--prune-rms", type=float, default=0.5)
    ap.add_argument("--no-minimize", action="store_true")
    ap.add_argument("--seed", type=int, default=2025)
    ap.add_argument("--limit", type=int, default=None, help="Only process first N ligands (debug)")
    args = ap.parse_args()

    cfg = PrepConfig(
        ph=args.ph,
        ph_tolerance=args.ph_tol,
        max_protomers=args.max_protomers,
        max_tautomers=args.max_tautomers,
        num_confs=args.num_confs,
        max_confs_written=args.max_confs_written,
        prune_rms=args.prune_rms,
        minimize=not args.no_minimize,
        seed=args.seed,
    )

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    if "Identifier" not in df.columns or "Canonical.SMILES" not in df.columns:
        raise ValueError("CSV must contain headers: Identifier,Canonical.SMILES")

    if args.limit is not None:
        df = df.head(int(args.limit))

    all_results = []
    for i, row in df.iterrows():
        lig_id = str(row["Identifier"])
        smi = str(row["Canonical.SMILES"])
        if not smi or smi.lower() == "nan":
            continue
        try:
            all_results.extend(process_one_ligand(lig_id, smi, args.outdir, cfg))
        except Exception as e:
            all_results.append({
                "Identifier": lig_id,
                "input_smiles": smi,
                "status": f"FAILED_LIGAND: {e}",
                "pdbqt": ""
            })

        if (i + 1) % 50 == 0:
            print(f"[INFO] Processed {i+1}/{len(df)} ligands...")

    out_table = os.path.join(args.outdir, "pdbqt_manifest.csv")
    pd.DataFrame(all_results).to_csv(out_table, index=False)
    print(f"[DONE] Wrote manifest: {out_table}")
    print(f"[DONE] PDBQT output root: {args.outdir}")


if __name__ == "__main__":
    main()
