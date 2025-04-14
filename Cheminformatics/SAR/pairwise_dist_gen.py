import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, is_aa
import sys


def residue_distance_matrix(pdb_path, chain_id=None, output_csv=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)

    # Extract all CA atoms from protein residues
    ca_atoms = []
    res_ids = []

    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for residue in chain:
                if is_aa(residue, standard=True) and 'CA' in residue:
                    ca_atoms.append(residue['CA'].coord)
                    res_ids.append(f"{chain.id}_{residue.get_id()[1]}")
            break  # Only use first model
        break

    coords = np.array(ca_atoms)
    dist_matrix = np.linalg.norm(coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=-1)

    df = pd.DataFrame(dist_matrix, index=res_ids, columns=res_ids)

    if output_csv:
        df.to_csv("residue_distance_matrix.csv")

    np.save("residue_distance_matrix.npy", dist_matrix)
    return df

# Example usage:
df = residue_distance_matrix("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake.pdb", chain_id="A", output_csv=True)
