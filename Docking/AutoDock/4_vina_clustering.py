import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import os

def extract_smiles_and_coords_from_pdb(pdb_file_path):
    """
    Extracts ligand names, SMILES codes, dock scores, distances, and coordinates from the PDB file.

    Args:
        pdb_file_path (str): Path to the PDB file.

    Returns:
        list: A list of tuples containing ligand names, SMILES codes, dock scores, distances, and coordinates.
    """
    data = []
    with open(pdb_file_path, 'r') as file:
        lines = file.readlines()

    current_name = None
    current_smiles = None
    current_score = None
    current_distance = None
    current_coords = []

    for line in lines:
        if "REMARK Name =" in line:
            match = re.search(r'REMARK Name = (\S+)', line)
            if match:
                current_name = match.group(1)
        elif "REMARK SMILES =" in line:
            match = re.search(r'REMARK SMILES = (\S+)', line)
            if match:
                current_smiles = match.group(1)
        elif "REMARK Vina score:" in line:
            match = re.search(r'REMARK Vina score: (-?\d+\.\d+)', line)
            if match:
                current_score = match.group(1)
        elif "REMARK Distance from crystal ligand centroid:" in line:
            match = re.search(r'REMARK Distance from crystal ligand centroid: (\d+\.\d+)', line)
            if match:
                current_distance = match.group(1)
        elif line.startswith("HETATM") or line.startswith("ATOM"):
            current_coords.append(line.strip())
        elif "ENDMDL" in line and current_name and current_smiles and current_score and current_distance:
            data.append((current_name, current_smiles, current_score, current_distance, current_coords))
            current_name = None
            current_smiles = None
            current_score = None
            current_distance = None
            current_coords = []

    return data

def generate_pharmacophore_features(smiles_list):
    """
    Generates pharmacophore features from a list of SMILES codes.

    Args:
        smiles_list (list): List of SMILES codes.

    Returns:
        list: List of pharmacophore fingerprints.
    """
    factory = Gobbi_Pharm2D.factory
    pharm_fps = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = Generate.Gen2DFingerprint(mol, factory)
            pharm_fps.append(fp)
    
    return pharm_fps

def convert_fingerprints_to_vectors(pharm_fps):
    """
    Converts pharmacophore fingerprints to binary vectors.

    Args:
        pharm_fps (list): List of pharmacophore fingerprints.

    Returns:
        list: List of binary vectors.
    """
    vectors = []
    for fp in pharm_fps:
        arr = np.zeros((1,))
        explicit_fp = DataStructs.ExplicitBitVect(len(fp))
        for bit in fp.GetOnBits():
            explicit_fp.SetBit(bit)
        DataStructs.ConvertToNumpyArray(explicit_fp, arr)
        vectors.append(arr)
    
    return np.array(vectors)

def cluster_ligands(fingerprint_vectors):
    """
    Clusters ligands based on pharmacophore feature vectors.

    Args:
        fingerprint_vectors (ndarray): Binary vectors of pharmacophore features.

    Returns:
        ndarray: Cluster labels for each ligand.
    """
    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=50).fit(fingerprint_vectors)
    return clustering.labels_

def save_clustered_ligands(data, cluster_labels, output_dir):
    """
    Saves clustered ligands into PDB files.

    Args:
        data (list): List of tuples containing ligand names, SMILES codes, dock scores, distances, and coordinates.
        cluster_labels (ndarray): Cluster labels for each ligand.
        output_dir (str): Directory to save the output files.
    """
    clusters = {}
    for (name, smiles, score, distance, coords), label in zip(data, cluster_labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append((name, smiles, score, distance, coords))
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for label, ligands in clusters.items():
        with open(os.path.join(output_dir, f'cluster_{label}.pdb'), 'w') as f:
            for ligand in ligands:
                name, smiles, score, distance, coords = ligand
                f.write("MODEL\n")
                f.write(f"REMARK Name = {name}\n")
                f.write(f"REMARK SMILES = {smiles}\n")
                f.write(f"REMARK Vina score: {score}\n")
                f.write(f"REMARK Distance from crystal ligand centroid: {distance} A\n")
                for coord in coords:
                    f.write(coord + "\n")
                f.write("ENDMDL\n")

def main():
    pdb_file_path = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/6IRT_lig_HTVS_sort_filter.pdb'
    output_dir = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/clusters/'

    # Extract ligand names, SMILES codes, dock scores, distances, and coordinates from the PDB file
    data = extract_smiles_and_coords_from_pdb(pdb_file_path)

    # Extract SMILES codes from the data
    smiles_list = [item[1] for item in data]

    # Generate pharmacophore features
    pharm_fps = generate_pharmacophore_features(smiles_list)

    # Convert fingerprints to binary vectors
    fingerprint_vectors = convert_fingerprints_to_vectors(pharm_fps)

    # Cluster ligands based on pharmacophore feature vectors
    cluster_labels = cluster_ligands(fingerprint_vectors)

    # Save clustered ligands to PDB files
    save_clustered_ligands(data, cluster_labels, output_dir)

if __name__ == "__main__":
    main()
