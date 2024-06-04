import re
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolsToGridImage

def extract_smiles_from_pdb(pdb_file_path):
    """
    Extracts ligand names, SMILES codes, dock scores, and distances from the PDB file.

    Args:
        pdb_file_path (str): Path to the PDB file.

    Returns:
        list: A list of tuples containing ligand names, SMILES codes, dock scores, and distances.
    """
    smiles_data = []
    with open(pdb_file_path, 'r') as file:
        lines = file.readlines()

    current_name = None
    current_smiles = None
    current_score = None
    current_distance = None

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
        elif "ENDMDL" in line and current_name and current_smiles and current_score and current_distance:
            smiles_data.append((current_name, current_smiles, current_score, current_distance))
            current_name = None
            current_smiles = None
            current_score = None
            current_distance = None

    return smiles_data

def generate_mols_grid_image(smiles_data, output_image_path):
    """
    Generates a grid image of molecules from a list of SMILES codes, ligand names, dock scores, and distances.

    Args:
        smiles_data (list): List of tuples containing ligand names, SMILES codes, dock scores, and distances.
        output_image_path (str): Path to save the output image.
    """
    mols = []
    legends = []
    for name, smiles, score, distance in smiles_data:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mols.append(mol)
            legends.append(f"{name}\n{smiles}\nScore: {score}\nDistance: {distance} A")

    img = MolsToGridImage(
        mols, 
        molsPerRow=4, 
        subImgSize=(300, 300),  # Increased size for higher resolution
        legends=legends
    )
    img.save(output_image_path)

def main():
    pdb_file_path = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/6IRT_lig_HTVS_sort_filter.pdb'
    output_image_path = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/ligands_grid_image.png'

    # Extract SMILES codes, ligand names, dock scores, and distances from the PDB file
    smiles_data = extract_smiles_from_pdb(pdb_file_path)

    # Generate and save the grid image
    generate_mols_grid_image(smiles_data, output_image_path)

if __name__ == "__main__":
    main()
