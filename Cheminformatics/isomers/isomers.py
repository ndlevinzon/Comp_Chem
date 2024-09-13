import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import subprocess

def enumerate_stereoisomers(mol):
    opts = StereoEnumerationOptions(onlyUnassigned=False)  # Consider all stereocenters
    isomers = EnumerateStereoisomers(mol, options=opts)
    return list(isomers)

def convert_mol_to_mol2(input_mol, output_mol2):
    # Convert .mol file to .mol2 using Open Babel
    subprocess.run(['obabel', input_mol, '-O', output_mol2], check=True)

def get_stereochemistry_info(mol):
    """Extracts R/S and E/Z stereochemistry info from the molecule."""
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    stereo_info = []

    # Get chirality (R/S)
    for atom_idx, chirality in chiral_centers:
        stereo_info.append(f'{chirality}{atom_idx + 1}')  # R1, S2, etc.

    # Get E/Z for double bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE:
            bond_stereo = bond.GetStereo()
            atom1_idx = bond.GetBeginAtomIdx() + 1
            atom2_idx = bond.GetEndAtomIdx() + 1
            if bond_stereo == Chem.rdchem.BondStereo.STEREOZ:
                stereo_info.append(f'Z{atom1_idx}-{atom2_idx}')
            elif bond_stereo == Chem.rdchem.BondStereo.STEREOE:
                stereo_info.append(f'E{atom1_idx}-{atom2_idx}')

    return stereo_info

def process_mol2_directory(input_dir):
    # Loop through each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith('.mol2'):
            # Construct the full file path
            file_path = os.path.join(input_dir, filename)

            # Read the mol2 file
            mol = Chem.MolFromMol2File(file_path)

            if mol is None:
                print(f"Warning: Could not read {filename}, skipping.")
                continue

            # Generate stereoisomers
            stereoisomers = enumerate_stereoisomers(mol)

            # Save each isomer back as MOL file and then convert to MOL2 format
            base_name = os.path.splitext(filename)[0]  # Remove '.mol2' from filename
            for i, iso in enumerate(stereoisomers):
                # Get stereochemistry information (R/S, E/Z)
                stereo_info = get_stereochemistry_info(iso)
                stereo_str = '_'.join(stereo_info) if stereo_info else f'isomer_{i + 1}'

                mol_filename = os.path.join(input_dir, f'{base_name}_{stereo_str}.mol')
                mol2_filename = os.path.join(input_dir, f'{base_name}_{stereo_str}.mol2')

                # Write MOL file using RDKit
                Chem.MolToMolFile(iso, mol_filename)

                # Convert the MOL file to MOL2 format using Open Babel
                convert_mol_to_mol2(mol_filename, mol2_filename)

                # Optionally, remove the intermediate MOL file if not needed
                os.remove(mol_filename)

            print(f"Processed {len(stereoisomers)} stereoisomers for {filename}")

# Example usage:
input_directory = ''  # Replace with your actual directory
process_mol2_directory(input_directory)
