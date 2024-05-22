import os
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

# Directory containing the docked ligands
ligands_dir = 'docked_ligands/'
output_file = 'ranked_ligands_by_rmsd.txt'
reference_file = 'crystallized_ligand.pdb'
reference_ligand_name = 'LIG'  # Change this to the residue name of the ligand in the reference PDB

# Ensure the output file does not exist
if os.path.exists(output_file):
    os.remove(output_file)

# Load the reference ligand
cmd.load(reference_file, 'reference')
cmd.remove('not resn ' + reference_ligand_name)  # Keep only the ligand
cmd.save('reference_ligand.pdb', 'reference')

# Function to calculate RMSD using PyMOL
def calculate_rmsd(ref_file, target_file):
    cmd.load(ref_file, 'ref')
    cmd.load(target_file, 'target')
    cmd.remove('not resn ' + reference_ligand_name)  # Keep only the ligand
    rmsd = cmd.align('target', 'ref')[0]
    cmd.delete('all')
    return rmsd

# Dictionary to store RMSD values
rmsd_values = {}

# Loop through each PDB file in the directory
for pdb_file in os.listdir(ligands_dir):
    if pdb_file.endswith('.pdb'):
        pdb_path = os.path.join(ligands_dir, pdb_file)
        try:
            rmsd = calculate_rmsd('reference_ligand.pdb', pdb_path)
            rmsd_values[pdb_file] = rmsd
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

# Sort the ligands by RMSD
sorted_rmsd = sorted(rmsd_values.items(), key=lambda x: x[1])

# Write the ranked ligands to the output file
with open(output_file, 'w') as out_file:
    for ligand, rmsd in sorted_rmsd:
        out_file.write(f'{ligand}: {rmsd}\n')

print(f'Ranked ligands have been written to {output_file}')
