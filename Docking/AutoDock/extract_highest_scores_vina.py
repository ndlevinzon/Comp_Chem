import os

# Directory containing the PDB files
pdb_dir = 'docked_ligands/'
output_file = 'highest_scored_ligands.txt'

# Ensure the output file does not exist
if os.path.exists(output_file):
    os.remove(output_file)

highest_scores = {}

# Function to parse PDBQT file and extract scores
def parse_pdb_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                score = float(parts[3])
                return score
    return None

# Loop through each PDB file in the directory
for pdb_file in os.listdir(pdb_dir):
    if pdb_file.endswith('.pdb'):
        pdb_path = os.path.join(pdb_dir, pdb_file)
        score = parse_pdb_file(pdb_path)
        if score is not None:
            highest_scores[pdb_file] = score

# Find the ligand with the highest score (lowest binding affinity)
highest_scored_ligands = sorted(highest_scores.items(), key=lambda x: x[1])

# Write the highest scored ligands to the output file
with open(output_file, 'w') as out_file:
    for ligand, score in highest_scored_ligands:
        out_file.write(f'{ligand}: {score}\n')

print(f'Highest scored ligands have been written to {output_file}')
