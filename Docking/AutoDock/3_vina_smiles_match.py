import re

def read_ligands_smiles(file_path):
    ligands_smiles = {}
    with open(file_path, 'r') as file:
        for line in file:
            smiles, ligand_name = line.strip().split()
            ligands_smiles[ligand_name] = smiles
    return ligands_smiles

def process_pdb_file(input_pdb_path, output_pdb_path, ligands_smiles):
    with open(input_pdb_path, 'r') as file:
        pdb_lines = file.readlines()

    updated_pdb_lines = []
    current_ligand_name = None
    for line in pdb_lines:
        updated_pdb_lines.append(line)
        if line.startswith("REMARK Name ="):
            match = re.search(r'REMARK Name = (\S+)', line)
            if match:
                current_ligand_name = match.group(1)
                if current_ligand_name in ligands_smiles:
                    smiles_remark = f"REMARK SMILES = {ligands_smiles[current_ligand_name]}\n"
                    updated_pdb_lines.append(smiles_remark)
    
    with open(output_pdb_path, 'w') as file:
        file.writelines(updated_pdb_lines)

def main():
    ligands_smiles_file = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test/tricyclics_analogs-analogs-i73-o16617.smi'
    input_pdb_file = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test/rmsd_sort.pdb'
    output_pdb_file = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test/smiles_pdb_match.pdb'

    ligands_smiles = read_ligands_smiles(ligands_smiles_file)
    process_pdb_file(input_pdb_file, output_pdb_file, ligands_smiles)

if __name__ == "__main__":
    main()
