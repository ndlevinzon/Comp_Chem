import os

# Function to parse the PDB file and extract relevant information
def parse_pdb(input_pdb_file):
    residues = {}

    with open(input_pdb_file, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                # Extract necessary information from each line
                residue_number = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                atom_number = line[6:11].strip()

                # Store the information in a nested dictionary
                if residue_number not in residues:
                    residues[residue_number] = {
                        'residue_name': residue_name,
                        'atoms': {},
                    }
                residues[residue_number]['atoms'][atom_name] = atom_number

    return residues

# Function to write nucleic acid information to separate output files
def write_nucleic_acid_files(residues, atom_orders, pdb_file_name):
    nucleic_acids = {'A': [], 'U': [], 'C': [], 'G': []}

    for residue_number, data in residues.items():
        residue_name = data['residue_name']
        ordered_atom_names = atom_orders.get(residue_name, '').split()
        atom_numbers = []

        for atom_name in ordered_atom_names:
            if atom_name in data['atoms']:
                atom_number = data['atoms'][atom_name]
                atom_numbers.append(atom_number)

        nucleic_acids[residue_name].append(' '.join(atom_numbers))

    # Write the results to separate output files for each residue name
    for residue_name, atom_numbers_list in nucleic_acids.items():
        output_file = f'{pdb_file_name}.list.{residue_name}.txt'
        with open(output_file, 'w') as output_file:
            first_line = '3' if residue_name in ['A', 'U'] else '2'
            second_line = atom_orders.get(residue_name, '')
            output_file.write(f'{first_line}\n{second_line}\n')
            output_file.write('\n'.join(atom_numbers_list)

# Dictionary specifying the atom orders for each nucleic acid
atom_orders = {
    'A': "N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4",
    'U': "N1 C6 H6 C5 H5 C4 O4 N3 H3 C2 O2",
    'C': "N1 C6 H6 C5 H5 C4 N4 H41 H42 N3 C2 O2",
    'G': "N9 C8 H8 N7 C5 C6 O6 N1 H1 C2 N2 H21 H22 N3 C4",
}

# Specify the input PDB file and obtain the base file name
input_pdb_file = 'UACCGU_0_amber.pdb'  # Replace with your PDB file
pdb_file_name = os.path.splitext(os.path.basename(input_pdb_file))[0]

# Parse the PDB file and write nucleic acid information to output files
residues = parse_pdb(input_pdb_file)
write_nucleic_acid_files(residues, atom_orders, pdb_file_name)
