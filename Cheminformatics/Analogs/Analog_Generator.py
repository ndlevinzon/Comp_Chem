import pickle
import math
import time
import re
import pandas as pd
from PIL import Image, ImageDraw
from collections import deque
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdchem import HybridizationType, BondType


# Get the current time before running the code
start_time = time.time()
print(f"Starting Time: {start_time}")


def trim_extremities(smiles):
    """Trim Parent Molecule Extremity Atoms Once"""
    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Create a copy of the molecule to trim
    trimmed_mol = Chem.RWMol(mol)

    def get_atom_degree(mol, atom_idx):
        # Count the number of neighboring bonds for the atom
        degree = 0
        for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
            if bond.GetBeginAtomIdx() == atom_idx or bond.GetEndAtomIdx() == atom_idx:
                degree += 1
        return degree

    # Get the atom indices for extremities
    extremity_atoms = [atom.GetIdx() for atom in trimmed_mol.GetAtoms()
                       if get_atom_degree(trimmed_mol, atom.GetIdx()) == 1]

    def trim_extremity(trimmed_mol, extremity_atoms):
        # Remove one atom from each extremity atom
        for atom_idx in extremity_atoms:
            new_mol = Chem.RWMol(trimmed_mol)
            new_mol.RemoveAtom(atom_idx)
            analogs.append(new_mol.GetMol())

    # Trim the molecule once on each extremity
    trim_extremity(trimmed_mol, extremity_atoms)

    # Convert the RWMol back to a Mol object for storage
    analogs.append(trimmed_mol.GetMol())

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def BO_stepup(smiles):
    """Recursively Increases Bond Order Until Maximum Conjugation"""
    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Create a copy of the molecule as an RWMol object
    rw_mol = Chem.RWMol(mol)

    # Find indices of carbon-carbon single bonds
    c_c_bond_indices = [bond.GetIdx() for bond in rw_mol.GetBonds() if
                        bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6
                        and bond.GetBondType() == Chem.BondType.SINGLE]

    def increase_bond_order_recursive(bond_indices, mol_copy):
        # Base case: no more bond indices to process
        if not bond_indices:
            analogs.append(Chem.RWMol(mol_copy))
            return

        # Process the first bond index
        bond_idx = bond_indices[0]
        bond = mol_copy.GetBondWithIdx(bond_idx)

        mol_copy_new = Chem.RWMol(mol_copy)

        if bond.GetBondType() == Chem.BondType.SINGLE:
            # Increase bond order from single to double
            mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(Chem.BondType.DOUBLE)
            increase_bond_order_recursive(bond_indices[1:], mol_copy_new)
        elif bond.GetBondType() == Chem.BondType.DOUBLE:
            # Increase bond order from double to triple
            mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(Chem.BondType.TRIPLE)
            increase_bond_order_recursive(bond_indices[1:], mol_copy_new)

        # Process the next bond index without changing the current molecule copy
        increase_bond_order_recursive(bond_indices[1:], Chem.RWMol(mol_copy))

    # Start the recursive function to increase bond order
    increase_bond_order_recursive(c_c_bond_indices, rw_mol)

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def BO_stepdown(smiles):
    """Recursively Changes Bond Order Until All Bonds Are Single"""
    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Create a copy of the molecule as an RWMol object
    rw_mol = Chem.RWMol(mol)

    # Find indices of carbon-carbon bonds
    c_c_bond_indices = [bond.GetIdx() for bond in rw_mol.GetBonds() if
                        bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6]

    def decrease_bond_order_recursive(bond_indices, mol_copy):
        # Base case: no more bond indices to process
        if not bond_indices:
            analogs.append(Chem.RWMol(mol_copy))
            return

        # Process the first bond index
        bond_idx = bond_indices[0]
        bond = mol_copy.GetBondWithIdx(bond_idx)

        for new_bond_order in range(int(bond.GetBondTypeAsDouble()) - 1, 0, -1):
            # Create a new molecule copy with decreased bond order
            mol_copy_new = Chem.RWMol(mol_copy)
            mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(Chem.BondType(new_bond_order))
            decrease_bond_order_recursive(bond_indices[1:], mol_copy_new)

        # Process the next bond index without changing the current molecule copy
        decrease_bond_order_recursive(bond_indices[1:], Chem.RWMol(mol_copy))

    # Start the recursive function to decrease bond order
    decrease_bond_order_recursive(c_c_bond_indices, rw_mol)

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def ring_breaker(smiles):
    """Enumerates Rings In Parent And Opens Rings"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Generate a 2D depiction of the molecule (optional)
    AllChem.Compute2DCoords(mol)

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

    # Disable kekulization
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Enumerate all cycles within the molecule
    cycles = Chem.GetSymmSSSR(mol)

    # Iterate over each cycle and break the bonds
    for cycle in cycles:
        for i in range(len(cycle)):
            atom1 = cycle[i]
            atom2 = cycle[(i + 1) % len(cycle)]

            # Create a copy of the molecule
            mol_copy = Chem.RWMol(mol)

            # Remove the bond between atom1 and atom2
            bond = mol_copy.GetBondBetweenAtoms(atom1, atom2)
            mol_copy.RemoveBond(atom1, atom2)

            # Add new atoms and bonds based on modified adjacency matrix
            adjacency_matrix = Chem.GetAdjacencyMatrix(mol_copy)
            for j in range(len(adjacency_matrix)):
                for k in range(j + 1, len(adjacency_matrix)):
                    if adjacency_matrix[j][k] and not mol_copy.GetBondBetweenAtoms(j, k):
                        mol_copy.AddBond(j, k, Chem.BondType.SINGLE)

            # Update the molecule properties and sanitize
            mol_copy.UpdatePropertyCache(strict=False)
            try:
                Chem.SanitizeMol(mol_copy, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
                Chem.Kekulize(mol_copy)

                # Process sanitized analogs
                sanitized_analogs = []
                sanitized_analogs.append(mol_copy)
                for analog in sanitized_analogs:
                    Chem.GetSymmSSSR(analog)
                    smiles = Chem.MolToSmiles(analog, isomericSmiles=True)
                    smiles = re.sub(r'\(\)', '', smiles)  # Remove empty parentheses
                    smiles = re.sub(r'\[.*?\]', '', smiles)  # Remove brackets denoting bond breakages
                    if smiles.startswith('*'):  # Check if SMILES starts with *
                        smiles = smiles[1:]  # Remove the asterisk character from the beginning
                    sanitized_analog = Chem.MolFromSmiles(smiles)
                    if sanitized_analog is not None:
                        analogs.append(sanitized_analog)
            except:
                pass

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def ring_maker(smiles):
    """Enumerates Terminal -CH3, Finds Paths Of Length (4, 5) And Forms Rings With sp3 Carbons On Path"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Find terminal -CH3 atoms
    terminal_CH3_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if
                          atom.GetSymbol() == 'C' and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3]

    # Define path lengths to search for, where the number is (ring_size - 1)
    path_lengths = [4, 5]

    # Store the modified molecules at each step
    modified_mols = [Chem.RWMol(mol)]

    # Iterate over the number of ring closures
    for num_closures in range(1, len(terminal_CH3_atoms) + 1):
        # Get the molecule at the previous step
        prev_mol = modified_mols[num_closures - 1]

        # Create a writable molecule for the current step
        curr_mol = Chem.RWMol(prev_mol)

        # Iterate over each terminal -CH3 atom
        for atom_idx in terminal_CH3_atoms:
            atom = curr_mol.GetAtomWithIdx(atom_idx)
            visited = set()
            queue = deque([(atom_idx, 0, [])])  # (atom_idx, current_path_length, current_path)

            # Perform breadth-first search
            while queue:
                curr_atom_idx, curr_path_length, curr_path = queue.popleft()
                visited.add(curr_atom_idx)

                # Check if the current path length matches desired lengths
                if curr_path_length in path_lengths and curr_atom_idx != atom_idx:
                    end_atom = curr_mol.GetAtomWithIdx(curr_atom_idx)
                    if end_atom.GetHybridization() == HybridizationType.SP3:
                        # print(f"Path length: {curr_path_length},
                        # Final atom: {end_atom.GetSymbol()},
                        # Hybridization: SP3")

                        # Create a new bond between the initial and final atom
                        curr_mol.AddBond(atom_idx, curr_atom_idx, BondType.SINGLE)

                        # Append the molecule to the list of analogs
                        analog = curr_mol.GetMol()
                        analogs.append(analog)

                        # Reset the molecule for the next iteration
                        curr_mol = Chem.RWMol(prev_mol)

                # Skip further exploration if current path length exceeds maximum desired length
                if curr_path_length >= max(path_lengths):
                    continue

                # Explore neighbors of the current atom
                neighbors = curr_mol.GetAtomWithIdx(curr_atom_idx).GetNeighbors()
                for neighbor in neighbors:
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:
                        new_path = curr_path + [neighbor_idx]
                        queue.append((neighbor_idx, curr_path_length + 1, new_path))

        # Store the modified molecule for the current step
        modified_mols.append(curr_mol.GetMol())

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def nitrogen_scanning(smiles, target_num=7):
    """Performs Nitrogen Scanning On Parent Molecule"""
    mol = Chem.MolFromSmiles(smiles)
    analogs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
            analog = Chem.RWMol(mol)
            a = analog.GetAtomWithIdx(atom.GetIdx())
            a_neighbors = [n for n in a.GetNeighbors()]
            if len(a_neighbors) == 2:
                a.SetAtomicNum(target_num)
                Chem.SanitizeMol(analog)
                analogs.append(analog)

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


def scanning(smiles, atomic_num, n_hs):
    """Scans Parent Molecule Aromatic + Aliphatic Carbons"""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    analogs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            n = atom.GetNeighbors()
            assert len(n) == 1
            n = n[0]
            analog = Chem.RWMol(mol)
            a = analog.GetAtomWithIdx(atom.GetIdx())
            a.SetAtomicNum(atomic_num)
            a.SetNumExplicitHs(n_hs)
            Chem.SanitizeMol(analog)
            analogs.append(analog)

    # Convert the list of RWMol objects to a set to remove duplicates
    analogs = list(set(analogs))
    return analogs


# Functions to perform specific scanning methods
def methyl_scanning(smiles):
    return scanning(smiles, 6, 3)


def amine_scanning(smiles):
    return scanning(smiles, 7, 2)


def hydroxyl_scanning(smiles):
    return scanning(smiles, 8, 1)


def fluorine_scanning(smiles):
    return scanning(smiles, 9, 0)


def main():
    # Define the analogue methods and their corresponding names
    analogue_methods = [
        [trim_extremities, "trim"],
        [BO_stepup, "increase-bond-order"],
        [BO_stepdown, "decrease-bond-order"],
        [ring_breaker, "ring-opening"],
        [ring_maker, "ring-closure"],
        [nitrogen_scanning, "n-scan"],
        [methyl_scanning, "ch3-scan"],
        [amine_scanning, "nh2-scan"],
        [hydroxyl_scanning, "oh-scan"],
        [fluorine_scanning, "f-scan"]
    ]

    # Specify the input and output file names
    path = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/'
    smiles_input_filename = 'smi-zn-ampc-all.smi'
    output_file_prefix = "total_benchmark"

    # Define the list of parent SMILES to take pictures of
    take_picture = [""]  # Add your desired parent SMILES here

    # Read the input file and store the smiles and zinc IDs in a DataFrame
    smiles_zinc_input = pd.read_csv(f'{path}{smiles_input_filename}', sep=' ', header=None, names=['Smiles', 'ZincID'])

    analogue_key = []
    lines_out = []
    analog_count = 0

    # Loop through each input molecule
    for _, row in smiles_zinc_input.iterrows():
        smiles = row['Smiles']
        zinc_id = row['ZincID']

        lines_out.append(f"{smiles} {zinc_id}")
        analogue_key_line = [[smiles, zinc_id], []]
        all_analogues = []

        # Loop through each analogue method
        for method in analogue_methods:
            analogue_key_scan = [method[1], []]

            # Generate analogues using the specified method
            for analogue in method[0](smiles):
                analogue_smiles = Chem.MolToSmiles(analogue, isomericSmiles=True)

                if analogue_smiles not in all_analogues:
                    fakezinc = "fakezinc" + str(analog_count).zfill(8)
                    lines_out.append(f"{analogue_smiles} {fakezinc}")
                    analogue_key_scan[1].append([analogue_smiles, fakezinc])
                    all_analogues.append(analogue_smiles)
                    analog_count += 1
                else:
                    print(analogue_smiles)

            analogue_key_line[1].append(analogue_key_scan)

        analogue_key.append(analogue_key_line)

    # Write the generated smiles and fake zincs to a file
    output_file = f"{path}{output_file_prefix}-analogues-i{len(smiles_zinc_input)}-o{len(lines_out)}.smi"
    with open(output_file, "w") as f2:
        for line in lines_out:
            f2.write(line + "\n")

    # Write the key specifying which molecules are analogues of which other molecules and by which analoging processes
    key_file = f"{path}{output_file_prefix}-key-i{len(smiles_zinc_input)}-o{len(lines_out)}"
    with open(key_file, "wb") as fp:
        pickle.dump(analogue_key, fp)

    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Runtime: {elapsed_time} seconds")

    # Calculate molecules generated per second
    benchmark = float(len(lines_out)) / elapsed_time
    print(f"Benchmark: {benchmark} molecules/second")

    # Check if the "take_picture" list is not empty
    if take_picture:
        # Generate and save the grid image for each parent SMILES in take_picture
        for analogue_line in analogue_key:
            parent_smiles = analogue_line[0][0]
            if parent_smiles in take_picture:
                for analogue_method in analogue_line[1]:
                    method = analogue_method[0]
                    analogue_data = analogue_method[1]

                    # Create a list to store analogue images
                    analogue_images = []

                    # Loop through the analogues and add them to the list
                    for analogue in analogue_data:
                        analogue_smiles = analogue[0]
                        fakezinc = analogue[1]
                        mol = Chem.MolFromSmiles(analogue_smiles)
                        img = Draw.MolToImage(mol, size=(500, 500))
                        analogue_images.append(img)

                    # Determine the dimensions of the grid
                    num_analogues = len(analogue_images)
                    num_rows = min(10, int(math.ceil(num_analogues / 10)))
                    num_cols = min(10, num_analogues)

                    # Create the grid image
                    grid_image = Image.new('RGB', (num_cols * 500, num_rows * 500))

                    # Loop through the analogues and paste them into the grid
                    for i, analogue_image in enumerate(analogue_images):
                        col_idx = i % 10
                        row_idx = i // 10
                        grid_image.paste(analogue_image, (col_idx * 500, row_idx * 500))

                    # Add the parent SMILES and method label to the top of the image
                    draw = ImageDraw.Draw(grid_image)
                    label = f"Parent SMILES: {parent_smiles}\nMethod: {method}"
                    draw.text((10, 10), label, fill="white")

                    # Save the grid image as a PNG file
                    grid_image.save(f"{path}{output_file_prefix}-{parent_smiles}-{method}.png")


if __name__ == '__main__':
    main()
