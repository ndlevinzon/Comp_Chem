# Jonathan Borowsky and Nathan Levinzon, UCSF 2023
import pickle
import math
import time
import re
import pandas as pd
from PIL import Image
from collections import deque
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import HybridizationType, BondType
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles

# Get the current time before running the code
start_time = time.time()
print(f"Starting Time: {start_time}")


# --------------------------Extremity Trimmer-------------------------- #


def trim_extremities(smiles):
    """Trim Parent Molecule Extremity Atoms Once"""
    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

    # Create a copy of the molecule to trim
    rw_mol = Chem.RWMol(mol)

    def get_atom_degree(mol, atom_idx):
        # Count the number of neighboring bonds for the atom
        degree = 0
        for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
            if bond.GetBeginAtomIdx() == atom_idx or bond.GetEndAtomIdx() == atom_idx:
                degree += 1
        return degree

    # Get the atom indices for extremities
    extremity_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if get_atom_degree(mol, atom.GetIdx()) == 1]

    def trim_extremity(trimmed_mol, extremity_atoms):
        # Remove one atom from each extremity atom
        for atom_idx in extremity_atoms:
            new_mol = Chem.RWMol(trimmed_mol)
            new_mol.RemoveAtom(atom_idx)
            analogs.append(new_mol.GetMol())

    # Trim the molecule once on each extremity
    trim_extremity(rw_mol, extremity_atoms)

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if mol != mol and Chem.MolToSmiles(mol) != smiles]
    return analogs

# --------------------------Bond Order Changes-------------------------- #

def BO_stepup(smiles):
    """Recursively Increases Bond Order Until Maximum Conjugation"""
    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

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

        n_neighbors1 = len([n for n in bond.GetBeginAtom().GetNeighbors()])
        n_neighbors2 = len([n for n in bond.GetEndAtom().GetNeighbors()])

        if bond.GetBondType() == Chem.BondType.SINGLE and n_neighbors1 < 3 and n_neighbors2 < 3:
            # Increase bond order from single to double
            try:
                mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(Chem.BondType.DOUBLE)
                analogs.append(Chem.RWMol(mol_copy_new))
            except:
                print("Error: Could not make Double Bond")
        elif bond.GetBondType() == Chem.BondType.DOUBLE and n_neighbors1 < 2 and n_neighbors2 < 2:
            # Increase bond order from double to triple
            try:
                mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(Chem.BondType.TRIPLE)
                analogs.append(Chem.RWMol(mol_copy_new))
            except:
                print("Error: Could not make Triple Bond")

        increase_bond_order_recursive(bond_indices[1:], Chem.RWMol(mol_copy))

    # Start the recursive function to increase bond order
    increase_bond_order_recursive(c_c_bond_indices, rw_mol)

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if mol != mol and Chem.MolToSmiles(mol) != smiles]
    return analogs


def BO_stepdown(smiles):
    """Decreases Bond Order of Non-Single Bonds"""

    def convert_to_explicit_bonds(smiles):
        """Converts SMILES into a representation with explicit bond notation"""
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol)  # Ensure explicit bond notation
        return Chem.MolToSmiles(mol, kekuleSmiles=True)

    def decrease_bond_order_recursive(smiles, bond_indices):
        """Recursively decreases bond order in the SMILES string"""
        analogs = []

        if len(bond_indices) == 0:
            return analogs

        bond_idx = bond_indices[0]
        bond_symbol = smiles[bond_idx]

        if bond_symbol == "#":
            # Decrease triple bond to double bond
            new_smiles = smiles[:bond_idx] + "=" + smiles[bond_idx + 1:]
            analogs.append(new_smiles)
            analogs.extend(decrease_bond_order_recursive(new_smiles, bond_indices[1:]))

        elif bond_symbol == "=":
            # Decrease double bond to single bond
            new_smiles = smiles[:bond_idx] + "-" + smiles[bond_idx + 1:]
            analogs.append(new_smiles)
            analogs.extend(decrease_bond_order_recursive(new_smiles, bond_indices[1:]))

        return analogs

    explicit_smiles = convert_to_explicit_bonds(smiles)
    bond_indices = [i for i, symbol in enumerate(explicit_smiles) if symbol in ("#", "=")]
    analogs_smiles = decrease_bond_order_recursive(explicit_smiles, bond_indices)
    # Remove duplicates and the initial SMILES from the analogs list
    analogs_smiles = [mol for mol in analogs_smiles if mol != smiles]
    analogs = [Chem.MolFromSmiles(smiles) for smiles in analogs_smiles]
    return analogs
# --------------------------Ring Changes-------------------------- #

def ring_breaker(smiles):
    """Enumerates Rings In Parent And Opens Rings"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

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
            if bond is not None:
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
                    Chem.SanitizeMol(mol_copy,
                                     sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                    Chem.Kekulize(mol_copy)

                    # Process sanitized analogs
                    sanitized_analogs = [mol_copy]
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
                except Chem.AtomValenceException:
                    # Skip if valence adjustment fails during sanitization
                    continue

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]
    return analogs


def ring_maker(smiles):
    """Enumerates Terminal -CH3, Finds Paths Of Length (4, 5) And Forms Rings With sp3 Carbons On Path"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

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

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]
    return analogs

# --------------------------Atom Walks-------------------------- #

def walks(smiles, target_num):
    """Performs Walks On Parent Molecule"""
    mol = Chem.MolFromSmiles(smiles)
    analogs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
            analog = Chem.RWMol(mol)
            a = analog.GetAtomWithIdx(atom.GetIdx())
            a_neighbors = [n for n in a.GetNeighbors()]
            if len(a_neighbors) == 2:
                a.SetAtomicNum(target_num)
                # Chem.SanitizeMol(analog)
                analogs.append(analog)

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]
    return analogs


def heterocycle_walks(smiles):
    """Performs Heterocycle Walks On Parent Molecule"""
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    # Enumerate heterocycles and append RWMol objects to analogs list
    for mol in EnumerateHeterocycles(mol, depth=1):
        analog = Chem.RWMol(mol)
        analogs.append(analog)

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]
    return analogs


# Functions to perform specific walking methods
def nitrogen_walks(smiles): return walks(smiles, target_num=7)
def oxygen_walks(smiles): return walks(smiles, target_num=8)
def sulfur_walks(smiles): return walks(smiles, target_num=16)

# --------------------------Atom Scans-------------------------- #

def scanning(smiles, atomic_num, n_hs):
    """Scans Parent Molecule Aromatic + Aliphatic Carbons"""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
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

            if atom.GetIsAromatic() and atomic_num == 6:
                aromatic_atoms_copy = aromatic_atoms.copy()
                aromatic_atoms_copy.remove(atom.GetIdx())  # Remove the current atom index from the list
                for position in aromatic_atoms_copy:
                    a = analog.GetAtomWithIdx(position)
                    a.SetAtomicNum(atomic_num)
                    a.SetNumExplicitHs(n_hs)
                    Chem.SanitizeMol(analog)
                    analogs.append(analog)

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]
    return analogs

# Functions to perform specific scanning methods
def methyl_scanning(smiles): return scanning(smiles, 6, 3)
def amine_scanning(smiles): return scanning(smiles, 7, 2)
def hydroxyl_scanning(smiles): return scanning(smiles, 8, 1)
def thiol_scanning(smiles): return scanning(smiles, 16, 1)
def fluorine_scanning(smiles): return scanning(smiles, 9, 0)
def chlorine_scanning(smiles): return scanning(smiles, 17, 0)
def bromine_scanning(smiles): return scanning(smiles, 35, 0)

# --------------------------Atom Neutralization-------------------------- #

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

# --------------------------Switchboard--------------------------#

def main():
    # Define the analogue methods and their corresponding names
    analog_methods = [
        # [trim_extremities, "trim"],
        # [BO_stepup, "increase-bond-order"],
        [BO_stepdown, "decrease-bond-order"],
        # [ring_breaker, "ring-opening"],
        # [ring_maker, "ring-closure"],
        [nitrogen_walks, "n-walk"],
        [oxygen_walks, "o-walk"],
        [sulfur_walks, "s-walk"],
        # [heterocycle_walks, "heterocycle-walk"],
        # [methyl_scanning, "ch3-scan"],
        # [amine_scanning, "nh2-scan"],
        # [hydroxyl_scanning, "oh-scan"],
        # [thiol_scanning, "sh-scan"],
        # [fluorine_scanning, "f-scan"],
        # [chlorine_scanning, "cl-scan"],
        # [bromine_scanning, "br-scan"],
    ]

    # Specify the input and output file names
    path = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/'
    smiles_input_filename = 'rings.smi'
    output_file_prefix = "phenol_walks_test"

    # Define the list of parent SMILES to take pictures of
    take_picture = []  # Add your desired parent SMILES here

    # Read the input file and store the smiles and zinc IDs in a DataFrame
    smiles_zinc_input = pd.read_csv(f'{path}{smiles_input_filename}', sep=' ', header=None, names=['Smiles', 'ZincID'])

    analog_key = []
    lines_out = []
    analog_count = 0

    # Loop through each input molecule
    for _, row in smiles_zinc_input.iterrows():
        smiles = row['Smiles']
        zinc_id = row['ZincID']

        lines_out.append(f"{smiles} {zinc_id}")
        analog_key_line = [[smiles, zinc_id], []]
        all_analogs = []

        # Loop through each analogue method
        for method in analog_methods:
            analog_key_scan = [method[1], []]

            # Generate analogues using the specified method
            for analog in method[0](smiles):
                analog_smiles = Chem.MolToSmiles(analog)

                if analog_smiles not in all_analogs:
                    fakezinc = "fakezinc" + str(analog_count).zfill(8)
                    lines_out.append(f"{analog_smiles} {fakezinc}")
                    analog_key_scan[1].append([analog_smiles, fakezinc])
                    all_analogs.append(analog_smiles)
                    analog_count += 1
                else:
                    print(analog_smiles)

            analog_key_line[1].append(analog_key_scan)

        analog_key.append(analog_key_line)

    # Write the generated smiles and fake zincs to a file
    output_file = f"{path}{output_file_prefix}-analogs-i{len(smiles_zinc_input)}-o{len(lines_out)}.smi"
    with open(output_file, "w") as f2:
        for line in lines_out:
            f2.write(line + "\n")

    # Write the key specifying which molecules are analogues of which other molecules and by which analoging processes
    key_file = f"{path}{output_file_prefix}-key-i{len(smiles_zinc_input)}-o{len(lines_out)}.out"
    with open(key_file, "wb") as fp:
        pickle.dump(analog_key, fp)

    # Calculate the number of analogs generated
    print(f"Molecules: {len(lines_out)} molecules")

    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Runtime: {elapsed_time} seconds")

    # Calculate molecules generated per second
    benchmark = float(len(lines_out)) / elapsed_time
    print(f"Benchmark: {benchmark} molecules/second")

    # Check if the "take_picture" list is not empty
    if take_picture:
        # Load the key file specifying which molecules are analogues
        with open(key_file, "rb") as fp:
            analog_key = pickle.load(fp)
            for analog_line in analog_key:
                parent_smiles = analog_line[0][0]
                if parent_smiles in take_picture:
                    for analog_method in analog_line[1]:

                        method = analog_method[0]
                        analog_data = analog_method[1]

                        # Create a list to store analog images
                        analog_images = []

                        # Loop through the analogues and add them to the list
                        for analog in analog_data:
                            analog_smiles = analog[0]
                            fakezinc = analog[1]
                            mol = Chem.MolFromSmiles(analog_smiles)
                            if mol is not None:  # Check if the molecule object is valid
                                img = Draw.MolToImage(mol, size=(50, 50))
                                analog_images.append(img)
                            else:
                                print(f"Invalid molecule: {analog_smiles}")

                        # Determine the dimensions of the grid
                        num_analogs = len(analog_images)
                        num_rows = min(10, int(math.ceil(num_analogs / 10)))
                        num_cols = min(10, num_analogs)

                        # Create the grid image
                        grid_image = Image.new('RGB', (num_cols * 500, num_rows * 500))

                        # Loop through the analogs and paste them into the grid
                        for i, analog_image in enumerate(analog_images):
                            col_idx = i % num_cols
                            row_idx = i // num_cols
                            offset_x = col_idx * 500
                            offset_y = row_idx * 500
                            resized_analog_image = analog_image.resize((400, 400))  # Resize the analog image if needed
                            grid_image.paste(resized_analog_image, (offset_x, offset_y))

                        # Save the grid image as a PNG file
                        try:
                            grid_image.save(f"{path}{parent_smiles}-{method}.png")
                        except Exception as e:
                            print(f"Error saving the grid image: {e}")


if __name__ == '__main__':
    main()
