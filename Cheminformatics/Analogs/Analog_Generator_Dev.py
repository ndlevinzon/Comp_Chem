# Nathan Levinzon, Jonathan Borowsky, and Olivier Mailhot | UCSF 2023
# Non-refactored Version
import time
import re
import pandas as pd
from collections import defaultdict
from collections import deque
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import HybridizationType, BondType
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles

# --------------------------Extremity Trimmer-------------------------- #

def trim_extremities(smiles):
    """Trim Parent Molecule Extremity Atoms One At A Time if M.W. > 500 Da"""

    # Convert the SMILES code to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    # Calculate the molecular weight of the parent molecule
    mw = Descriptors.MolWt(mol)

    # Exit the function if the molecular weight is less than or equal to 500 Da
    if mw <= 500:
        return []

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
    if not mol:
        return []

    analogs = []

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

    # Create a copy of the molecule as an RWMol object
    rw_mol = Chem.RWMol(mol)

    # Find indices of carbon-carbon single bonds
    c_c_bond_indices = [bond.GetIdx() for bond in rw_mol.GetBonds() if
                        bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6
                        and bond.GetBondType() == BondType.SINGLE]

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

        if bond.GetBondType() == BondType.SINGLE and n_neighbors1 < 3 and n_neighbors2 < 3:
            # Increase bond order from single to double
            try:
                mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(BondType.DOUBLE)
                Chem.AdjustBondStereo(mol_copy_new, bond_idx)  # Adjust bond stereochemistry
                analogs.append(Chem.RWMol(mol_copy_new))
            except:
                print("Error: Could not make Double Bond")
        elif bond.GetBondType() == BondType.DOUBLE and n_neighbors1 < 2 and n_neighbors2 < 2:
            # Increase bond order from double to triple
            try:
                mol_copy_new.GetBondWithIdx(bond_idx).SetBondType(BondType.TRIPLE)
                Chem.AdjustBondStereo(mol_copy_new, bond_idx)  # Adjust bond stereochemistry
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
        if not mol:
            return []

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
    analogs_smiles = list(set(analogs_smiles) - {smiles})
    analogs = [Chem.MolFromSmiles(smiles) for smiles in analogs_smiles]

    return analogs


# --------------------------Ring Changes-------------------------- #

def ring_breaker(smiles):
    """Enumerates Rings In Parent And Opens Rings"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

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

            # Remove the bond between atom1 and atom2 if it exists
            bond = mol_copy.GetBondBetweenAtoms(atom1, atom2)
            if bond is not None:
                mol_copy.RemoveBond(atom1, atom2)

                # Add new atoms and bonds based on the modified adjacency matrix
                adjacency_matrix = Chem.GetAdjacencyMatrix(mol_copy)
                for j in range(len(adjacency_matrix)):
                    for k in range(j + 1, len(adjacency_matrix)):
                        if adjacency_matrix[j][k] and not mol_copy.GetBondBetweenAtoms(j, k):
                            mol_copy.AddBond(j, k, BondType.SINGLE)

                # Update the molecule properties and sanitize with kekulization
                mol_copy.UpdatePropertyCache(strict=False)
                try:
                    Chem.SanitizeMol(mol_copy, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                    Chem.SanitizeMol(mol_copy, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE)

                    # Process sanitized analogs
                    sanitized_analogs = [mol_copy]
                    for analog in sanitized_analogs:
                        Chem.GetSymmSSSR(analog)
                        smiles = Chem.MolToSmiles(analog)
                        smiles = re.sub(r'\(\)', '', smiles)  # Remove empty parentheses
                        smiles = re.sub(r'\[.*?\]', '', smiles)  # Remove brackets denoting bond breakages
                        if smiles.startswith('*'):  # Check if SMILES starts with *
                            smiles = smiles[1:]  # Remove the asterisk character from the beginning
                        sanitized_analog = Chem.MolFromSmiles(smiles)
                        if sanitized_analog is not None:
                            analogs.append(sanitized_analog)
                except (Chem.AtomValenceException, Chem.KekulizeException):
                    # Skip if valence adjustment or kekulization fails during sanitization
                    continue

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]

    return analogs


def ring_maker(smiles):
    """Enumerates Terminal -CH3, Finds Paths Of Length (4, 5) And Forms Rings With sp3 Carbons On Path"""
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

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

# --------------------------Walks-------------------------- #

def walks(smiles, target_num):
    """Performs Walks On Parent Molecule"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    analogs = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
            analog = Chem.RWMol(mol)
            a = analog.GetAtomWithIdx(atom.GetIdx())
            a_neighbors = [n for n in a.GetNeighbors()]
            if len(a_neighbors) == 2:
                a.SetAtomicNum(target_num)
                try:
                    analog_mol = Chem.MolFromSmiles(Chem.MolToSmiles(analog, isomericSmiles=True))
                    if analog_mol is not None:
                        analogs.append(analog_mol)
                except Chem.KekulizeException:
                    # Skip the molecule if kekulization fails
                    print(f"Kekulization error for molecule: {Chem.MolToSmiles(analog)}")
                    continue

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if mol is not None and Chem.MolToSmiles(mol) != smiles]

    return analogs

def heterocycle_walks(smiles):
    """Performs Heterocycle Walks On Parent Molecule"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    analogs = []

    # Enumerate heterocycles and append RWMol objects to analogs list
    for mol in EnumerateHeterocycles(mol, depth=1):
        try:
            analog = Chem.RWMol(mol)
            analogs.append(analog)
        except Chem.KekulizeException:
            # Skip the molecule if kekulization fails
            print(f"Kekulization error for molecule: {Chem.MolToSmiles(mol)}")
            continue

    # Remove duplicates and the initial SMILES from the analogs list
    analogs = [mol for mol in analogs if Chem.MolToSmiles(mol) != smiles]

    return analogs

# Functions to perform specific walking methods

def nitrogen_walks(smiles): return walks(smiles, target_num=7)
def oxygen_walks(smiles): return walks(smiles, target_num=8)
def sulfur_walks(smiles): return walks(smiles, target_num=16)

# --------------------------Scans-------------------------- #

def scanning(smiles, fragments):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)

    analogs = set()  # Use a set to store unique analogs

    for hydrogen in mol.GetAtoms():
        if hydrogen.GetSymbol() == 'H':
            parent = Chem.RWMol(mol)
            hydrogen_idx = hydrogen.GetIdx()

            # Find the carbon atom connected to the hydrogen
            carbon = None
            for neighbor in hydrogen.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    carbon = neighbor
                    break

            if carbon is None:
                continue

            for fragment in fragments:
                analog = Chem.RWMol(parent)

                # Replace the hydrogen with the fragment
                analog.RemoveAtom(hydrogen_idx)
                fragment_mol = Chem.MolFromSmiles(fragment)

                if fragment_mol is not None:
                    atom_map = {}  # Map atom indices from the fragment to the analog
                    for atom in fragment_mol.GetAtoms():
                        new_atom_idx = analog.AddAtom(atom)
                        atom_map[atom.GetIdx()] = new_atom_idx
                    for bond in fragment_mol.GetBonds():
                        analog.AddBond(atom_map[bond.GetBeginAtomIdx()], atom_map[bond.GetEndAtomIdx()], bond.GetBondType())

                    # Add the bond between the carbon and the fragment
                    analog.AddBond(carbon.GetIdx(), atom_map[0], Chem.BondType.SINGLE)

                try:
                    # Sanitize the generated analog molecule
                    Chem.SanitizeMol(analog, catchErrors=True)
                    analog_smiles = Chem.MolToSmiles(analog)
                    if '.' not in analog_smiles:
                        analogs.add(analog_smiles)
                except (Chem.MolSanitizeException, Chem.KekulizeException) as e:
                    print("Molecule sanitization failed:", str(e))
                    continue

    analogs = [Chem.MolFromSmiles(analog_smiles) for analog_smiles in analogs]

    return analogs

# Functions to perform specific scanning methods

def methyl_scanning(smiles): return scanning(smiles, ['C'])
def amine_scanning(smiles): return scanning(smiles, ['N'])
def hydroxyl_scanning(smiles): return scanning(smiles, ['O'])
def thiol_scanning(smiles): return scanning(smiles, ['S'])
def fluorine_scanning(smiles): return scanning(smiles, ['F'])
def chlorine_scanning(smiles): return scanning(smiles, ['Cl'])
def bromine_scanning(smiles): return scanning(smiles, ['Br'])
def iodine_scanning(smiles): return scanning(smiles, ['I'])

def bioisosters_scanning(smiles):
    fragments = [
        'C(=O)O',  # Carboxylic acid
        'C(=O)N',  # Amide
        'C(=O)N(=O)=O',  # Nitro group
        'C#N',  # Cyanide
        'C(=O)Cl',  # Acid chloride
        'C(=O)F',  # Fluoride
        '[C](F)(F)F',  # Trifluoro
        'C(=O)Br',  # Bromide
        'C(=O)I',  # Iodide
        'C(=S)N',  # Thioamide
        'C(=O)S',  # Thioester
        'C(=N)N',  # Azide
        'C(=O)N(=O)N',  # Nitramide
        'C(=O)OCC',  # Ester
        'C(=O)OC',  # Ether
        'C(=O)NCC',  # Carbamate
        'C(=O)OCC(=O)O',  # Anhydride
        'C(c1ccccc1)',   # -CH2-Benzene
        'C(n1cccc1)' # -CH2-Pyrrole
    ]
    return scanning(smiles, fragments)

# --------------------------Stereoisomers--------------------------#

def enumerate_stereoisomers(mol):
    """Enumerate stereoisomers of the input RWMol object.

    Parameters:
        mol (RWMol): The input RWMol object.

    Returns:
        List[RWMol] or None: A list of RWMol objects representing each unique stereoisomer,
                             or None if there are no stereoisomers to enumerate.
    """
    if mol is None:
        return []

    # Find all stereocenters in the molecule
    stereo_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('_ChiralityPossible')]

    # If there are no stereocenters, return None
    if not stereo_centers:
        return []

    # Enumerate all possible combinations of stereo configurations
    num_stereo_centers = len(stereo_centers)
    num_combinations = 2 ** num_stereo_centers

    stereoisomers = []
    for i in range(num_combinations):
        isomer = Chem.RWMol(mol)
        for j in range(num_stereo_centers):
            atom_idx = stereo_centers[j]
            # Get the original atom chirality
            original_chirality = mol.GetAtomWithIdx(atom_idx).GetChiralTag()
            # Invert the chirality for this isomer based on the combination
            isomer.GetAtomWithIdx(atom_idx).InvertChirality()
            # Check if this isomer already exists in the list
            if isomer not in stereoisomers:
                stereoisomers.append(isomer)

    return stereoisomers

# --------------------------Switchboard--------------------------#

def main():

    # Get the current time before running the code
    start_time = time.time()
    print(f"Starting Time: {start_time}")

    # Define the analogue methods and their corresponding names
    analog_methods = [
        [trim_extremities, "trim"],
        [BO_stepup, "increase-bond-order"],
        [BO_stepdown, "decrease-bond-order"],
        [ring_breaker, "ring-opening"],
        [ring_maker, "ring-closure"],
        [nitrogen_walks, "n-walk"],
        [oxygen_walks, "o-walk"],
        [sulfur_walks, "s-walk"],
        [heterocycle_walks, "heterocycle-walk"],
        [methyl_scanning, "ch3-scan"],
        [amine_scanning, "nh2-scan"],
        [hydroxyl_scanning, "oh-scan"],
        [thiol_scanning, "sh-scan"],
        [fluorine_scanning, "f-scan"],
        [chlorine_scanning, "cl-scan"],
        [bromine_scanning, "br-scan"],
        [iodine_scanning, "i-scan"],
        [bioisosters_scanning, "bioisoster-scan"],
    ]

    # Specify the input and output file names
    path = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/'
    smiles_input_filename = 'alpha_2a_ligands.smi'
    output_file_prefix = "jk_test"

    # Read the input file and store the smiles and zinc IDs in a DataFrame
    smiles_zinc_input = pd.read_csv(f'{path}{smiles_input_filename}', sep=' ', header=None, names=['Smiles', 'ZincID'])

    # Create a dictionary to store the analogs for each input SMILES code and their zinc IDs
    all_analogs_dict = defaultdict(set)  # Use sets for analogs

    # Loop through each input molecule
    for _, row in smiles_zinc_input.iterrows():
        smiles = row['Smiles']
        zinc_id = row['ZincID']

        # Initialize the counter for each parent compound (reset to zero for each parent)
        analogs_count = 0

        # Loop through each analogue method for the current input molecule
        for method in analog_methods:
            # Generate analogues using the specified method
                for analog in method[0](smiles):
                    try:
                        # Enumerate stereoisomers for the current analog
                        stereoisomers = enumerate_stereoisomers(analog)
                        if stereoisomers is None:
                            analog_smiles = Chem.MolToSmiles(analog, isomericSmiles=True)
                            if not analog_smiles:
                                raise ValueError(f"Invalid SMILES: {analog_smiles}")
                            # Add analog to the set and store zinc ID
                            all_analogs_dict[smiles].add((analog_smiles, zinc_id))
                            analogs_count += 1
                            # Print analog SMILES and create fakezinc inside the loop
                            fakezinc = str(zinc_id) + "_analog" + str(analogs_count).zfill(4)
                            print(analog_smiles, fakezinc)
                        else:
                            for isomer in stereoisomers:
                                analog_smiles = Chem.MolToSmiles(isomer, isomericSmiles=True)
                                if not analog_smiles:
                                    raise ValueError(f"Invalid SMILES: {analog_smiles}")
                                # Add analog to the set and store zinc ID
                                all_analogs_dict[smiles].add((analog_smiles, zinc_id))
                                analogs_count += 1
                                # Print analog SMILES and create fakezinc inside the loop
                                fakezinc = str(zinc_id) + "_analog" + str(analogs_count).zfill(4)
                                print(analog_smiles, fakezinc)
                    except ValueError as e:
                        print(str(e))  # Skip the entry and print the error message
                        continue  # Skip to the next iteration

    # Write the generated smiles and fake zincs to a file
    total_analogs_count = 0
    output_file = f"{path}{output_file_prefix}-analogs-i{len(smiles_zinc_input)}-o{(sum(len(analogs) for _, analogs in all_analogs_dict.items()))}.smi"
    print(f"Writing to {output_file}")
    with open(output_file, "w") as f2:
        for _, row in smiles_zinc_input.iterrows():
            smiles = row['Smiles']
            zinc_id = row['ZincID']
            f2.write(f"{smiles} {zinc_id}\n")
            parent_analogs = all_analogs_dict[smiles]
            for i, (analog_smiles, parent_zinc_id) in enumerate(parent_analogs, start=1):
                fakezinc = str(parent_zinc_id) + "_analog" + str(i).zfill(4)
                f2.write(f"{analog_smiles} {fakezinc}\n")
                total_analogs_count += 1

    # Calculate the number of analogs generated
    print(f"Molecules: {total_analogs_count} molecules")

    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Runtime: {elapsed_time} seconds")

    # Calculate molecules generated per second
    benchmark = total_analogs_count / elapsed_time
    print(f"Benchmark: {benchmark} molecules/second")


if __name__ == '__main__':
    main()
