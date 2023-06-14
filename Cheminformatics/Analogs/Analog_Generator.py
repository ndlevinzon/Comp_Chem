import pickle
import time
import re
import pandas as pd
from rdkit import Chem
from collections import deque
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import HybridizationType, BondType


# Get the current time before running the code
start_time = time.time()


# Function to remove duplicates from a list of molecules
def remove_duplicates(analogs):
    unique_molecules = []
    for molecule in analogs:
        smiles = Chem.MolToSmiles(molecule)
        if smiles not in unique_molecules:
            unique_molecules.append(smiles)
    return unique_molecules


def ring_breaker(smiles):
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)

    # Generate a 2D depiction of the molecule (optional)
    AllChem.Compute2DCoords(mol)

    # Adjust valence of the molecule
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

    # Disable kekulization
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Enumerate all cycles within the molecule
    cycles = Chem.GetSymmSSSR(mol)

    # Break each bond within each cycle
    analogs = []
    for cycle in cycles:
        for i in range(len(cycle)):
            atom1 = cycle[i]
            atom2 = cycle[(i + 1) % len(cycle)]

            # Create a copy of the molecule
            mol_copy = Chem.RWMol(mol)

            # Remove the bond
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

    # Remove duplicates from the list of analogs
    remove_duplicates(analogs)

    return analogs


def ring_maker(smiles):
    # Create the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)

    # Find terminal -CH3 atoms
    terminal_CH3_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if
                          atom.GetSymbol() == 'C' and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3]

    # Define path lengths to search for, where the number is (ring_size - 1)
    path_lengths = [4, 5]

    # Store the modified molecules at each step
    modified_mols = [Chem.RWMol(mol)]

    # Create an empty list to store the output molecules of newly closed structures
    analogs = []

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

    # Remove duplicates from the list of analogs
    remove_duplicates(analogs)

    # Return the list of analogs
    return analogs


# Function to perform nitrogen scanning on a molecule
def nitrogen_scanning(smiles, target_num=7):
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
    return analogs


# Function to perform aromatic scanning with different atomic numbers and hydrogen counts
def aromatic_scanning(smiles, atomic_num, n_hs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    analogs = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            n = atom.GetNeighbors()
            assert len(n) == 1
            n = n[0]
            if n.GetSymbol() == 'C' and n.GetIsAromatic():
                analog = Chem.RWMol(mol)
                a = analog.GetAtomWithIdx(atom.GetIdx())
                a.SetAtomicNum(atomic_num)
                a.SetNumExplicitHs(n_hs)
                Chem.SanitizeMol(analog)
                analogs.append(analog)
    return analogs


# Functions to perform specific scanning methods
def methyl_scanning(smiles):
    return aromatic_scanning(smiles, 6, 3)


def amine_scanning(smiles):
    return aromatic_scanning(smiles, 7, 2)


def hydroxyl_scanning(smiles):
    return aromatic_scanning(smiles, 8, 1)


def fluorine_scanning(smiles):
    return aromatic_scanning(smiles, 9, 0)


def main():

    # Define the analogue methods and their corresponding names
    analogue_methods = [
        [ring_breaker, "ring_opening"],
        [ring_maker, "ring_closure"],
        [nitrogen_scanning, "n-scan"],
        [methyl_scanning, "ch3-scan"],
        [amine_scanning, "nh2-scan"],
        [hydroxyl_scanning, "oh-scan"],
        [fluorine_scanning, "f-scan"]
    ]

    # Specify the input and output file names
    path = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/'
    smiles_input_filename = 'smi-zn-ampc-all.smi'
    output_file_prefix = "with_rings_o_c"

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

    # Get the current time after running the code
    end_time = time.time()

    # Calculate the elapsed time
    elapsed_time = end_time - start_time

    # Print the runtime
    print(f"Runtime: {elapsed_time} seconds")


if __name__ == '__main__':
    main()
