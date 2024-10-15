from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolfiles import MolToPDBFile
import numpy as np
import random


def load_ligand_from_pdb(pdb_file):
    print(f"Loading ligand1 from PDB: {pdb_file}")
    return Chem.rdmolfiles.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)

def load_ligand_from_mol2(mol2_file):
    print(f"Loading ligand2 from MOL2: {mol2_file}")
    return Chem.MolFromMol2File(mol2_file, removeHs=False)

def align_ligands_by_mcs(ligand1_mol, ligand2_mol):
    print("Finding Maximum Common Substructure (MCS) between ligands...")
    mcs_result = rdFMCS.FindMCS([ligand1_mol, ligand2_mol])
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    
    ligand1_match = ligand1_mol.GetSubstructMatch(mcs_mol)
    ligand2_match = ligand2_mol.GetSubstructMatch(mcs_mol)
    
    if not ligand1_match or not ligand2_match:
        raise ValueError("Could not find MCS match in one or both ligands.")
    
    atom_map = list(zip(ligand2_match, ligand1_match))
    
    print("Aligning ligand2 to ligand1 based on MCS...")
    rmsd = AllChem.AlignMol(ligand2_mol, ligand1_mol, atomMap=atom_map)
    print(f"RMSD after alignment: {rmsd:.4f}")
    return rmsd

def optimize_ligand_geometry(mol, max_iters=50000):
    print(f"Optimizing geometry using MMFF/UFF force field (max iterations: {max_iters})...")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return mol

def minimize_ligand_energy(mol, max_iters=50000):
    print(f"Performing energy minimization using MMFF (max iterations: {max_iters})...")
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    ff.Minimize(maxIts=max_iters)
    return mol

def planarize_ligand(mol):
    """Make the ligand as planar as possible by setting dihedrals to 0 or 180 degrees for conjugated bonds."""
    conf = mol.GetConformer()

    for bond in mol.GetBonds():
        if bond.IsInRing():  # Skip bonds that are part of a ring system
            continue

        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.AROMATIC]:
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()

            atom1_neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(atom1).GetNeighbors() if a.GetIdx() != atom2]
            atom2_neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(atom2).GetNeighbors() if a.GetIdx() != atom1]

            if len(atom1_neighbors) > 0 and len(atom2_neighbors) > 0:
                dihedral_atoms = [atom1_neighbors[0], atom1, atom2, atom2_neighbors[0]]
                # Set dihedral to 0 or 180 degrees to planarize the molecule
                AllChem.SetDihedralDeg(conf, dihedral_atoms[0], dihedral_atoms[1], dihedral_atoms[2], dihedral_atoms[3], 180)
    
    return mol


def straighten_bonds(mol):
    """Straighten bonds by adjusting bond angles."""
    conf = mol.GetConformer()

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()

            pos1 = np.array(conf.GetAtomPosition(atom1))
            pos2 = np.array(conf.GetAtomPosition(atom2))
            direction = pos2 - pos1
            direction /= np.linalg.norm(direction)

            new_pos2 = pos1 + direction * np.linalg.norm(pos2 - pos1)
            conf.SetAtomPosition(atom2, new_pos2)
    
    return mol

def get_rotatable_bonds(mol):
    """Identify rotatable bonds in the molecule."""
    rotatable_bonds = []
    for bond in mol.GetBonds():
        # A rotatable bond is a single, non-ring bond that is not connected to a terminal atom
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetDegree() > 1 and atom2.GetDegree() > 1:  # Avoid terminal atoms
                rotatable_bonds.append(bond)
    return rotatable_bonds

def manipulate_bonds(mol, bond, stretch_factor=0.02, angle_step=np.pi/36):
    """Rotate and stretch bonds with a focus on avoiding excessive distortion."""
    conf = mol.GetConformer()

    atom1 = bond.GetBeginAtomIdx()
    atom2 = bond.GetEndAtomIdx()

    # Get the neighbors of the bond atoms to define a dihedral
    atom1_neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(atom1).GetNeighbors() if a.GetIdx() != atom2]
    atom2_neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(atom2).GetNeighbors() if a.GetIdx() != atom1]

    if len(atom1_neighbors) == 0 or len(atom2_neighbors) == 0:
        print(f"Cannot define dihedral for bond {bond.GetIdx()}. Skipping rotation.")
        return mol

    # Define atoms for the dihedral: [atom1_neighbor, atom1, atom2, atom2_neighbor]
    dihedral_atoms = [atom1_neighbors[0], atom1, atom2, atom2_neighbors[0]]

    # Stretch the bond slightly with a smaller factor
    pos1 = np.array(conf.GetAtomPosition(atom1))
    pos2 = np.array(conf.GetAtomPosition(atom2))
    direction = pos2 - pos1
    direction /= np.linalg.norm(direction)

    # Stretch by a smaller factor to avoid large distortions
    new_pos2 = pos1 + direction * (np.linalg.norm(pos2 - pos1) + stretch_factor)
    conf.SetAtomPosition(atom2, new_pos2)

    # Rotate the bond around the axis formed by atom1 and atom2 in smaller steps
    for angle in np.arange(0, 2 * np.pi, angle_step):  # Use finer steps for smoother adjustments
        AllChem.SetDihedralDeg(conf, dihedral_atoms[0], dihedral_atoms[1], dihedral_atoms[2], dihedral_atoms[3], np.degrees(angle))

    return mol

def find_distal_atoms(mol):
    """Find the two most distal atoms in a molecule based on the maximum distance."""
    conf = mol.GetConformer()
    max_distance = 0
    distal_atoms = (0, 0)  # Placeholder for the most distal atom pair

    for i in range(mol.GetNumAtoms()):
        for j in range(i + 1, mol.GetNumAtoms()):
            pos_i = np.array(conf.GetAtomPosition(i))
            pos_j = np.array(conf.GetAtomPosition(j))
            distance = np.linalg.norm(pos_j - pos_i)
            if distance > max_distance:
                max_distance = distance
                distal_atoms = (i, j)
    
    return distal_atoms

def align_by_distal_atoms(ligand1, ligand2):
    """Align ligand2 to ligand1 by pulling the distal atoms of ligand2 to match ligand1."""
    conf1 = ligand1.GetConformer()
    conf2 = ligand2.GetConformer()

    # Find the most distal atoms in both ligands
    distal_atoms_lig1 = find_distal_atoms(ligand1)
    distal_atoms_lig2 = find_distal_atoms(ligand2)

    # Get the positions of the distal atoms in ligand1 and ligand2
    pos1_lig1 = np.array(conf1.GetAtomPosition(distal_atoms_lig1[0]))
    pos2_lig1 = np.array(conf1.GetAtomPosition(distal_atoms_lig1[1]))
    pos1_lig2 = np.array(conf2.GetAtomPosition(distal_atoms_lig2[0]))
    pos2_lig2 = np.array(conf2.GetAtomPosition(distal_atoms_lig2[1]))

    # Compute the translation needed to align the distal atoms
    translation_vector1 = pos1_lig1 - pos1_lig2
    translation_vector2 = pos2_lig1 - pos2_lig2

    # Move ligand2's distal atoms to align with ligand1's distal atoms
    for i in range(ligand2.GetNumAtoms()):
        pos = np.array(conf2.GetAtomPosition(i))

        # Linearly interpolate the movement between distal atoms 1 and 2
        t = (np.linalg.norm(pos2_lig1 - pos2_lig2)) / (np.linalg.norm(pos1_lig1 - pos1_lig2) + np.linalg.norm(pos2_lig1 - pos2_lig2))
        new_pos = pos + t * translation_vector1 + (1 - t) * translation_vector2
        conf2.SetAtomPosition(i, new_pos)

    return ligand2

def improve_alignment_with_bond_manipulation(ligand1, ligand2, iterations=5, stretch_factor=0.02, angle_step=np.pi/36, peg_atoms=None):
    """Conservatively stretch and rotate bonds, ensuring the PEG linker is linear, to improve ligand2's alignment."""
    rotatable_bonds = get_rotatable_bonds(ligand2)
    best_rmsd = float('inf')  # Track the best RMSD
    best_ligand2 = None  # Track the best conformation of ligand2

    for i in range(iterations):
        print(f"Iteration {i+1}: Manipulating bonds and optimizing ligand2...")
        for bond in rotatable_bonds:
            manipulate_bonds(ligand2, bond, stretch_factor, angle_step)

        # Align again and calculate RMSD
        rmsd = align_ligands_by_mcs(ligand1, ligand2)
        print(f"RMSD after bond manipulation iteration {i+1}: {rmsd:.4f}")

        # Track the best conformation based on the lowest RMSD
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_ligand2 = Chem.Mol(ligand2)  # Create a copy of the current best ligand2 conformation

    print(f"Best RMSD: {best_rmsd:.4f}")
    return best_ligand2

def align_by_distal_atoms_with_straightening(ligand1, ligand2):
    """Align ligand2 to ligand1 by pulling the distal atoms of ligand2 to match ligand1 and ensure straight PEG linker."""
    conf1 = ligand1.GetConformer()
    conf2 = ligand2.GetConformer()

    # Find the most distal atoms in both ligands
    distal_atoms_lig1 = find_distal_atoms(ligand1)
    distal_atoms_lig2 = find_distal_atoms(ligand2)

    # Get the positions of the distal atoms in ligand1 and ligand2
    pos1_lig1 = np.array(conf1.GetAtomPosition(distal_atoms_lig1[0]))
    pos2_lig1 = np.array(conf1.GetAtomPosition(distal_atoms_lig1[1]))
    pos1_lig2 = np.array(conf2.GetAtomPosition(distal_atoms_lig2[0]))
    pos2_lig2 = np.array(conf2.GetAtomPosition(distal_atoms_lig2[1]))

    # Compute direction vector for the line between distal atoms in ligand1
    direction_vector_lig1 = pos2_lig1 - pos1_lig1
    direction_vector_lig1 /= np.linalg.norm(direction_vector_lig1)

    # Align ligand2's distal atoms to ligand1's distal atoms
    for i in range(ligand2.GetNumAtoms()):
        pos = np.array(conf2.GetAtomPosition(i))
        
        # Linearly interpolate between the distal atoms, pulling towards a straight line
        t = np.dot(pos - pos1_lig2, direction_vector_lig1) / np.linalg.norm(pos2_lig1 - pos1_lig1)
        new_pos = pos1_lig1 + t * direction_vector_lig1 * np.linalg.norm(pos2_lig1 - pos1_lig1)
        conf2.SetAtomPosition(i, new_pos)

    return ligand2

def apply_random_fluctuation(mol, fluctuation_magnitude=0.1):
    """Apply random fluctuations to the atomic positions of the molecule."""
    conf = mol.GetConformer()
    
    for i in range(mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i))
        # Add random fluctuation to each coordinate
        pos += np.random.uniform(-fluctuation_magnitude, fluctuation_magnitude, size=3)
        conf.SetAtomPosition(i, pos)
    
    return mol


def minimize_ligand_to_align_with_substructure(ligand1, ligand2, max_iters=1000):
    """
    Minimize ligand2's geometry to better align with ligand1, focusing on the substructure
    match between the two ligands.
    """
    # Use the force field to optimize the ligand2 structure relative to ligand1
    print(f"Starting minimization to improve substructure alignment...")

    # Check if both ligands have substructure matches
    try:
        rmsd = align_ligands_by_mcs(ligand1, ligand2)
        print(f"Initial RMSD before minimization: {rmsd:.4f}")
    except ValueError as e:
        print(f"Error during MCS alignment: {e}")
        return ligand2

    # Perform MMFF minimization of ligand2
    minimized_ligand2 = minimize_ligand_energy(ligand2, max_iters)

    # Calculate RMSD after minimization
    final_rmsd = align_ligands_by_mcs(ligand1, minimized_ligand2)
    print(f"RMSD after minimization: {final_rmsd:.4f}")

    return minimized_ligand2

# Main flow, loading ligands, and aligning by distal atoms
ligand1 = load_ligand_from_pdb('barrios/P2/md/alkyne/SHP2_modelling/source_ligand.pdb')
ligand2 = load_ligand_from_mol2('barrios/P2/md/alkyne/SHP2_modelling/17-alkyne-PEG8-IA3.mol2')

rmsd_threshold = 2.0  # Set the RMSD threshold to 2 Å
max_iterations = 10   # Maximum iterations
iteration = 0
current_rmsd = float('inf')
fluctuation_magnitude = 0.05  # Fluctuation magnitude

# Variables to track best RMSD and conformation
best_rmsd = float('inf')
best_ligand2 = None

# Iteratively align, refine, and fluctuate until RMSD < 2 Å or max iterations
while current_rmsd >= rmsd_threshold and iteration < max_iterations:
    iteration += 1
    print(f"Iteration {iteration}: Starting alignment and refinement...")
    
    # Align ligand2 to ligand1 by distal atoms and straighten the PEG linker
    ligand2 = align_by_distal_atoms_with_straightening(ligand1, ligand2)

    # Planarize and straighten ligand2
    print("Making ligand2 planar and straight...")
    planar_ligand2 = planarize_ligand(ligand2)
    straight_ligand2 = straighten_bonds(planar_ligand2)

    # Refine the alignment by stretching and rotating bonds
    refined_ligand2 = improve_alignment_with_bond_manipulation(ligand1, straight_ligand2)

    # Apply rigorous minimization for substructure alignment
    refined_ligand2 = minimize_ligand_to_align_with_substructure(ligand1, refined_ligand2)

    # Calculate RMSD after minimization
    current_rmsd = align_ligands_by_mcs(ligand1, refined_ligand2)
    print(f"RMSD after refinement and minimization: {current_rmsd:.4f}")

    # Track the best RMSD and conformation
    if current_rmsd < best_rmsd:
        best_rmsd = current_rmsd
        best_ligand2 = Chem.Mol(refined_ligand2)  # Store the best conformation

    # Check if RMSD is below the threshold
    if current_rmsd < rmsd_threshold:
        print(f"RMSD is below threshold ({current_rmsd:.4f} Å). Stopping refinement.")
        break

    # Apply random fluctuations to ligand2 to escape local minima
    print("Applying random fluctuations to ligand2 structure...")
    ligand2 = apply_random_fluctuation(refined_ligand2, fluctuation_magnitude)

# Perform final minimization on the best refined ligand2 structure
print(f"Best RMSD found: {best_rmsd:.4f}")
print("Performing final minimization on best-aligned ligand2...")
minimized_ligand2 = minimize_ligand_energy(best_ligand2)

# Write the minimized, best-aligned ligand2 to PDB file
print("Writing minimized best-aligned ligand2 to PDB file...")
MolToPDBFile(minimized_ligand2, 'barrios/P2/md/alkyne/SHP2_modelling/minimized_best_refined_ligand2.pdb')
