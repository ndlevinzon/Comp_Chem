from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolfiles import MolToPDBFile
import numpy as np

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

def improve_alignment_with_bond_manipulation(ligand1, ligand2, iterations=5, stretch_factor=0.02, angle_step=np.pi/36):
    """Conservatively stretch and rotate bonds, including PEG linker, to improve ligand2's alignment to ligand1."""
    rotatable_bonds = get_rotatable_bonds(ligand2)

    for i in range(iterations):
        print(f"Iteration {i+1}: Manipulating bonds (including PEG flexibility) and optimizing ligand2...")
        for bond in rotatable_bonds:
            # Stretch and rotate each bond, with smaller steps and a focus on the PEG linker
            manipulate_bonds(ligand2, bond, stretch_factor, angle_step)

        # Optimize geometry after manipulation
        optimize_ligand_geometry(ligand2)

        # Align again and calculate RMSD
        rmsd = align_ligands_by_mcs(ligand1, ligand2)
        print(f"RMSD after bond manipulation iteration {i+1}: {rmsd:.4f}")

    return ligand2


def load_ligand_from_pdb(pdb_file):
    print(f"Loading ligand1 from PDB: {pdb_file}")
    return Chem.rdmolfiles.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)

def load_ligand_from_mol2(mol2_file):
    print(f"Loading ligand2 from MOL2: {mol2_file}")
    return Chem.MolFromMol2File(mol2_file, removeHs=False)

# Load ligands
ligand1 = load_ligand_from_pdb('barrios/P2/md/alkyne/SHP2_modelling/source_ligand.pdb')
ligand2 = load_ligand_from_mol2('barrios/P2/md/alkyne/SHP2_modelling/13-alkyne-PEG8-IA3.mol2')

# Align ligand2 to ligand1
try:
    rmsd = align_ligands_by_mcs(ligand1, ligand2)
except ValueError as e:
    print(f"Error in alignment: {e}")

# Optimize the geometry of ligand2 after alignment
optimized_ligand2 = optimize_ligand_geometry(ligand2)

# Planarize and straighten ligand2
print("Making ligand2 planar and straight...")
planar_ligand2 = planarize_ligand(optimized_ligand2)
straight_ligand2 = straighten_bonds(planar_ligand2)

# Refine the alignment by stretching and rotating bonds
refined_ligand2 = improve_alignment_with_bond_manipulation(ligand1, straight_ligand2)

# Write the refined ligand2 to PDB file
print("Writing refined ligand2 to PDB file...")
MolToPDBFile(refined_ligand2, 'barrios/P2/md/alkyne/SHP2_modelling/refined_ligand2.pdb')

# Perform energy minimization on the refined ligand2
minimized_ligand2 = minimize_ligand_energy(refined_ligand2)

# Write minimized ligand2 to PDB file
print("Writing minimized ligand2 to PDB file...")
MolToPDBFile(minimized_ligand2, 'barrios/P2/md/alkyne/SHP2_modelling/minimized_ligand2.pdb')

print("Process completed successfully.")
