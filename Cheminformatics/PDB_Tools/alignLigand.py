import numpy as np

def parse_pdb(pdb_file):
    """Extract atomic coordinates of the ligand (residue LIG) and protein from a PDB file."""
    ligand_coords = []
    protein_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                residue_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if residue_name == "lig":
                    # Store ligand coordinates
                    ligand_coords.append([x, y, z])
                else:
                    # Store protein atoms for output later
                    protein_lines.append(line)
    return np.array(ligand_coords), protein_lines

def parse_mol2(mol2_file):
    """Extract atomic coordinates of the ligand from a MOL2 file."""
    coords = []
    read_atoms = False
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                read_atoms = True
            elif line.startswith("@<TRIPOS>BOND"):
                read_atoms = False
            elif read_atoms:
                parts = line.split()
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                coords.append([x, y, z])
    return np.array(coords)

def find_mcs(ligand1_coords, ligand2_coords):
    """
    Identify the maximum common substructure (MCS) between ligand1 and ligand2.
    For simplicity, this method assumes that the first N atoms that overlap
    constitute the MCS, where N is the smaller number of atoms in the two ligands.
    """
    # Use the minimum length as a proxy for the MCS size
    mcs_size = min(len(ligand1_coords), len(ligand2_coords))
    return ligand1_coords[:mcs_size], ligand2_coords[:mcs_size]

def kabsch(P, Q):
    """Perform Kabsch algorithm to find the optimal rotation matrix."""
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U

def rmsd(V, W):
    """Calculate the root-mean-square deviation between two sets of vectors."""
    return np.sqrt(np.mean(np.sum((V - W)**2, axis=1)))

def fit_ligands(ligand1_coords, ligand2_coords, max_iter=10000, tolerance=1e-4):
    """Iteratively fit ligand2 to ligand1 using MCS to minimize RMSD."""
    # Identify the maximum common substructure
    mcs_ligand1, mcs_ligand2 = find_mcs(ligand1_coords, ligand2_coords)
    
    # Center the MCS coordinates for both ligands
    mcs_ligand1_center = np.mean(mcs_ligand1, axis=0)
    mcs_ligand2_center = np.mean(mcs_ligand2, axis=0)
    mcs_ligand1 -= mcs_ligand1_center
    mcs_ligand2 -= mcs_ligand2_center
    
    # Compute the optimal rotation using the Kabsch algorithm
    rotation_matrix = kabsch(mcs_ligand2, mcs_ligand1)
    
    # Rotate and translate the entire ligand2 using the calculated transformation
    ligand2_coords = np.dot(ligand2_coords - mcs_ligand2_center, rotation_matrix) + mcs_ligand1_center
    
    # Calculate RMSD between the MCS coordinates
    current_rmsd = rmsd(mcs_ligand1, np.dot(mcs_ligand2, rotation_matrix))
    print(f"Final RMSD after alignment: {current_rmsd:.4f}")

    return ligand2_coords

def save_pdb_with_ligand(protein_lines, ligand_coords, output_file="protein_with_ligand2.pdb"):
    """Save the protein structure along with the aligned ligand2 to a new PDB file, adding a TER card."""
    with open(output_file, 'w') as f:
        # Write the protein atoms
        for line in protein_lines:
            f.write(line)
        
        # Add TER card to separate protein from ligand2
        f.write("TER\n")
        
        # Write the aligned ligand2 coordinates
        for i, (x, y, z) in enumerate(ligand_coords):
            f.write(f"HETATM{i + 1:5d}  C   LIG     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")

# Example usage
ligand1_coords, protein_lines = parse_pdb("C:/Users/ndlev/OneDrive/Documents/Research/barrios/P2/md/alkyne/lowest_energy_struct.pdb")
ligand2_coords = parse_mol2("C:/Users/ndlev/OneDrive/Documents/Research/barrios/P2/md/alkyne/17-alkyne-PEG8-IA3.mol2")

# Fit ligand2 to ligand1 using their MCS
fitted_ligand2_coords = fit_ligands(ligand1_coords, ligand2_coords)

# Save the protein structure along with the aligned ligand2 coordinates
save_pdb_with_ligand(protein_lines, fitted_ligand2_coords, output_file="C:/Users/ndlev/OneDrive/Documents/Research/barrios/P2/md/alkyne/protein_with_aligned_ligand2.pdb")
