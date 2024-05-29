import os
import re
from Bio.PDB import PDBParser
import numpy as np

def get_crystal_ligand_centroid(pdb_file, ligand_resname):
    """
    Calculate the centroid of the ligand from the source PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('source', pdb_file)

    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == ligand_resname:
                    for atom in residue:
                        ligand_atoms.append(atom.get_coord())

    if ligand_atoms:
        ligand_atoms = np.array(ligand_atoms)
        return np.mean(ligand_atoms, axis=0)

    return None

def calculate_centroid(coords):
    """
    Calculate the centroid of a set of coordinates.
    """
    if len(coords) == 0:
        return None
    return np.mean(coords, axis=0)

def parse_docked_poses(pdb_file):
    """
    Parse the PDB file containing docked poses to extract scores and coordinates.
    """
    poses = []
    scores = []
    names = []  # Store the docked pose names
    current_pose = []
    current_score = None
    current_name = None

    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", line)
            if match:
                current_score = float(match.group(1))
        elif line.startswith("REMARK Name = "):
            match = re.search(r'REMARK\s*Name\s*=\s*\s*(\S+)', line)
            if match:
                current_name = line.strip()
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            current_pose.append((atom_name, coords))
        elif line.startswith("ENDMDL"):
            if current_pose:
                poses.append(current_pose)
                scores.append(current_score)
                names.append(current_name)
                current_pose = []
                current_score = None
                current_name = None

    return poses, scores, names

def filter_poses_by_distance(crystal_centroid, poses, scores, names, threshold=10.0):
    """
    Filter poses by distance to the crystal ligand centroid.
    """
    filtered_poses = []
    filtered_scores = []
    filtered_distances = []
    filtered_names = []

    for idx, (pose, score, name) in enumerate(zip(poses, scores, names)):
        pose_coords = np.array([coord for _, coord in pose])
        pose_centroid = calculate_centroid(pose_coords)
        distance = calculate_distance(crystal_centroid, pose_centroid)
        print(f"Pose {idx+1}: Distance from crystal ligand centroid = {distance:.3f} Å")

        if distance <= threshold:
            filtered_poses.append(pose)
            filtered_scores.append(score)
            filtered_distances.append(distance)
            filtered_names.append(name)
            
    return filtered_poses, filtered_scores, filtered_distances, filtered_names

def calculate_distance(centroid, pose_centroid):
    """
    Calculate the distance between two centroids.
    """
    return np.linalg.norm(centroid - pose_centroid)

def write_filtered_poses(output_file, filtered_poses, filtered_scores, filtered_distances, filtered_names, ligand_name, cutoff_distance):
    """
    Write filtered poses to a PDB file with additional information in the header.
    """
    with open(output_file, 'w') as out_f:
        out_f.write(f"REMARK Filtered poses generated using distance cutoff of {cutoff_distance} A from the crystal ligand centroid\n")
        out_f.write(f"REMARK Source ligand: {ligand_name}\n")
        out_f.write("\n")

        for idx, (pose, score, distance, name) in enumerate(zip(filtered_poses, filtered_scores, filtered_distances, filtered_names)):
            out_f.write(f"MODEL {idx+1}\n")
            out_f.write(f"REMARK Name = {name}\n")  # Include the name remark
            out_f.write(f"REMARK Model number: {idx+1}\n")
            out_f.write(f"REMARK Vina score: {score}\n")
            out_f.write(f"REMARK Distance from crystal ligand centroid: {distance:.3f} A\n")
            for atom_name, coord in pose:
                out_f.write(f"HETATM{atom_name:<4} {'':4}{'LIG A   1    '}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00           C\n")
            out_f.write("ENDMDL\n")

def main(source_pdb, ligand_resname, docked_pdb, output_pdb):
    """
    Main function to filter docked poses based on distance to the crystal ligand centroid.
    """
    crystal_centroid = get_crystal_ligand_centroid(source_pdb, ligand_resname)
    if crystal_centroid is None:
        print("Could not find ligand in the source PDB file.")
        return

    poses, scores, names = parse_docked_poses(docked_pdb)
    total_input_poses = len(poses)  # Count of input docked poses
    filtered_poses, filtered_scores, filtered_distances, filtered_names = filter_poses_by_distance(crystal_centroid, poses, scores, names)
    total_filtered_poses = len(filtered_poses)  # Count of filtered docked poses

    output_path = os.path.dirname(source_pdb)
    output_file = os.path.join(output_path, output_pdb)
    write_filtered_poses(output_file, filtered_poses, filtered_scores, filtered_distances, filtered_names, ligand_resname, cutoff_distance=10.0)
    
    # Print the counts to the console
    print(f"Total docked poses inputted: {total_input_poses}")
    print(f"Total docked poses written to output: {total_filtered_poses}")
    print(f"Filtered poses written to: {output_file}")

if __name__ == "__main__":
    source_pdb = "C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test/6IRT_lig.pdb"
    ligand_resname = "AUU"  # Specify the residue name of the source ligand
    docked_pdb = "C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test/best_docked.pdb"
    output_pdb = "rmsd_sort.pdb"
    main(source_pdb, ligand_resname, docked_pdb, output_pdb)
