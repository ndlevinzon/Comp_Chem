import os
import re

def extract_scores_and_poses(pdb_file):
    """
    Extract all docking scores and their corresponding poses from a PDBQT file.
    """
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    scores_and_poses = []
    score = None
    pose = []
    name_remark = None

    for line in lines:
        if line.startswith("REMARK VINA RESULT:") or line.startswith("REMARK  Name = "):
            if score is not None and pose:
                scores_and_poses.append((score, name_remark, pose))
            if line.startswith("REMARK VINA RESULT:"):
                match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", line)
                if match:
                    score = float(match.group(1))
            elif line.startswith("REMARK  Name = "):
                name_remark = line.strip()
            pose = []
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ENDMDL"):
            pose.append(line)
    
    if score is not None and pose:
        scores_and_poses.append((score, name_remark, pose))

    return scores_and_poses

def combine_poses_to_pdb(output_pdb_file, sorted_poses):
    """
    Combine sorted poses into a single PDB file, including REMARK lines containing "Name = ".
    """
    with open(output_pdb_file, 'w') as output_file:
        for idx, (score, name_remark, pose) in enumerate(sorted_poses):
            output_file.write(f"MODEL {idx+1}\n")
            if name_remark:
                output_file.write(f"{name_remark}\n")
            output_file.write(f"REMARK VINA RESULT: {score}\n")
            output_file.writelines(line for line in pose if not line.startswith("REMARK VINA RESULT:") and not line.startswith("REMARK Name = "))
            output_file.write("ENDMDL\n")


def main(input_directory):
    """
    Main function to combine PDBQT files into a single sorted PDB file.
    """
    pdbqt_files = [f for f in os.listdir(input_directory) if f.endswith('.pdbqt')]

    all_scores_and_poses = []
    for pdbqt_file in pdbqt_files:
        file_path = os.path.join(input_directory, pdbqt_file)
        scores_and_poses = extract_scores_and_poses(file_path)
        all_scores_and_poses.extend(scores_and_poses)

    sorted_poses = sorted(all_scores_and_poses, key=lambda x: x[0])

    output_pdb_file = os.path.join(input_directory, "best_docked.pdb")
    combine_poses_to_pdb(output_pdb_file, sorted_poses)
    print(f"Combined PDB file created: {output_pdb_file}")

if __name__ == "__main__":
    input_directory = "C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/test"
    main(input_directory)
