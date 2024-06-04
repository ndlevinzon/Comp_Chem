import re
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform

# Function to extract SMILES codes and Vina scores from PDB file
def extract_smiles_and_scores_from_pdb(file_path):
    smiles_scores_list = []
    current_smiles = None
    current_score = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('REMARK'):
                smiles_match = re.search(r'SMILES\s*=\s*(\S+)', line)
                score_match = re.search(r'Vina score:\s*(-?\d+\.\d+)', line)
                if smiles_match:
                    current_smiles = smiles_match.group(1)
                if score_match:
                    current_score = float(score_match.group(1))
                    if current_smiles is not None:
                        smiles_scores_list.append((current_smiles, current_score))
                        current_smiles = None
                        current_score = None
    return smiles_scores_list

# Function to convert SMILES to fingerprints
def smiles_to_fingerprints(smiles_list):
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            arr = np.zeros((1,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fingerprints.append(arr)
    return np.array(fingerprints)

# Function to perform clustering and write output to file
def perform_clustering_and_write_output(smiles_scores, output_file):
    smiles_list = [smiles for smiles, _ in smiles_scores]
    fingerprints = smiles_to_fingerprints(smiles_list)
    
    # Calculate the distance matrix
    distance_matrix = squareform(pdist(fingerprints, metric='jaccard'))
    
    # Perform Agglomerative Clustering
    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=0.5, metric='precomputed', linkage='average')
    clusters = clustering.fit_predict(distance_matrix)
    
    # Organize smiles codes by cluster and sort by Vina scores
    cluster_dict = {}
    for (smiles, score), cluster in zip(smiles_scores, clusters):
        if cluster not in cluster_dict:
            cluster_dict[cluster] = []
        cluster_dict[cluster].append((smiles, score))
    
    for cluster in cluster_dict:
        cluster_dict[cluster].sort(key=lambda x: x[1])
    
    # Write clusters to output file and prepare images for each cluster
    cluster_images = []
    legends = []
    with open(output_file, 'w') as f:
        for cluster in sorted(cluster_dict.keys()):
            f.write(f'Cluster {cluster}:\n')
            for smiles, score in cluster_dict[cluster]:
                f.write(f'SMILES: {smiles}, Vina score: {score}\n')
            f.write('\n')
            # Get the SMILES code with the most negative Vina score in the cluster
            best_smiles = cluster_dict[cluster][0][0]
            best_score = cluster_dict[cluster][0][1]
            mol = Chem.MolFromSmiles(best_smiles)
            if mol:
                cluster_images.append(mol)
                legends.append(f'Cluster {cluster}\nSMILES: {best_smiles}\nVina score: {best_score}')
    
    # Draw grid of cluster images with legends
    img = Draw.MolsToGridImage(cluster_images, molsPerRow=4, subImgSize=(300, 300), legends=legends)
    img.save('C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/dissimilarity/cluster_molecules.png')

# Path to the PDB file
pdb_file_path = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/dissimilarity/lig_HTVS_sort_filter.pdb'
output_file_path = 'C:/Users/ndlev/OneDrive/Documents/Research/babst/htvs/dissimilarity/output_clusters.txt'

# Extract SMILES codes and Vina scores
smiles_scores = extract_smiles_and_scores_from_pdb(pdb_file_path)

# Perform clustering and write output
perform_clustering_and_write_output(smiles_scores, output_file_path)
