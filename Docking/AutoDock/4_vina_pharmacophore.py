import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Cluster import Butina
from collections import defaultdict

# Load your data (assuming a CSV file with 'SMILES' and 'Score' columns)
data = pd.read_csv('ligands_scores.csv')

# Convert SMILES to RDKit Molecules
data['Molecule'] = data['SMILES'].apply(Chem.MolFromSmiles)

# Calculate molecular fingerprints
data['Fingerprint'] = data['Molecule'].apply(lambda x: rdMolDescriptors.GetMorganFingerprintAsBitVect(x, 2))

# Cluster molecules based on fingerprints
fps = list(data['Fingerprint'])
dists = []
for i in range(1, len(fps)):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1-x for x in sims])

# Cluster using Butina algorithm
clusters = Butina.ClusterData(dists, len(fps), 0.7, isDistData=True)
data['Cluster'] = -1
for i, cluster in enumerate(clusters):
    for idx in cluster:
        data.at[idx, 'Cluster'] = i

# Analyze clusters
cluster_analysis = defaultdict(list)
for cluster_id in set(data['Cluster']):
    cluster_data = data[data['Cluster'] == cluster_id]
    avg_score = cluster_data['Score'].mean()
    cluster_analysis[cluster_id] = {
        'avg_score': avg_score,
        'size': len(cluster_data),
        'ligands': cluster_data[['SMILES', 'Score']].to_dict('records')
    }

# Print analysis
for cluster_id, info in cluster_analysis.items():
    print(f"Cluster {cluster_id}:")
    print(f"  Average Score: {info['avg_score']}")
    print(f"  Number of Ligands: {info['size']}")
    for ligand in info['ligands']:
        print(f"    {ligand['SMILES']}: {ligand['Score']}")

# Save cluster analysis to a file
with open('cluster_analysis.txt', 'w') as file:
    for cluster_id, info in cluster_analysis.items():
        file.write(f"Cluster {cluster_id}:\n")
        file.write(f"  Average Score: {info['avg_score']}\n")
        file.write(f"  Number of Ligands: {info['size']}\n")
        for ligand in info['ligands']:
            file.write(f"    {ligand['SMILES']}: {ligand['Score']}\n")
