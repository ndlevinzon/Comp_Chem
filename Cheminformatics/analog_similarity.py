import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors


def readCSV(csv_file):
    """Read CSV file from ChEMBL Query"""
    data = pd.read_csv(csv_file, lineterminator='\n', sep=';', dtype={13: str})

    # Iterate over each row in the CSV to enumerate Actives in a lookup table
    ligand_lookup = {}
    for _, row in data.iterrows():
        # Store the values of row[7] (SMILES) and row[10] (POTENCY_VALUE) in the lookup table
        ligand_lookup[row[7]] = row[10]
    return ligand_lookup


def generateFingerprint(smiles):
    """Generate RDKit fingerprint for a given SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return fingerprint
    return None


def fingerprint(ligand_lookup):
    """Calculate fingerprints from ChEMBL SMILES"""
    df = pd.DataFrame(list(ligand_lookup.items()), columns=['SMILES', 'Potency'])
    # Create a new column to store the fingerprints
    df['Fingerprint'] = ''

    for index, row in df.iterrows():
        # Get the SMILES code from the row
        smiles = str(row['SMILES'])

        # Generate the RDKit fingerprint
        fingerprint = generateFingerprint(smiles)

        # Store the fingerprint in the DataFrame
        df.at[index, 'Fingerprint'] = fingerprint

    # Return the updated DataFrame
    return df


def tanimotoDistanceMatrix(fp_list):
    """Calculate distance matrix for RDKit fingerprints"""
    num_structures = len(fp_list)
    dissimilarity_matrix = np.zeros((num_structures, num_structures))

    for i in range(1, num_structures):
        for j in range(i):
            if fp_list[i] is not None and fp_list[j] is not None:
                similarity = DataStructs.TanimotoSimilarity(fp_list[i], fp_list[j])
                dissimilarity = 1 - similarity
                dissimilarity_matrix[i][j] = dissimilarity
                dissimilarity_matrix[j][i] = dissimilarity

    return dissimilarity_matrix


def clusterData(data, nPts, distThresh, isDistData=False):
    """Cluster a set of data points
    Parameters:
        data: distance matrix or pre-calculated distance list
        nPts: number of data points
        distThresh: distance threshold for clustering
        isDistData: True if data is a distance matrix, False if data is a pre-calculated distance list
    Returns:
        list of clusters, where each cluster is represented by a list of indices
    """
    if not isDistData:
        # Convert pre-calculated distance list to a distance matrix
        distance_matrix = np.zeros((nPts, nPts))
        tri_indices = np.triu_indices(nPts, k=1)
        distance_matrix[tri_indices] = data
        distance_matrix += distance_matrix.T

        # Assign the distance matrix to the data variable
        data = distance_matrix

    # Create a list of clusters
    clusters = []
    # Create a list to store whether a data point has been assigned to a cluster
    assigned = np.zeros(nPts, bool)

    for pt in range(nPts):
        if not assigned[pt]:
            # Add a new cluster
            new_cluster = [pt]
            clusters.append(new_cluster)
            assigned[pt] = 1
            # Grow the cluster
            growCluster(pt, data[pt], new_cluster, assigned, distThresh)

    return clusters


def calculateDistanceMatrix(fingerprints):
    num_fingerprints = len(fingerprints)
    distance_matrix = np.zeros((num_fingerprints, num_fingerprints))

    for i in range(1, num_fingerprints):
        for j in range(i):
            similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            distance = 1 - similarity
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    return distance_matrix


def clusterFingerprints(fingerprints, cutoff=0.2):
    distance_matrix = calculateDistanceMatrix(fingerprints)
    clusters = clusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def growCluster(pt, neighbors, cluster, assigned, distThresh):
    """Recursively grow a cluster
    Parameters:
        pt: current point being processed
        neighbors: distance or similarity values of pt with other points
        cluster: current cluster being grown
        assigned: array indicating if a point has been assigned to a cluster
        distThresh: distance threshold for clustering
    """
    # Convert neighbors to a NumPy array
    neighbors = np.array(neighbors)

    # Find the indices of the unassigned points within the distance threshold
    indices = np.where(np.logical_and(neighbors <= distThresh, np.logical_not(assigned)))[0]

    # Assign the unassigned points to the cluster
    assigned[indices] = 1
    cluster.extend(indices.tolist())

    # Recursively grow the cluster for each unassigned point
    for index in indices:
        growCluster(index, neighbors[index], cluster, assigned, distThresh)


def main():
    # Read the CSV and calculate fingerprint representations
    ligand_lookup = readCSV(csv_file='Q9NUW8_lig.csv')
    df = fingerprint(ligand_lookup)
    print("CSV Reading Complete! Generating Fingerprints...")

    # Extract RDKit fingerprints from DataFrame
    fingerprints = df['Fingerprint'].tolist()
    print("Fingerprint Generation Complete! Calculating Distance Matrix...")

    # Calculate Tanimoto distance matrix
    dissimilarity_matrix = tanimotoDistanceMatrix(fingerprints)
    print("Tanimoto Distance Matrix Complete! Generating Histogram...")

    # Generate histogram of Tanimoto distances
    plt.hist(dissimilarity_matrix, bins=50, edgecolor='black')
    plt.xlabel('Tanimoto Distance')
    plt.ylabel('Frequency')
    plt.title('Distribution of Tanimoto Distances')
    plt.show()

    # Run the clustering procedure for the dataset
    clusters = clusterFingerprints(fingerprints, cutoff=0.3)
    print("Clustering Complete!")

    # Give a short report about the numbers of clusters and their sizes
    num_clust_g1 = sum(1 for c in clusters if len(c) == 1)
    num_clust_g5 = sum(1 for c in clusters if len(c) > 5)
    num_clust_g25 = sum(1 for c in clusters if len(c) > 25)
    num_clust_g100 = sum(1 for c in clusters if len(c) > 100)

    print("total # clusters: ", len(clusters))
    print("# clusters with only 1 compound: ", num_clust_g1)
    print("# clusters with >5 compounds: ", num_clust_g5)
    print("# clusters with >25 compounds: ", num_clust_g25)
    print("# clusters with >100 compounds: ", num_clust_g100)

    # Plot the size of the clusters
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.set_xlabel("Cluster index")
    ax.set_ylabel("Number of molecules")
    ax.bar(range(1, len(clusters) + 1), [len(c) for c in clusters], lw=5)


if __name__ == '__main__':
    main()
#
############################################
#( \( )(  _ \(  )    (__ \ / _ \(__ \ (__ )#
# )  (  )(_) ))(__    / _/( (_) )/ _/  (_ \#
#(_)\_)(____/(____)  (____)\___/(____)(___/#
############################################
