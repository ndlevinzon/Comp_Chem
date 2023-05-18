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
    print("CSV Lookup:\n" + str(ligand_lookup))
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

        if fingerprint is not None:
            # Store the fingerprint in the DataFrame
            df.at[index, 'Fingerprint'] = fingerprint

    # Filter out rows with None fingerprints
    df = df[df['Fingerprint'] != '']

    # Return the updated DataFrame
    print("PANDAS Dataframe:\n" + str(df))
    return df


def calculateDistanceMatrix(fingerprints):
    """Calculate Tanimoto similarity matrix between fingerprints"""
    # Get the number of fingerprints
    num_fingerprints = len(fingerprints)
    # Create an empty distance matrix
    distance_matrix = np.zeros((num_fingerprints, num_fingerprints))

    # Iterate over the fingerprints to calculate the distance matrix
    for i in range(1, num_fingerprints):
        for j in range(i):
            # Calculate the Tanimoto similarity between the fingerprints
            similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            # Calculate the distance as 1 minus the similarity
            distance = 1 - similarity
            # Store the distance in both positions of the distance matrix
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    # Return the distance matrix
    print("Distance (Tanimoto) Matrix:\n" + str(distance_matrix))
    return distance_matrix


def growCluster(distance_matrix, pt, neighbors, cluster, assigned, cutoff):
    """Recursively grow a cluster
    Parameters:
        distance_matrix: Tanimoto distance matrix
        pt: current point being processed
        neighbors: distance or similarity values of pt with other points
        cluster: current cluster being grown
        assigned: array indicating if a point has been assigned to a cluster
        cutoff: distance threshold for clustering
    """
    # Convert neighbors to a NumPy array
    neighbors = np.array(neighbors)

    # Find the indices of the unassigned points within the distance threshold
    indices = np.where(np.logical_and(neighbors <= cutoff, np.logical_not(assigned)))[0]

    # Assign the unassigned points to the cluster
    assigned[indices] = 1
    cluster.extend(indices.tolist())

    # Recursively grow the cluster for each unassigned point
    for index in indices:
        growCluster(distance_matrix, index, distance_matrix[index], cluster, assigned, cutoff)


def clusterData(distance_matrix, nPts, cutoff):
    """Cluster a set of data points
    Parameters:
        distance_matrix: distance matrix or pre-calculated distance list
        nPts: number of data points
        cutoff: distance threshold for clustering
    """

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
            growCluster(distance_matrix, pt, distance_matrix[pt], new_cluster, assigned, cutoff)

    return clusters


def clusterFingerprints(distance_matrix, cutoff):
    clusters = clusterData(distance_matrix, len(distance_matrix), cutoff)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def main():
    # Read the CSV and calculate fingerprint representations
    ligand_lookup = readCSV(csv_file='Q9NUW8_lig.csv')
    df = fingerprint(ligand_lookup)

    # Extract RDKit fingerprints from DataFrame
    fingerprints = df['Fingerprint'].tolist()

    # Calculate Tanimoto distance matrix
    distance_matrix = calculateDistanceMatrix(fingerprints)

    # # Generate histogram of Tanimoto distances
    # plt.hist(distance_matrix, bins=50, edgecolor='black', density=True)
    # plt.xlabel('Tanimoto Distance')
    # plt.ylabel('Density')
    # plt.title('Normalized Distribution of Tanimoto Distances')
    # plt.show()

    # Change clustering cutoff and run the clustering procedure for the dataset
    clusters = clusterFingerprints(distance_matrix, cutoff=0.5)
    print(clusters)

    # Give a short report about the numbers of clusters and their sizes
    print("total # clusters: ", len(clusters))
    print("# clusters with only 1 compound: ", sum(1 for c in clusters if len(c) == 1))
    print("# clusters with >5 compounds: ",    sum(1 for c in clusters if len(c) > 5))
    print("# clusters with >25 compounds: ",   sum(1 for c in clusters if len(c) > 25))
    print("# clusters with >100 compounds: ",  sum(1 for c in clusters if len(c) > 100))

    # # Plot the size of the clusters
    # plt.figure(figsize=(15, 4))
    # plt.xlabel("Cluster index")
    # plt.ylabel("Number of molecules")
    # plt.bar(range(1, len(clusters) + 1), [len(c) for c in clusters], lw=5)
    # plt.show()

########################################################################################################################

    # 1. Create a dictionary to map molecule ID numbers to their corresponding cluster number
    cluster_mapping = {}
    cluster_counter = 1

    # 2. Iterate over each group in the 'clusters' list
    for group in clusters:
        # Assign the same cluster number to all molecules in the group
        for molecule_id in group:
            cluster_mapping[molecule_id] = cluster_counter
        cluster_counter += 1

    # 3. Create the 'cluster' column in the DataFrame based on the cluster mapping
    df['Cluster'] = df.index.map(cluster_mapping)

    print(df)

    # 4. Group the DataFrame by the 'Cluster' column
    grouped = df.groupby('Cluster')

    # 5. Iterate over each cluster
    for cluster, group in grouped:
        
        # Exclude clusters with less than 5 members
        if len(group) < 5:
            continue

        # Sort the group by potency in ascending order
        sorted_group = group.sort_values('Potency')

        # Calculate the cumulative sum of molecules with lower potency within the cluster
        cumulative_sum = np.arange(1, len(sorted_group) + 1) / len(sorted_group)

        # Plot the cumulative sum against the -log(potency) values
        plt.plot(-np.log(sorted_group['Potency']), cumulative_sum, label=f'Cluster {cluster}')

    # 6. Set labels for the axes, show the plot
    plt.xlabel('-log(Potency)')
    plt.ylabel('Proportion of Molecules with Improvement within the Same Cluster')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()


##############################################################
#______  _____   _          ______   ______ ______  ________ #
#|  ___ \(____ \ | |        (_____ \ / __   (_____ \(_______/#
#| |   | |_   \ \| |          ____) ) | //| | ____) )  ____  #
#| |   | | |   | | |         /_____/| |// | |/_____/  (___ \ #
#| |   | | |__/ /| |_____    _______|  /__| |_______ _____) )#
#|_|   |_|_____/ |_______)  (_______)\_____/(_______|______/ #
##############################################################    
