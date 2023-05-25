# TODO: lower cutoffs for clustering and use the biggest cluster;
# TODO: Output structures from the most populated clusters


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors, Draw
from scipy.stats import t
from scipy.optimize import curve_fit


def readCSV(csv_file):
    """Read CSV file from ChEMBL Query and create a lookup table"""
    data = pd.read_csv(csv_file, lineterminator='\n', sep=';')

    # Create a dictionary to store the SMILES and POTENCY_VALUE
    ligand_lookup = {}

    for _, row in data.iterrows():
        SMILES = row['Smiles']
        POTENCY= row['Standard Value']
        ligand_lookup[SMILES] = POTENCY

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
    df = pd.DataFrame(list(ligand_lookup.items()), columns=['SMILES', 'POTENCY'])

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

########################################################################################################################

def calculate_linkage_distance(cluster1, cluster2, distance_matrix):
    """Calculate the linkage distance between two clusters using complete linkage"""
    distances = distance_matrix[np.ix_(cluster1, cluster2)]
    return np.max(distances)


def merge_clusters(clusters, i, j):
    """Merge two clusters"""
    clusters[i].extend(clusters[j])
    del clusters[j]


def clusterData(distance_matrix, nPts, cutoff):
    """Cluster a set of data points using complete linkage clustering
    Parameters:
        distance_matrix: distance matrix or pre-calculated distance list
        nPts: number of data points
        cutoff: distance threshold for clustering
    """

    # Create a list of clusters
    clusters = [[i] for i in range(nPts)]

    while True:
        max_distance = 0
        merge_indices = None

        # Find the pair of clusters with the maximum distance
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                distance = calculate_linkage_distance(clusters[i], clusters[j], distance_matrix)
                if distance > max_distance:
                    max_distance = distance
                    merge_indices = (i, j)

        # Stop if the maximum distance is below the cutoff
        if max_distance <= cutoff:
            break

        # Merge the clusters with the maximum distance
        merge_clusters(clusters, merge_indices[0], merge_indices[1])

    return clusters

def clusterFingerprints(distance_matrix, cutoff):
    clusters = clusterData(distance_matrix, len(distance_matrix), cutoff)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters

########################################################################################################################

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

    # Convert units of potency from nM to M
    df['POTENCY_Converted'] = df['POTENCY'] * 1e-9

    # Create a dictionary to map molecule ID numbers to their corresponding cluster number
    cluster_mapping = {}
    cluster_counter = 1

    # Iterate over each group in the 'clusters' list
    for group in clusters:
        # Assign the same cluster number to all molecules in the group
        for molecule_id in group:
            cluster_mapping[molecule_id] = cluster_counter
        cluster_counter += 1

    # Create the 'cluster' column in the DataFrame based on the cluster mapping
    df['CLUSTER'] = df.index.map(cluster_mapping)

    # Group the DataFrame by the 'Cluster' column
    grouped = df.groupby('CLUSTER')

    # Print PANDAS DF
    print(df)

    # Iterate over each cluster
    for cluster, group in grouped:

        # Exclude clusters with less than 5 members
        if len(group) < 5:
            continue

        # Drop entries with NaN potency values
        group = group.dropna(subset=['POTENCY_Converted'])

        # Sort the group by potency in ascending order
        sorted_group = group.sort_values('POTENCY_Converted')

        # Calculate the proportion for each molecule in the cluster
        proportion = np.arange(1, len(sorted_group) + 1) / len(sorted_group)

        # Assign the proportion values to the 'Proportion' column
        df.loc[sorted_group.index, 'Proportion'] = proportion

    # Exclude entries with NaN values from subsequent calculations
    df = df.dropna()

    # Perform linear regression using all data points
    x_data = df['POTENCY_Converted']
    y_data = df['Proportion']
    slope, intercept = np.polyfit(x_data, y_data, deg=1)
    x_fit = np.linspace(min(x_data), max(x_data), 100)
    y_fit_linear = slope * x_fit + intercept

    # Estimate initial parameters for sigmoid function
    x_min = min(x_data)
    a_init = 1
    b_init = x_min
    c_init = 1

    # Fit a sigmoid function to the data
    def sigmoid(x, a, b, c):
        return c / (1 + np.exp(-a * (x - b)))

    params, _ = curve_fit(sigmoid, x_data, y_data, p0=[a_init, b_init, c_init])
    a, b, c = params
    y_fit_sigmoid = sigmoid(x_fit, a, b, c)

    # Calculate standard error of the estimate
    y_pred = slope * x_data + intercept
    error = y_data - y_pred
    mse = np.mean(error ** 2)
    se = np.sqrt(mse)
    
    # Calculate upper and lower bounds of the confidence band
    confidence = 0.95
    n = len(x_data)
    t_value = t.ppf((1 + confidence) / 2, n - 2)
    confidence_band = t_value * se


    # Set the x-axis values based on isLogarithmic variable
    isLogarithmic = True
    if isLogarithmic:
        x_data_plot = -np.log10(x_data)
        x_fit_plot = -np.log10(x_fit)
        x_label = '-log(Potency) (M)'
    else:
        x_data_plot = x_data
        x_fit_plot = x_fit
        x_label = 'Potency (M)'

    # Create a new figure and axes for the graph
    fig, ax_graph = plt.subplots()

    # Plot the cumulative sum against the potency values
    ax_graph.scatter(x_data_plot, y_data, c=df['CLUSTER'], label='Clusters')

    # Plot the regression line based on isSigmoid variable
    isSigmoid = True
    if isSigmoid:
        # Plot the sigmoid curve and include the formula for the sigmoid curve
        formula: str = f'Sigmoid Curve: y = {c:.3f} / (1 + exp(-{a:.3f} * (x - {b:.3f})))'
        ax_graph.plot(x_fit_plot, y_fit_sigmoid, 'r--', label=formula)

        # Calculate correlation coefficient
        correlation_coefficient = np.corrcoef(x_data, y_data)[0, 1]

        ax_graph.legend(title=f'Correlation Coefficient: {correlation_coefficient:.3f}')
    else:
        # Plot the linear regression line and include the formula for the curve
        ax_graph.plot(x_fit_plot, y_fit_linear, 'r--', label=f'Linear Regression: y = {slope:.3g}x + {intercept:.3g}')

        # Calculate the confidence curve for linear regression
        y_fit_upper = y_fit_linear + confidence_band
        y_fit_lower = y_fit_linear - confidence_band

        # Plot the confidence curve for regression
        ax_graph.plot(x_fit_plot, y_fit_upper, 'b--', label='95% Confidence')
        ax_graph.plot(x_fit_plot, y_fit_lower, 'b--')

        # Calculate correlation coefficient
        correlation_coefficient = np.corrcoef(x_data, y_data)[0, 1]

        # Include the correlation coefficient in the legend
        ax_graph.legend(title=f'Correlation Coefficient: {correlation_coefficient:.3f}')

    # Set the x-axis limits for the graph
    ax_graph.set_xlim(min(x_data_plot), max(x_data_plot))

    # Set labels for the axes of the graph
    ax_graph.set_xlabel(x_label)
    ax_graph.set_ylabel('Prop. Molecules with Improvement within the Same Cluster')

    # Show the plot
    plt.show()

    # Save the plot as a PNG file
    output_file = "output_graph.png"
    plt.savefig(output_file, dpi=300)
    print(f"Graph saved as {output_file}")


if __name__ == '__main__':
    main()
