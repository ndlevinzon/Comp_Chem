import heapq
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import stats
from scipy.sparse import lil_matrix
from scipy.optimize import curve_fit


def read_csv(csv_file):
    """
    Read CSV file from ChEMBL Query and create a pandas DataFrame

    Parameters:
    csv_file (str): Path to the CSV file

    Returns:
    df (pandas.DataFrame): DataFrame containing the data from the CSV file
    """

    data = pd.read_csv(csv_file, lineterminator='\n', sep=';', dtype=str)                # Read the CSV file into a pandas DataFrame
    df = pd.DataFrame({'SMILES': data['Smiles'], 'POTENCY': data['Standard Value']})     # Create a new DataFrame with 'SMILES' and 'POTENCY' columns from the CSV data
    df.dropna(subset=['SMILES', 'POTENCY'], inplace=True)                                # Drop rows with missing values in 'SMILES' and 'POTENCY' columns
    df['POTENCY'] = pd.to_numeric(df['POTENCY'], errors='coerce')                        # Convert the 'POTENCY' column to numeric type
    df['POTENCY_Converted'] = df['POTENCY'] * 1e-9                                       # Convert the potency units from nM to M by multiplying by 1e-9

    return df


def generate_fingerprint(smiles):
    """Generate RDKit fingerprint for a given SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        new_fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return new_fingerprint
    return None


def fingerprint(df):
    """Calculate fingerprints from ChEMBL SMILES"""

    df['Fingerprint'] = ''                                                                 # Create a new column to store the fingerprints
    for index, row in df.iterrows():
        smiles = str(row['SMILES'])                                                        # Get the SMILES code from the row
        fingerprint = generate_fingerprint(smiles)                                         # Generate the RDKit fingerprint

        if fingerprint is not None:
            df.at[index, 'Fingerprint'] = fingerprint                                      # Store the fingerprint in the DataFrame

    df = df[df['Fingerprint'] != '']                                                       # Filter out rows with None fingerprints

    print("PANDAS Dataframe:\n" + str(df))                                                 # Return the updated DataFrame
    return df


def calculate_distance_matrix(fingerprints):
    """Calculate Tanimoto similarity matrix between fingerprints"""
    num_fingerprints = len(fingerprints)                                                   # Get the number of fingerprints
    distance_matrix = lil_matrix((num_fingerprints, num_fingerprints))                     # Create an empty sparse distance matrix

    # Iterate over the fingerprints to calculate the distance matrix
    for i in range(1, num_fingerprints):
        for j in range(i):
            similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])  # Calculate the Tanimoto similarity between the fingerprints
            distance = 1 - similarity                                                      # Calculate the distance as 1 minus the similarity
            # Store the distance in both positions of the distance matrix
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance
            
    return distance_matrix.tocsr()                                                         # Convert to a Compressed Sparse Row (CSR) matrix# Convert to a Compressed Sparse Row (CSR) matrix


def calculate_linkage_distance(cluster1, cluster2, distance_matrix):
    """Calculate the linkage distance between two clusters using complete linkage"""
    distances = distance_matrix[cluster1][:, cluster2]
    return np.max(distances)


def merge_clusters(clusters, index1, index2):
    """Merge two clusters in the list of clusters"""
    if index1 >= len(clusters) or index2 >= len(clusters):
        return
    
    clusters[index1].extend(clusters[index2])
    del clusters[index2]


def cluster_data(distance_matrix, nPts, cutoff):
    """Cluster a set of data points using complete linkage clustering
    Parameters:
        distance_matrix: distance matrix or pre-calculated distance list
        nPts: number of data points
        cutoff: distance threshold for clustering
    """

    # Create a list of clusters and a max heap of distances
    clusters = [[i] for i in range(nPts)]
    max_heap = []
    for i in range(len(clusters)):
        for j in range(i + 1, len(clusters)):
            distance = calculate_linkage_distance(clusters[i], clusters[j], distance_matrix)
            heapq.heappush(max_heap, (-distance, (i, j)))

    while len(max_heap) > 0:
        max_distance, merge_indices = heapq.heappop(max_heap)
        
        if len(max_heap) == 0:                                                          # Check if the heap is empty
            break
            
        if -max_distance <= cutoff:                                                     # Stop if the maximum distance is below the cutoff
            break

        merge_clusters(clusters, merge_indices[0], merge_indices[1])                    # Merge the clusters with the maximum distance
    return clusters


def cluster_fingerprints(distance_matrix, cutoff):
    clusters = cluster_data(distance_matrix, len(distance_matrix), cutoff)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def add_scaffold_and_group(df):
    """Add 'SCAFFOLD' and 'GROUP' columns to the DataFrame"""

    # Add 'SCAFFOLD' column to the DataFrame
    df['SCAFFOLD'] = df['SMILES'].apply(
        lambda x: Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(str(x)))))

    # Initialize dictionary and counter for grouping
    group_mapping = {}
    group_counter = 1

    def assign_group(scaffold):
        nonlocal group_counter

        # Check if scaffold is already in the group mapping
        if scaffold in group_mapping:
            return group_mapping[scaffold]                                             # Return the existing group number
        else:
            group_mapping[scaffold] = group_counter                                    # Add the scaffold to the group mapping and assign a new group number
            group_counter += 1
            return group_mapping[scaffold]                                             # Return the newly assigned group number for the scaffold

    df['GROUP'] = df['SCAFFOLD'].apply(assign_group)                                   # Add 'GROUP' column to the DataFrame by applying the assign_group function to each scaffold

    # Print the total number of groups found
    print("Total number of groups found:", group_counter - 1)
    
    # Print the number of groups with more than 1 member
    print("Number of groups with more than 1 member:",                       
          len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 1]))

    # Print the number of groups with more than 5 members
    print("Number of groups with more than 5 members:",
          len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 5]))

    # Print the number of groups with more than 10 members
    print("Number of groups with more than 10 members:",
          len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 10]))

    # Print the number of groups with more than 25 members
    group_counts = df['GROUP'].value_counts()
    groups_over_25 = group_counts[group_counts > 25].index.tolist()
    print("Number of groups with more than 25 members:", len(groups_over_25))

    # Print the scaffold and group index for all groups over 25 members
    for group in groups_over_25:
        scaffold = df[df['GROUP'] == group]['SCAFFOLD'].iloc[0]
        print(f"Group {group}: Scaffold - {scaffold}")
    return df                                                                         # Return the modified DataFrame with added 'SCAFFOLD' and 'GROUP' columns


def rank_entries(df):
    """Rank the entries within each group based on POTENCY values"""
    
    df['RANK'] = df.groupby('GROUP')['POTENCY'].rank(ascending=True, method='min')    # Rank the entries within each group based on POTENCY values
    group_counts = df['GROUP'].value_counts().reset_index()                           # Count the number of entries in each group
    group_counts.columns = ['GROUP', 'COUNT']

    df = pd.merge(df, group_counts, on='GROUP')                                       # Merge the group counts back into the DataFrame
    df['PROP_IMPROV'] = df['RANK'] / df['COUNT']                                      # Calculate the proportion of molecules with improvement within the same scaffold

    return df                                                                         # Return the modified DataFrame with added 'RANK', 'COUNT', and 'PROP_IMPROV' columns


def plot(df, isLogarithmic):
    """Plot the clusters and logarithmic curve"""
    plt.figure()
    
    # Fit a logarithmic function
    def log_func(x, a, b):
    """Logarithmic function"""
        return a * np.log(x) + b

    # Fit a sigmoid function
    def sigmoid(x, a, b, c):
        return c / (1 + np.exp(-a * (x - b)))

    y_data = df['PROP_IMPROV']
    y_label = 'Prop. Molecules with Improvement within the Same Scaffold'
    x_data = df['POTENCY_Converted']

    if isLogarithmic:
        x_label = '-log(Potency) (M)'
        x_data = -np.log10(x_data)  # Convert x_data to -log(x)
        plt.xscale("log")  # Scale x-axis to logarithmic if isLogarithmic is True
    else:
        x_label = 'log(Potency) (M)'

    # Select the top 5 groups with the most members
    group_counts = df['GROUP'].value_counts()
    top_groups = group_counts[group_counts.index != ''].nlargest(5).index.tolist()
    print("Top 5 Most Populous Groups:", top_groups)

    # Combine data points from the top groups
    combined_x = []
    combined_y = []
    for group in top_groups:
        group_data = df[df['GROUP'] == group]
        combined_x.extend(x_data[group_data.index])
        combined_y.extend(y_data[group_data.index])

        # Each group will have a unique color
        color = plt.cm.Set1(group % plt.cm.Set1.N)
        plt.scatter(x_data[group_data.index], y_data[group_data.index], color=color, alpha=0.6, s=10)

        # Print the scaffold for the group
        scaffold = group_data['SCAFFOLD'].iloc[0]
        print(f"Group {group}: Scaffold - {scaffold}")

        lowest_potency_index = group_data['POTENCY'].idxmin()
        lowest_potency_smiles = group_data.loc[lowest_potency_index, 'SMILES']
        print(f"Group: {group}, Lowest POTENCY SMILES: {lowest_potency_smiles}")

    # Convert combined_x and combined_y to NumPy arrays
    combined_x = np.array(combined_x)
    combined_y = np.array(combined_y)

    # Fit the logarithmic curve to the combined data
    popt_log, _ = curve_fit(log_func, combined_x, combined_y)

    # Generate the line coordinates for the logarithmic fit
    line_x_log = np.linspace(min(combined_x), max(combined_x), 100)
    line_y_log = log_func(line_x_log, *popt_log)

    # Calculate R-squared value for logarithmic fit
    residuals_log = combined_y - log_func(combined_x, *popt_log)
    ss_residuals_log = np.sum(residuals_log ** 2)
    ss_total_log = np.sum((combined_y - np.mean(combined_y)) ** 2)
    r_squared_log = 1 - (ss_residuals_log / ss_total_log)
    r_squared_text_log = f"R-squared (Logarithmic): {r_squared_log:.2f}"

    # Plot the logarithmic fit and its legend entry
    plt.plot(line_x_log, line_y_log, color='black', linestyle='--',
             label=f'Logarithmic Fit: {popt_log[0]:.3f} * log(x) + {popt_log[1]:.3f}')
    plt.annotate(r_squared_text_log, xy=(0, 0.92), xycoords='axes fraction', ha='left', va='top')

    # Fit the line of best fit to the combined data
    slope, _, r_value, p_value, std_err = stats.linregress(combined_x, combined_y)

    # Generate the line coordinates for the line of best fit
    line_x_lin = np.linspace(min(combined_x), max(combined_x), 100)
    line_y_lin = slope * combined_x

    # Convert x-axis to logarithmic scale and invert the sign
    plt.xscale("log")

    # Calculate residuals and standard error of estimate
    residuals = combined_y - (slope * combined_x)
    s_residual = np.sqrt(np.sum(residuals ** 2) / (len(combined_x) - 2))

    # Calculate confidence interval bounds
    t_value = stats.t.ppf(0.975, len(combined_x) - 2)
    conf_interval = t_value * s_residual * np.sqrt(
        1 / len(combined_x) + (line_x_lin - np.mean(combined_x)) ** 2 / np.sum((combined_x - np.mean(combined_x)) ** 2))

    # Calculate confidence band bounds
    conf_band = t_value * s_residual * np.sqrt(
        1 + 1 / len(combined_x) + (line_x_lin - np.mean(combined_x)) ** 2 / np.sum(
            (combined_x - np.mean(combined_x)) ** 2))

    # Plot the line of best fit and confidence band
    plt.plot(line_x_lin, line_y_lin, color='red', linestyle='--', label='Line of Best Fit')
    plt.fill_between(line_x_lin, line_y_lin - conf_interval, line_y_lin + conf_interval, color='gray', alpha=0.2)
    plt.fill_between(line_x_lin, line_y_lin - conf_band, line_y_lin + conf_band, color='gray', alpha=0.1)

    # Calculate R-squared value for linear fit
    residuals_lin = combined_y - (slope * combined_x)
    ss_residuals_lin = np.sum(residuals_lin ** 2)
    ss_total_lin = np.sum((combined_y - np.mean(combined_y)) ** 2)
    r_squared_lin = 1 - (ss_residuals_lin / ss_total_lin)
    r_squared_text_lin = f"R-squared (Linear): {r_squared_lin:.2f}"

    # Plot the line of best fit and its legend entry
    plt.plot(line_x_lin, line_y_lin, color='red', linestyle='--',
             label=f'Line of Best Fit: y = {slope:.3f}x')
    plt.annotate(r_squared_text_lin, xy=(0, 0.85), xycoords='axes fraction', ha='left', va='top')
    plt.ylim(0, 1)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.suptitle('Geminin, UniProt: O75496')
    plt.title('The Proportion of Analogs as a Function of Potency')
    plt.show()


def main():
    # # Initialize parser
    # parser = argparse.ArgumentParser()
    #
    # # Adding optional argument
    # parser.add_argument("-i", "--input", help="Input CSV from ChEMBL")
    # parser.add_argument("-B", "--Bemis", help="Group by Bemis-Murko Scaffold")
    # parser.add_argument("-F", "--Fingerprint", help="Cluster by SMILES Fingerprint")
    #
    # # Read arguments from command line
    # args = parser.parse_args()
    #
    # if args.Bemis:
        df = read_csv(csv_file='ligands/O75496/DOWNLOAD-s4J2I-ywcUMLMswllk75v0xJVZYFnVXldsps_qCh0_A=.csv')
        df = add_scaffold_and_group(df)
        df = rank_entries(df)
        plot(df, isLogarithmic=True)
    # if args.Fingerprint:
    #       # Read the CSV and calculate fingerprint representations
    #       df = read_csv(csv_file='Q9NUW8_lig_full.csv')
    #       df = fingerprint(df)
    #
    #       # Extract RDKit fingerprints from DataFrame
    #       fingerprints = df['Fingerprint'].tolist()
    #
    #       # Calculate Tanimoto distance matrix
    #       distance_matrix = calculateDistanceMatrix(fingerprints)


if __name__ == '__main__':
    main()
