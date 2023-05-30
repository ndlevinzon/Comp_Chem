import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def readCSV(csv_file):
    """Read CSV file from ChEMBL Query, create a pandas DataFrame, and add 'SCAFFOLD' and 'GROUP' columns"""
    data = pd.read_csv(csv_file, lineterminator='\n', sep=';', dtype=str)

    # Create a DataFrame with SMILES and POTENCY_VALUE columns
    df = pd.DataFrame({'SMILES': data['Smiles'], 'POTENCY': data['Standard Value']})

    # Drop rows with missing SMILES values
    df.dropna(subset=['SMILES'], inplace=True)

    # Convert POTENCY column to numeric type
    df['POTENCY'] = pd.to_numeric(df['POTENCY'], errors='coerce')

    # Add 'SCAFFOLD' column
    df['SCAFFOLD'] = df['SMILES'].apply(lambda x:
                                        Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(str(x)))))

    # Add 'GROUP' column
    group_mapping = {}
    group_counter = 1

    def assign_group(scaffold):
        nonlocal group_counter
        if scaffold in group_mapping:
            return group_mapping[scaffold]
        else:
            group_mapping[scaffold] = group_counter
            group_counter += 1
            return group_mapping[scaffold]

    df['GROUP'] = df['SCAFFOLD'].apply(assign_group)

    # Add 'PROP_IMPROV' column
    def rank_within_group(group, potency):
        group_purities = df[df['GROUP'] == group]['POTENCY'].values
        ranks = np.argsort(group_purities)
        return ranks.tolist().index(np.where(group_purities == potency)[0][0])

    df['PROP_IMPROV'] = df.apply(lambda row: rank_within_group(row['GROUP'], row['POTENCY']), axis=1)

    # Print the total number of groups found
    print("Total number of groups found:", group_counter - 1)

    # Print the number of groups with more than 1 member
    print("Number of groups with more than 1 member:", len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 1]))

    # Print the number of groups with more than 5 members
    print("Number of groups with more than 5 members:", len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 5]))

    # Print the number of groups with more than 10 members
    print("Number of groups with more than 10 members:",
          len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 10]))

    # Print the number of groups with more than 25 members
    print("Number of groups with more than 25 members:",
          len(df['GROUP'].value_counts()[df['GROUP'].value_counts() > 25]))

########################################################################################################################

    # Store clusters with more than 25 members
    large_clusters = df['GROUP'].value_counts()[df['GROUP'].value_counts() > 25].index.tolist()

    # Create a graph for all clusters
    plt.figure()

    # Plot each cluster as a separate curve
    for cluster in large_clusters:
        cluster_data = df[df['GROUP'] == cluster]
        x = cluster_data['POTENCY']
        y = cluster_data['PROP_IMPROV']

        plt.plot(x, y, label=f'Cluster {cluster}')  # Each cluster has a different label

    plt.xlabel('Potency')
    plt.ylabel('PROP_IMPROV')
    plt.title('Group Over 25 Members')
    plt.legend()  # Show the legend with cluster labels
    plt.show()

    print(df)

    return df


def main():
    # Read the CSV and calculate Bemis-Murko scaffolds and Proportion values
    readCSV(csv_file='Q9NUW8_lig_full.csv')


if __name__ == '__main__':
    main()
