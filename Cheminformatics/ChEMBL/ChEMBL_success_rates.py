# ND Levinzon, UCSF 2023
#  __     ________ _   _ _______
#  \ \   / /  ____| \ | |__   __|/\
#   \ \_/ /| |__  |  \| |  | |  /  \
#    \   / |  __| | . ` |  | | / /\ \
#     | |  | |____| |\  |  | |/ ____ \
#     |_|  |______|_| \_|  |_/_/    \_\


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sqlalchemy import create_engine


def read_csv(csv_file):
    """
    Read CSV file from ChEMBL Query and create a pandas DataFrame

    Parameters:
    csv_file (str): Path to the CSV file

    Returns:
    df (pandas.DataFrame): DataFrame containing the data from the CSV file
    """
    
    data = pd.read_csv(csv_file, lineterminator='\n', sep=';', dtype=str)             # Read the CSV file into a pandas DataFrame

    # Create a new DataFrame with 'SMILES' and 'POTENCY' columns from the CSV data
    df = pd.DataFrame({'TARGET': data['Target ChEMBL ID'], 'SMILES': data['Smiles'],
                       'POTENCY': data['Standard Value']})

    df.dropna(subset=['SMILES', 'POTENCY'], inplace=True)                             # Drop rows with missing values in 'SMILES' and 'POTENCY' columns
    df['POTENCY'] = pd.to_numeric(df['POTENCY'], errors='coerce')                     # Convert the 'POTENCY' column to numeric type
    df['POTENCY_Converted'] = df['POTENCY'] * 1e-9                                    # Convert the potency units from nM to M by multiplying by 1e-9

    return df


def add_scaffold_and_group(df):
    """Add 'SCAFFOLD', 'GROUP', and 'Prop_10X' columns to the DataFrame"""

    # Add 'SCAFFOLD' column to the DataFrame
    df['SCAFFOLD'] = df['SMILES'].apply(
        lambda x: Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(str(x)))))

    # Initialize dictionary and counter for grouping
    group_mapping = {}
    group_counter = 1
    def assign_group(scaffold):
        nonlocal group_counter

        if scaffold in group_mapping:                                          # Check if scaffold is already in the group mapping
            return group_mapping[scaffold]                                     # Return the existing group number
        else:
            group_mapping[scaffold] = group_counter                            # Add the scaffold to the group mapping and assign a new group number
            group_counter += 1
            return group_mapping[scaffold]                                     # Return the newly assigned group number for the scaffold

    # Add 'GROUP' column to the DataFrame by applying the assign_group function to each scaffold
    df['GROUP'] = df['SCAFFOLD'].apply(assign_group)

    # Calculate the proportion of entries with a 10-fold reduction in POTENCY_Converted within each group
    df['Prop_10X'] = 0.0  # Initialize the 'Prop_10X' column with 0.0

    # Iterate through each group over 25 members
    group_counts = df['GROUP'].value_counts()
    groups_over_25 = group_counts[group_counts > 25].index.tolist()

    for group in groups_over_25:
        group_indices = df[df['GROUP'] == group].index

        for idx in group_indices:
            potency = df.loc[idx, 'POTENCY_Converted']
            other_potencies = df.loc[group_indices, 'POTENCY_Converted']
            reduction_count = sum(other_potencies <= (potency/10))
            prop_10x = reduction_count / (len(group_indices)) if len(group_indices) > 1 else 0.0

            df.loc[idx, 'Prop_10X'] = prop_10x

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
    print("Number of groups with more than 25 members:",
          len(groups_over_25))

    # Print the scaffold and group index for all groups over 25 members
    for group in groups_over_25:
        scaffold = df[df['GROUP'] == group]['SCAFFOLD'].iloc[0]
        print(f"Group {group}: Scaffold - {scaffold}")
        
    return df  # Return the modified DataFrame with added 'SCAFFOLD', 'GROUP', and 'Prop_10X' columns


def rank_entries(df):
    """Rank the entries within each B-M family based on POTENCY values"""

    df['RANK'] = df.groupby('GROUP')['POTENCY'].rank(ascending=True, method='min')     # Rank the entries within each group based on POTENCY values

    # Count the number of entries in each group
    group_counts = df['GROUP'].value_counts().reset_index()
    group_counts.columns = ['GROUP', 'COUNT']

    df = pd.merge(df, group_counts, on='GROUP')                                        # Merge the group counts back into the DataFrame
    df['PROP_IMPROV'] = df['RANK'] / df['COUNT']                                       # Calculate the proportion of molecules with improvement within the same scaffold

    return df                                                                          # Return the modified DataFrame with added 'RANK', 'COUNT', and 'PROP_IMPROV' columns


def clean_df(df):
    """Filters DataFrame for small groups and deviations between SMILES and B-M Scaffold"""

      # Define a function to calculate the heavy atom count
    def calculate_heavy_atom_count(smiles):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            Chem.Kekulize(molecule)  # Kekulize the molecule to ensure proper atom types
            heavy_atom_count = rdMolDescriptors.CalcNumHeavyAtoms(molecule)
            return heavy_atom_count
        else:
            return None
    
    df_copy = df.copy()  # Create a copy of the DataFrame to avoid the SettingWithCopyWarning

    # Calculate the heavy atom count for each SMILES and Bemis-Murcko scaffold
    df_copy['HeavyAtomCount'] = df_copy['SMILES'].apply(calculate_heavy_atom_count)
    df_copy['BemisMurckoScaffoldHeavyAtomCount'] = df_copy['SCAFFOLD'].apply(calculate_heavy_atom_count)

    # Iterate through unique targets
    for target in df['TARGET'].unique():
        target_df = df[df['TARGET'] == target]                                                                # Filter the DataFrame for the current target
        group_counts = target_df.groupby('GROUP').size()                                                      # Group by the 'GROUP' column and count the number of members in each group                        
        small_groups = group_counts[group_counts < 25].index                                                  # Get the groups with fewer than 25 members
        df_copy = df_copy[~((df_copy['TARGET'] == target) & (df_copy['GROUP'].isin(small_groups)))]           # Remove the small groups from the DataFrame

    df_filtered = df_copy[abs(df_copy['HeavyAtomCount'] - df_copy['BemisMurckoScaffoldHeavyAtomCount']) <= 2] # Filter the DataFrame based on the heavy atom deviation
    df_filtered.drop(['HeavyAtomCount', 'BemisMurckoScaffoldHeavyAtomCount'], axis=1, inplace=True)           # Remove the extra columns

    return df_filtered                                                                                        # Return the filtered DataFrame


def graph(df, subtitle):
    """Plot the families sharing the same B-M Scaffold"""
    plt.figure()

    y_data = df['PROP_IMPROV']
    y_label = 'Prop. Molecules with Improvement within the Same Scaffold'

    x_label = '-log(Potency) (M)'
    x_data = df['POTENCY_Converted']
    x_data = -np.log10(x_data)  # Convert x_data to -log(x)
    plt.xscale("log")  # Scale x-axis to logarithmic if isLogarithmic is True

    # Select the groups with more than 25 members
    group_counts = df['GROUP'].value_counts()

    # Print from the top five groups
    top_groups = group_counts[group_counts.index != ''].nlargest(5).index.tolist()
    print("Top 5 Most Populous Groups:", top_groups)
    for group in top_groups:
        group_data = df[df['GROUP'] == group]

        # Print the scaffold for the group
        scaffold = group_data['SCAFFOLD'].iloc[0]
        print(f"Group {group}: Scaffold - {scaffold}")

        lowest_potency_index = group_data['POTENCY'].idxmin()
        lowest_potency_smiles = group_data.loc[lowest_potency_index, 'SMILES']
        print(f"Group: {group}, Lowest POTENCY SMILES: {lowest_potency_smiles}")

        highest_potency_index = group_data['POTENCY'].idxmax()
        highest_potency_smiles = group_data.loc[highest_potency_index, 'SMILES']
        print(f"Group: {group}, Highest POTENCY SMILES: {highest_potency_smiles}")

    # Combine data points from the groups over 25 members
    combined_x = []
    combined_y = []
    groups_over_25 = group_counts[(group_counts > 25) & (group_counts.index != '')].index.tolist()
    for group in groups_over_25:
        group_data = df[df['GROUP'] == group]
        combined_x.extend(x_data[group_data.index])
        combined_y.extend(y_data[group_data.index])

        # Each group will have a unique color
        color = plt.cm.Set1(group % plt.cm.Set1.N)
        plt.scatter(x_data[group_data.index], y_data[group_data.index], color=color, alpha=0.6, s=10)

    plt.ylim(0, 1)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.suptitle(subtitle)
    plt.title('The Proportion of Analogs as a Function of Potency')
    plt.show()


def create_histogram(df, subtitle):
    """Create a histogram of 'Prop_10X' values for groups with more than 25 members"""

    # Select groups with more than 25 members
    group_counts = df['GROUP'].value_counts()
    valid_groups = group_counts[group_counts > 25].index.tolist()

    # Set up the histogram plot
    plt.figure()
    plt.hist(df[df['GROUP'].isin(valid_groups)]['Prop_10X'], bins=20, range=(0, 1), edgecolor='black')
    plt.xlabel('Success Rate \n (10-fold Increase in Potency from Parent with Same B-M Scaffold)')
    plt.ylabel('Frequency')
    plt.suptitle(subtitle)
    plt.title('Distribution of Success Rates for Analogs')
    plt.show()


def main():
    # Source CSV from ChEMBL
    csv_file = 'ligands/Q9UNA4/DOWNLOAD-rTvoiuJRKQ_g5uD2vsXC9Wni1wpZE_I730Ak9kgxFis=.csv'
    subtitle = 'Glucagon-like peptide 1 receptor, UniProt: P43220'

    # Connect to the SQLite database
    database_file = 'ligands/histogram_data_filtered.db'
    engine = create_engine(f'sqlite:///{database_file}')

    # Check if the database file exists
    if not os.path.isfile(database_file):
        df = read_csv(csv_file=csv_file)
        df = add_scaffold_and_group(df)
        df = rank_entries(df)
        df = clean_df(df)

        # Graph data
        graph(df=df, subtitle=subtitle)
        create_histogram(df=df, subtitle=subtitle)

        # Store the DataFrame in the database
        df.to_sql('histogram_data', con=engine, index=False)
        print("Database created and data inserted.")
    else:
        # Read the new data from CSV
        df = read_csv(csv_file=csv_file)
        df = add_scaffold_and_group(df)
        df = rank_entries(df)
        df = clean_df(df)

        # Graph data
        graph(df=df, subtitle=subtitle)
        create_histogram(df=df, subtitle=subtitle)

        # Retrieve the existing data from the database
        query = "SELECT TARGET FROM histogram_data"
        existing_targets = pd.read_sql_query(query, con=engine)['TARGET'].tolist()

        # Find new entries not present in the database
        new_entries = df[~df['TARGET'].isin(existing_targets)]

        if not new_entries.empty:
            # Store the new entries in the database
            new_entries.to_sql('histogram_data', con=engine, index=False, if_exists='append')
            print("New data inserted into the database.")
        else:
            print("No new data to insert.")

    # Retrieve the updated data from the database
    query = "SELECT * FROM histogram_data"
    print(pd.read_sql_query(query, con=engine))


if __name__ == '__main__':
    main()
