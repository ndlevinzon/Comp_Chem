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
import argparse


def read_csv(csv_file):
    """
    Read CSV file from ChEMBL Query and create a pandas DataFrame

    Parameters:
    csv_file (str): Path to the CSV file

    Returns:
    df (pandas.DataFrame): DataFrame containing the data from the CSV file
    """

    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file, lineterminator='\n', sep=';', dtype=str)

    # Create a new DataFrame with 'SMILES' and 'POTENCY' columns from the CSV data
    df = pd.DataFrame({'TARGET': data['Target ChEMBL ID'], 'SMILES': data['Smiles'],
                       'POTENCY': data['Standard Value']})

    # Drop rows with missing values in 'SMILES' and 'POTENCY' columns
    df.dropna(subset=['SMILES', 'POTENCY'], inplace=True)

    # Convert the 'POTENCY' column to numeric type
    df['POTENCY'] = pd.to_numeric(df['POTENCY'], errors='coerce')

    # Convert the potency units from nM to M by multiplying by 1e-9
    df['POTENCY_Converted'] = df['POTENCY'] * 1e-9

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

        # Check if scaffold is already in the group mapping
        if scaffold in group_mapping:
            # Return the existing group number
            return group_mapping[scaffold]
        else:
            # Add the scaffold to the group mapping and assign a new group number
            group_mapping[scaffold] = group_counter
            group_counter += 1
            # Return the newly assigned group number for the scaffold
            return group_mapping[scaffold]

    # Add 'GROUP' column to the DataFrame by applying the assign_group function to each scaffold
    df['GROUP'] = df['SCAFFOLD'].apply(assign_group)

    # Calculate the proportion of entries with a 10-fold reduction in POTENCY_Converted within each group
    df['Prop_10X'] = 0.0  # Initialize the 'Prop_10X' column with 0.0

    # Iterate through each group over 5 members
    group_counts = df['GROUP'].value_counts()
    groups_over_5 = group_counts[group_counts > 5].index.tolist()

    for group in groups_over_5:
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
    groups_over_25 = group_counts[group_counts > 25].index.tolist()
    print("Number of groups with more than 25 members:",
          len(groups_over_25))

    # Print the scaffold and group index for all groups over 25 members
    for group in groups_over_25:
        scaffold = df[df['GROUP'] == group]['SCAFFOLD'].iloc[0]
        print(f"Group {group}: Scaffold - {scaffold}")

    # Return the modified DataFrame with added 'SCAFFOLD', 'GROUP', and 'Prop_10X' columns
    return df


def rank_entries(df):
    """Rank the entries within each B-M family based on POTENCY values"""

    # Rank the entries within each group based on POTENCY values
    df['RANK'] = df.groupby('GROUP')['POTENCY'].rank(ascending=False, method='min')

    # Count the number of entries in each group
    group_counts = df['GROUP'].value_counts().reset_index()
    group_counts.columns = ['GROUP', 'COUNT']

    # Merge the group counts back into the DataFrame
    df = pd.merge(df, group_counts, on='GROUP')

    # Calculate the proportion of molecules with improvement within the same scaffold
    df['PROP_IMPROV'] = df['RANK'] / df['COUNT']

    # Return the modified DataFrame with added 'RANK', 'COUNT', and 'PROP_IMPROV' columns
    return df


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

    # Create a copy of the DataFrame to avoid the SettingWithCopyWarning
    df_copy = df.copy()

    # Calculate the heavy atom count for each SMILES and Bemis-Murcko scaffold
    df_copy['HeavyAtomCount'] = df_copy['SMILES'].apply(calculate_heavy_atom_count)
    df_copy['BemisMurckoScaffoldHeavyAtomCount'] = df_copy['SCAFFOLD'].apply(calculate_heavy_atom_count)

    # Iterate through unique targets
    for target in df_copy['TARGET'].unique():
        # Filter the DataFrame for the current target
        target_df = df_copy[df_copy['TARGET'] == target].copy()

        # Group by the 'GROUP' column and count the number of members in each group
        group_counts = target_df.groupby('GROUP').size()

        # Get the groups with fewer than 5 members
        small_groups = group_counts[group_counts < 5].index

        # Remove the small groups from the DataFrame
        df_copy = df_copy[~((df_copy['TARGET'] == target) & (df_copy['GROUP'].isin(small_groups)))]

        # Filter the DataFrame based on the heavy atom deviation
        df_cleaned = df_copy[abs(df_copy['HeavyAtomCount'] - df_copy['BemisMurckoScaffoldHeavyAtomCount']) <= 2].copy()

        # Remove the extra columns
        df_cleaned.drop(['HeavyAtomCount', 'BemisMurckoScaffoldHeavyAtomCount'], axis=1, inplace=True)

        # Return the filtered DataFrame
        return df_cleaned


def create_scatter(df, subtitle, group_number=None):
    """Plot the families sharing the same B-M Scaffold"""
    plt.figure()

    y_data = df['PROP_IMPROV']
    y_label = 'Prop. Molecules with Improvement within the Same Scaffold'

    x_data = df['POTENCY_Converted']
    x_data = -np.log10(x_data)  # Convert x_data to -log(x_data)
    x_label = '-log(Potency) (M)'
    plt.xscale("log")  # Scale x-axis to logarithmic

    # Select the groups with more than 25 members
    group_counts = df['GROUP'].value_counts()

    if group_number is not None:
        # Plot a specific group if group_number is specified
        group_data = df[df['GROUP'] == group_number]

        # Print the scaffold for the group
        scaffold = group_data['SCAFFOLD'].iloc[0]
        print(f"Group {group_number}: Scaffold - {scaffold}")

        lowest_potency_index = group_data['POTENCY'].idxmin()
        lowest_potency_smiles = group_data.loc[lowest_potency_index, 'SMILES']
        print(f"Group: {group_number}, Lowest POTENCY SMILES: {lowest_potency_smiles}")

        highest_potency_index = group_data['POTENCY'].idxmax()
        highest_potency_smiles = group_data.loc[highest_potency_index, 'SMILES']
        print(f"Group: {group_number}, Highest POTENCY SMILES: {highest_potency_smiles}")

        # Each group will have a unique color
        color = plt.cm.Set1(group_number % plt.cm.Set1.N)
        plt.scatter(x_data[group_data.index], y_data[group_data.index], color=color, alpha=0.6, s=10)

    else:
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

        # Combine data points from the groups over 5 members
        groups_over_5 = group_counts[(group_counts > 5) & (group_counts.index != '')].index.tolist()
        for group in groups_over_5:
            group_data = df[df['GROUP'] == group]
            plt.scatter(x_data[group_data.index], y_data[group_data.index], alpha=0.6, s=10)

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
    valid_groups = group_counts[group_counts > 5].index.tolist()

    # Set up the histogram plot
    plt.figure()
    plt.hist(df[df['GROUP'].isin(valid_groups)]['Prop_10X'], bins=20, range=(0, 1), edgecolor='black')
    plt.xlabel('Success Rate \n (10-fold Increase in Potency from Parent with Same B-M Scaffold)')
    plt.ylabel('Observations')
    plt.suptitle(subtitle)
    plt.title('Distribution of Success Rates for Analogs')
    plt.show()


def main(csv_file, database_file):
    subtitle = 'Target Name'

    # Connect to the SQLite database
    engine = create_engine(f'sqlite:///{database_file}')

    df = read_csv(csv_file=csv_file)
    df = add_scaffold_and_group(df)
    df = rank_entries(df)
    df_cleaned = clean_df(df)

    # Graph data
    create_scatter(df=df_cleaned, subtitle=subtitle, group_number=375)
    create_histogram(df=df_cleaned, subtitle=subtitle)

    # Check if the database file exists
    if not os.path.isfile(database_file):
        # Store the DataFrame in the database
        df.to_sql('histogram_data', con=engine, index=False)
        print("Database created and data inserted.")
        
    else:
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
    # Usage: python script.py --csv input.csv --database data.db
    parser = argparse.ArgumentParser(description='Process ChEMBL data and create histograms')
    parser.add_argument('--csv', required=True, help='Path to the input CSV file')
    parser.add_argument('--database', required=True, help='Path to the database file')
    args = parser.parse_args()

    main(csv_file=args.csv, database_file=args.database)

