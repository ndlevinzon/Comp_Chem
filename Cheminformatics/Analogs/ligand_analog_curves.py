import re
import random
import pandas as pd
import numpy as np
from fuzzywuzzy import fuzz
import matplotlib.pyplot as plt


def extract_group(input_string):
    result = re.search(r'(?:ZINC|INC|NC|C)([A-Za-z0-9]+)(?:_analog|\.\d+|$)', input_string)
    if result:
        return result.group(1)
    else:
        return None


def find_best_parent_index(group_data):
    best_similarity = 0
    best_parent_index = None

    for index, row in group_data.iterrows():
        for parent_index, parent_row in group_data.iterrows():
            if index != parent_index:
                similarity = fuzz.ratio(row['id_num'], parent_row['id_num'])
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_parent_index = parent_index

    return best_parent_index


def find_group_parents(group_df):
    parent_indices = {}

    for group_id, group_data in group_df.groupby('GROUP'):
        non_analog_entries = group_data[~group_data['id_num'].str.contains('_analog')]

        if len(non_analog_entries) > 1:
            parent_index = non_analog_entries['Total_Converted'].idxmin()
            parent_indices[group_id] = parent_index
        elif len(non_analog_entries) == 1:
            parent_indices[group_id] = non_analog_entries.index[0]
        else:
            # Call the find_best_parent_index function for groups without a clear parent
            best_parent_index = find_best_parent_index(group_data)
            if best_parent_index is not None:
                parent_indices[group_id] = best_parent_index

    return parent_indices


def read_csv(csv_file, scalar=1):
    # Scalars:
    # JK: 9.41E53
    # Fangyu:
    # Alpha2a: 8.14E18
    """
    Read CSV file from formatted OUTDOCK and create a pandas DataFrame

    Parameters:
    csv_file (str): Path to the CSV file

    Returns:
    df (pandas.DataFrame): DataFrame containing the data from the CSV file
    """
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file, sep=',')

    # Convert the 'POTENCY' column to numeric type
    df['Total'] = pd.to_numeric(df['Total'], errors='coerce')

    # Convert the potency units from nM to M by multiplying by 1e-9 and a scalar
    df['Total_Converted'] = (np.exp((df['Total'] * scalar) / (300 * 0.001987)))

    # Apply the function to the 'id_num' column and create a new column 'group'
    df['GROUP'] = df['id_num'].apply(extract_group)

    # Assign group numbers sequentially
    group_mapping = {group: group_count for group, group_count in
                     zip(df['GROUP'].unique(), range(1, len(df['GROUP'].unique()) + 1))}

    # Map the group numbers to the 'group' column
    df['GROUP'] = df['GROUP'].map(group_mapping)

    # Convert the 'GROUP' column to integers
    df['GROUP'] = df['GROUP'].astype('Int64')

    # Group the DataFrame by 'GROUP' and apply the find_group_parents function
    group_parents = df.groupby('GROUP').apply(find_group_parents)

    # Create a dictionary that maps group IDs to parent indices
    group_to_parent = {group: parent_index for group, parent_dict in group_parents.items() for group, parent_index in
                       parent_dict.items()}

    # Map group IDs to their parent indices and assign to all entries
    df['Parent_Index'] = df['GROUP'].map(group_to_parent)

    # Calculate PROP_IMPROV values for each entry in the DataFrame
    df['PROP_IMPROV'] = None  # Initialize PROP_IMPROV column

    for group_number in df['GROUP'].unique():
        group_data = df[df['GROUP'] == group_number]
        min_total_converted = group_data['Total_Converted'].min()
        max_total_converted = group_data['Total_Converted'].max()

        for index in group_data.index:
            total_converted = group_data.loc[index, 'Total_Converted']

            if total_converted == min_total_converted:
                prop_improv = 0
            elif total_converted == max_total_converted:
                prop_improv = 1
            else:
                other_entries_lower = group_data[group_data['Total_Converted'] < total_converted]
                prop_improv = abs((len(other_entries_lower)) / (len(group_data)-1))

            df.at[index, 'PROP_IMPROV'] = float(prop_improv)  # Format to three decimal places

    print(df.head(10))

    df.to_csv('C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_grouped.csv', index=False)
    # JK: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/jk/jk_grouped.csv'
    # Fangyu: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/fangyu/fangyu_grouped.csv'
    # Alpha2a: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_grouped.csv'

    return df


def group_meets_criteria(group_data):
    # Define your criteria here
    # Check if the group has a PROP_IMRPOV with a value of 1,
    # a PROP_IMPROV with a value of 0, and a single parent
    # Return True if the group meets the criteria, else False
    # Example:
    return (1 in group_data['PROP_IMPROV'].values) and \
           (0 in group_data['PROP_IMPROV'].values) and \
           (len(group_data['Parent_Index'].dropna().unique()) == 1)


def group_has_valid_parent(group_data, df):
    # Define your criteria for a valid parent here
    # Example: Check if the parent index is not null and exists in the DataFrame
    parent_index = group_data['Parent_Index'].iloc[0]
    return parent_index is not np.nan and parent_index in df.index


def create_scatter(df, subtitle):
    plt.figure()

    y_data = pd.to_numeric(df['PROP_IMPROV'], errors='coerce')  # Convert to numeric
    y_label = 'Proportion of Improvement'

    x_data = df['Total_Converted']
    x_label = 'Calculated p(Ki) From UCSF DOCK 3.8'

    # # Scatter plot for all data points
    # plt.scatter(-np.log(x_data), y_data, color='gray', marker='o', alpha=0.4, s=20, label='All Data Points')

    # Select the top five most populous groups
    top_groups = df['GROUP'].value_counts().nlargest(5).index

    # Generate a random color for each group
    unique_colors = ['#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for _ in range(len(top_groups))]

    # Loop through each group
    for i, group_number in enumerate(top_groups):
        group_data = df[df['GROUP'] == group_number]

        # Check if the group meets the criteria and has a valid parent
        if group_meets_criteria(group_data) and group_has_valid_parent(group_data, df):
            color = unique_colors[i]
            label = f'Group {group_number}'

            parent_indices = group_data['Parent_Index'].values

            if parent_indices.size > 0:
                parent_min_potency = float('inf')
                min_potency_index = None

                # Find the instance of the parent with the lowest POTENCY_CONVERTED value
                for parent_index in parent_indices:
                    parent_data = group_data[group_data.index == parent_index]

                    if not parent_data.empty:
                        parent_potency = parent_data['Total_Converted'].values[0]

                        if parent_potency < parent_min_potency:
                            parent_min_potency = parent_potency
                            min_potency_index = parent_index

                if min_potency_index is not None:
                    # Plot the parent with the lowest POTENCY_CONVERTED value as a star marker
                    plt.scatter(-np.log(x_data[min_potency_index]), y_data[min_potency_index], color=color, marker='*',
                                alpha=1, s=100, label='Parent')

            # Plot other group members with 'o' marker
            non_parent_indices = group_data.index.difference(parent_indices)
            plt.scatter(-np.log(x_data[non_parent_indices]), y_data[non_parent_indices], color=color, marker='o', alpha=0.6,
                        s=10, label=label)

    plt.xscale('symlog')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title("The Proportion of Improved Analogs as a Function of p(Ki) \n (Top 5 Most Populated Groups)")
    plt.suptitle(subtitle)
    plt.ylim([0, 1])
    plt.show()


def main():
    # Source CSV (extract_all.py -> formatting in Excel ->)
    csv_file = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/fangyu/fangyu_extract_all_fuzzy.csv'
    # JK: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/jk/jk_extracted.csv'
    # Fangyu: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/fangyu/fangyu_extract_all_fuzzy.csv'
    # Alpha2a: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_extract_all.csv'

    subtitle = 'Liu, F. (2023) (unpublished)'
    # JK: 'Lyu, J. et al. (2019) DOI: 10.1038/s41586-019-0917-9'
    # Fangyu: 'Liu, F. (2023) (unpublished)'
    # Alpha2a: 'Fink, E. A. et al. (2022) DOI:10.1126/science.abn7065'

    # Read the new data from CSV
    df=read_csv(csv_file=csv_file)

    # Graph data
    create_scatter(df=df, subtitle=subtitle)


if __name__ == '__main__':
    main()
