import re
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def extract_group(input_string):
    result = re.search(r'(?:ZINC|INC|NC|C)([A-Za-z0-9]+)(?:_analog|\.\d+|$)', input_string)
    if result:
        return result.group(1)
    else:
        return None


def find_group_parents(group_df):
    parent_indices = {}

    for group_id, group_data in group_df.groupby('GROUP'):
        non_analog_entries = group_data[~group_data['id_num'].str.contains('_analog')]

        if len(non_analog_entries) > 1:
            parent_index = non_analog_entries['Total_Converted'].idxmin()
            parent_indices[group_id] = parent_index
        elif len(non_analog_entries) == 1:
            parent_indices[group_id] = non_analog_entries.index[0]

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
    df['Total_Converted'] = (np.exp((df['Total'] *scalar) / (300 * 0.001987)))

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
            other_entries_lower = group_data[group_data['Total_Converted'] < total_converted]

            if len(group_data) > 2:  # Modified condition to handle highest and lowest entries
                if total_converted == min_total_converted:
                    prop_improv = 0
                elif total_converted == max_total_converted:
                    prop_improv = 1
                else:
                    prop_improv = abs((len(other_entries_lower)) / (len(group_data)-2))
            else:
                prop_improv = np.nan

            df.at[index, 'PROP_IMPROV'] = format(prop_improv, ".3f")  # Format to three decimal places

    print(df.head(10))

    df.to_csv('C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_grouped.csv', index=False)
    # JK: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/jk/jk_grouped.csv'
    # Fangyu: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/fangyu/fangyu_grouped.csv'
    # Alpha2a: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_grouped.csv'

    return df

def create_scatter(df, subtitle):
    plt.figure()

    y_data = pd.to_numeric(df['PROP_IMPROV'], errors='coerce')  # Convert to numeric
    y_label = 'Proportion of Improvement'

    x_data = df['Total_Converted']
    x_label = 'Calculated p(Ki) From UCSF DOCK 3.8'

    # Generate a random color for each group
    unique_colors = ['#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for _ in range(df['GROUP'].nunique())]

    # Assign markers to denote parent ('*') and non-parent ('o') of a group
    for i, group_number in enumerate(df['GROUP'].unique()):
        group_data = df[df['GROUP'] == group_number]

        color = unique_colors[i]
        label = f'Group {group_number}'

        plt.scatter(-np.log(x_data[group_data.index]), y_data[group_data.index], color=color, marker='o', alpha=0.6, s=30,
                    label=label)

        parent_index = group_data['Parent_Index'].values[0] if group_data['Parent_Index'].count() > 0 else None
        marker = '*' if parent_index in group_data.index else 'o'  # Use '*' marker if the entry is a parent, 'o' otherwise

        if marker == '*' and parent_index is not None:
            parent_data = group_data[group_data.index == parent_index]
            plt.scatter(-np.log(x_data[parent_data.index]), y_data[parent_data.index], color=color, marker='*', alpha=0.8, s=100,
                        label='Parent')

    plt.xscale("log")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title("The Proportion of Improved Analogs as a Function of p(Ki)")
    plt.suptitle(subtitle)
    plt.ylim([0, 1])
    # plt.gca().invert_xaxis()  # Invert the x-axis

    plt.show()


def main():
    # Source CSV (extract_all.py -> formatting in Excel ->)
    csv_file = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_extract_all.csv'
    # JK: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/jk/jk_extracted.csv'
    # Fangyu: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/ampc/fangyu/fangyu_extract_all_fuzzy.csv'
    # Alpha2a: 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/alpha2a/alpha2a_extract_all.csv'

    subtitle = 'Fink, E. A. et al. (2022) DOI:10.1126/science.abn7065'
    # JK: 'Lyu, J. et al. (2019) DOI: 10.1038/s41586-019-0917-9'
    # Fangyu: 'Liu, F. (2023) (unpublished)'
    # Alpha2a: 'Fink, E. A. et al. (2022) DOI:10.1126/science.abn7065'

    # Read the new data from CSV
    df=read_csv(csv_file=csv_file)

    # Graph data
    create_scatter(df=df, subtitle=subtitle)


if __name__ == '__main__':
    main()
