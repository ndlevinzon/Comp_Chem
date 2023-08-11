import pandas as pd
from fuzzywuzzy import fuzz, process
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description="Completes truncated ZINC IDs from OUTDOCK using key from analogs.py")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input SMI file name.")
    parser.add_argument("-k", "--key", type=str, required=True, help=".SMI output from analogs.py")
    return parser.parse_args()

# Get command-line arguments
args = parse_arguments()

# Get the current working directory
current_directory = os.getcwd()

# Load the first CSV file into a DataFrame
df1 = pd.read_csv(os.path.join(current_directory, args.input))

# Load the second CSV file into a DataFrame
df2 = pd.read_csv(os.path.join(current_directory, args.key))

# Create an empty list to store the matched id_nums
matched_id_nums = []

# Iterate through each row in the first DataFrame
for index, row in df1.iterrows():
    id_num = row['id_num']

    if "_analog" in id_num:
        id_num_prefix = id_num.split("_analog")[0]
        id_num_suffix_length = len(id_num_prefix)
        best_match_suffix = \
        process.extractOne(id_num_prefix, df2['id_num'].str[-id_num_suffix_length:], scorer=fuzz.ratio)[0]
        best_match = df2.loc[df2['id_num'].str.endswith(best_match_suffix), 'id_num'].values[0]
        id_num_suffix = id_num.split("_analog")[1]
        best_match = best_match + "_analog" + id_num_suffix
    else:
        best_match = process.extractOne(id_num, df2['id_num'], scorer=fuzz.ratio)[0]

    matched_id_nums.append(best_match)
    print(f"Match for '{id_num}': '{best_match}'")

# Add the matched_id_nums list as a new column in df1
df1['id_num_new'] = matched_id_nums

# Save the updated DataFrame with the new column
df1.to_csv((os.path.join(current_directory, args.input, '_fuzzy')), index=False)
