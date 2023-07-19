import re
import pandas as pd
import subprocess

# Define the input file path
input_file = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/phenol_benchmark-analogs-i1-o96.smi'

# Define the first output file path
output_file1 = 'C:/Users/ndlev/PycharmProjects/shoichet/analogs/phenol/smiles_test.txt'

# Define the second output file path
output_file2 = '/mnt/nfs/home/nlevinzon/benchmarks/output.txt'

# Define the regular expression pattern to capture both smile_codes and zinc_codes
pattern = r'^(\S+)\s+(\S+)$'

# Open the input file for reading
with open(input_file, 'r') as file:
    # Read the contents of the file
    contents = file.read()

    # Use regex to extract the SMILES_CODE and ZINC_CODE columns
    matches = re.findall(pattern, contents, re.MULTILINE)

# Open the first output file for writing (creating if it doesn't exist)
with open(output_file1, 'w') as file:
    # Write the extracted SMILES_CODE and ZINC_CODE values to the first output file
    for match in matches:
        file.write(f"{match[0]}\t{match[1]}\n")

# Run the command in the terminal
command = f"./cxcalc -i {output_file1} logP pKa > {output_file2}"
subprocess.run(command, shell=True)

# Read the first output file into a pandas DataFrame
df1 = pd.read_csv(output_file1, sep='\t', header=None, names=['smiles_codes', 'zinc_codes'])

# Read the second output file into a pandas DataFrame
df2 = pd.read_csv(output_file2, sep='\t', header=None, names=['smiles', 'logP', 'apKa1', 'apKa2', 'bpKa1', 'bpKa2', 'atoms'])

# Combine the two DataFrames
df_combined = pd.concat([df1, df2[['logP', 'apKa1', 'apKa2', 'bpKa1', 'bpKa2', 'atoms']]], axis=1)

# Print the combined DataFrame
print(df_combined)
