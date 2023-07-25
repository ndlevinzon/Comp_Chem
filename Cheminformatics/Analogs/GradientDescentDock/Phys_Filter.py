import os
import pandas as pd
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate pKa/b, logP, and Lipinski Rules for .smi")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input SMILES file name.")
    return parser.parse_args()

def main():
    # Get command-line arguments
    args = parse_arguments()
    
    # Get the current working directory
    current_directory = os.getcwd()

    # Define the input file path
    smi_filename = args.input  # Replace this with your actual file name or path
    input_file = os.path.join(current_directory, smi_filename)
    print(f'Input File: {input_file}')

    # Define the output of the cxcalc calculation
    cxcalc_output = os.path.join(current_directory, 'cxcalc_temp')
    print(f'CXCALC Output File: {cxcalc_output}')

    # Define the output of the cxcalc calculation
    RO5_output = os.path.join(current_directory, 'RO5_temp')
    print(f'RO5 Output File: {RO5_output}')

    # Run cxcalc operation and read the 'RO5_output' file once
    df_cxcalc = pKa(input_file, cxcalc_output)
    df_combined = RO5(df_cxcalc, input_file, RO5_output)

    # Output the combined DataFrame to a CSV file
    df_combined.to_csv(f'{smi_filename}_filtered.csv', index=False)

    # Delete _temp files
    os.remove(cxcalc_output)
    os.remove(RO5_output)


def pKa(input_file, cxcalc_output):
    # Read the input file into a pandas DataFrame
    df_smi = pd.read_csv(input_file, sep=' ', header=None, names=['smiles_codes', 'zinc_codes'])
    print(df_smi)

    # Run cxcalc operation
    cxcalc_cmd = f'/mnt/nfs/home/nlevinzon/freechem/freechem-19.15.r4/bin/cxcalc -i smiles {input_file} logP pKa > {cxcalc_output}'
    subprocess.run(cxcalc_cmd, shell=True)

    # Read the second output file into a pandas DataFrame
    df_cx = pd.read_csv(cxcalc_output, sep='\t', header=None,
                        names=['smiles', 'logP', 'apKa1', 'apKa2', 'bpKa1', 'bpKa2', 'atoms'])

    # Drop the first row
    df_cx.drop(0, inplace=True)

    # Drop the redundant 'smiles' column
    df_cx.drop(columns='smiles', inplace=True)

    print(df_cx)

    # Reset the index of both DataFrames before concatenation
    df_smi.reset_index(drop=True, inplace=True)
    df_cx.reset_index(drop=True, inplace=True)

    # Merge the DataFrames based on 'smiles_codes' to match
    # 'logP', 'apKa1', 'apKa2', 'bpKa1', 'bpKa2' to the correct smiles_codes
    df_cxcalc = pd.concat([df_smi, df_cx], axis=1, ignore_index=True)

    # Rename the columns
    df_cxcalc.columns = ['smiles_code',
                           'zinc_code',
                           'logP',
                           'apKa1',
                           'apKa2',
                           'bpKa1',
                           'bpKa2',
                           'atoms']

    # Print the cxcalc dataframe
    print(df_cxcalc)
    return df_cxcalc


def RO5(df_cxcalc, input_file, RO5_output):
    # Run the RO5 evaluation
    RO5_cmd = \
        f'/mnt/nfs/home/nlevinzon/freechem/freechem-19.15.r4/bin/evaluate -e "mass() <= 500 && logP() <= 5 && donorCount() <= 5 && acceptorCount() <= 10 && rotatableBondCount() <= 10 && PSA() <= 140" {input_file} > {RO5_output}'
    subprocess.run(RO5_cmd, shell=True)

    # Read the RO5 output file line by line into a list
    with open(RO5_output, 'r') as ro5_file:
        ro5_data = [int(line.strip()) for line in ro5_file if line.strip().isdigit()]

    # Create a DataFrame with the RO5 data
    df_RO5 = pd.DataFrame(ro5_data, columns=['RO5'])

    # If the DataFrame is empty, fill it with the first encountered integer value
    if df_RO5.empty and ro5_data:
        df_RO5['RO5'] = ro5_data[0]

    # Merge the DataFrames to the correct smiles_codes and concatenate the DataFrames along the columns axis
    df_combined = pd.concat([df_cxcalc, df_RO5], axis=1, ignore_index=True)

    # Rename the columns
    df_combined.columns = ['smiles_code',
                           'zinc_code',
                           'logP',
                           'apKa1',
                           'apKa2',
                           'bpKa1',
                           'bpKa2',
                           'atoms',
                           'RO5?']

    # Print the combined DataFrame
    print(df_combined)
    return df_combined


if __name__ == '__main__':
    main()
