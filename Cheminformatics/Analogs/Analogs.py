# Nathan Levinzon, Jonathan Borowsky, and Olivier Mailhot | UCSF 2023
# Version 1.0
import argparse
import time
import os
import pandas as pd
from rdkit import RDLogger
from collections import defaultdict
from Analog_Methods import *


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate analogs from a .smi file.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input SMI file name.")
    return parser.parse_args()


def main():
    # Get command-line arguments
    args = parse_arguments()

    # Get the current time before running the code
    start_time = time.time()
    print(f"Starting Time: {start_time}")

    # Silence RDKit error messages
    RDLogger.DisableLog('rdApp.error')

    # Define the analogue methods and their corresponding names
    analog_methods = [
        [trim_extremities, "trim"],
        [BO_stepup, "increase-bond-order"],
        [BO_stepdown, "decrease-bond-order"],
        [ring_breaker, "ring-opening"],
        [ring_maker, "ring-closure"],
        [nitrogen_walks, "n-walk"],
        [oxygen_walks, "o-walk"],
        [sulfur_walks, "s-walk"],
        [heterocycle_walks, "heterocycle-walk"],
        [methyl_scanning, "ch3-scan"],
        [amine_scanning, "nh2-scan"],
        [hydroxyl_scanning, "oh-scan"],
        [thiol_scanning, "sh-scan"],
        [fluorine_scanning, "f-scan"],
        [chlorine_scanning, "cl-scan"],
        [bromine_scanning, "br-scan"],
        [iodine_scanning, "i-scan"],
        [bioisosters_scanning, "bioisoster-scan"],
    ]

    # Specify the input and output file names
    path = os.getcwd()
    smi_input_filename = args.input

    # Use os.path.splitext() to split the filename and its extension
    filename_without_extension, file_extension = os.path.splitext(smi_input_filename)

    # Read the input file and store the smiles and zinc IDs in a DataFrame
    smiles_zinc_input = pd.read_csv(os.path.join(path, smi_input_filename),
                                    sep=' ', header=None, names=['Smiles', 'ZincID'])

    # Create a dictionary to store the analogs for each input SMILES code and their zinc IDs
    all_analogs_dict = defaultdict(set)  # Use sets for analogs

    # Loop through each input molecule
    for _, row in smiles_zinc_input.iterrows():
        smiles = row['Smiles']
        zinc_id = row['ZincID']

        # Initialize the counter for each parent compound (reset to zero for each parent)
        analogs_count = 0

        # Loop through each analogue method for the current input molecule
        for method in analog_methods:
            # Generate analogues using the specified method
            for analog in method[0](smiles):
                try:
                    # Enumerate stereoisomers for the current analog
                    stereoisomers = enumerate_stereoisomers(analog)
                    if stereoisomers is None:
                        analog_smiles = Chem.MolToSmiles(analog, isomericSmiles=True)
                        if not analog_smiles:
                            raise ValueError(f"Invalid SMILES: {analog_smiles}")
                        # Add analog to the set and store zinc ID
                        all_analogs_dict[smiles].add((analog_smiles, zinc_id))
                        analogs_count += 1
                        # Print analog SMILES and create fakezinc inside the loop
                        fakezinc = str(zinc_id) + "_analog" + str(analogs_count).zfill(4)
                        print(analog_smiles, fakezinc)
                    else:
                        for isomer in stereoisomers:
                            analog_smiles = Chem.MolToSmiles(isomer, isomericSmiles=True)
                            if not analog_smiles:
                                raise ValueError(f"Invalid SMILES: {analog_smiles}")
                            # Add analog to the set and store zinc ID
                            all_analogs_dict[smiles].add((analog_smiles, zinc_id))
                            analogs_count += 1
                            # Print analog SMILES and create fakezinc inside the loop
                            fakezinc = str(zinc_id) + "_analog" + str(analogs_count).zfill(4)
                            print(analog_smiles, fakezinc)
                except ValueError as e:
                    print(str(e))  # Skip the entry and print the error message
                    continue  # Skip to the next iteration

    # Write the generated smiles and fake zincs to a file
    total_analogs_count = 0
    output_file = os.path.join(path,
                               f"{filename_without_extension}_analogs-i{len(smiles_zinc_input)}-o{(sum(len(analogs) for _, analogs in all_analogs_dict.items()))}.smi")
    print(f"\nWriting to {output_file}")
    with open(output_file, "w") as f2:
        for _, row in smiles_zinc_input.iterrows():
            smiles = row['Smiles']
            zinc_id = row['ZincID']
            f2.write(f"{smiles} {zinc_id}\n")
            parent_analogs = all_analogs_dict[smiles]
            for i, (analog_smiles, parent_zinc_id) in enumerate(parent_analogs, start=1):
                fakezinc = str(parent_zinc_id) + "_analog" + str(i).zfill(4)
                f2.write(f"{analog_smiles} {fakezinc}\n")
                total_analogs_count += 1

    # Calculate the number of analogs generated
    print(f"Molecules: {total_analogs_count} molecules")

    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Runtime: {elapsed_time} seconds")

    # Calculate molecules generated per second
    benchmark = float(sum(len(analogs) for _, analogs in all_analogs_dict.items())) / elapsed_time
    print(f"Benchmark: {benchmark} molecules/second")

if __name__ == '__main__':
    main()
