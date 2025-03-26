import json
import os
from itertools import product


class Protein:
    """Represents a protein with an ID and sequence."""

    def __init__(self, protein_id, sequence):
        self.protein_id = protein_id
        self.sequence = list(sequence)  # Store sequence as a mutable list

    def get_sequence(self):
        """Returns the protein sequence as a string."""
        return "".join(self.sequence)

    def get_residues(self, residue_numbers):
        """
        Returns the one-letter codes of specified residues.

        Parameters:
        residue_numbers (list of int): List of residue positions (1-based index).

        Returns:
        list of str: Residue codes corresponding to the positions.
        """
        residues = []
        for num in residue_numbers:
            if 1 <= num <= len(self.sequence):  # Ensure valid index
                residues.append(self.sequence[num - 1])  # Convert to 0-based indexing
            else:
                residues.append("X")  # Placeholder for invalid positions
        return residues

    def edit_residue(self, position, new_amino_acid):
        """
        Edits a specific residue in the sequence.

        Parameters:
        position (int): 1-based index of the residue to be changed.
        new_amino_acid (str): The new amino acid (must be one-letter code).
        """
        if not isinstance(new_amino_acid, str) or len(new_amino_acid) != 1:
            print(f"Error: '{new_amino_acid}' is not a valid one-letter amino acid.")
            return

        if 1 <= position <= len(self.sequence):  # Ensure position is valid
            self.sequence[position - 1] = new_amino_acid  # Modify sequence in place
            print(f"Residue at position {position} changed to {new_amino_acid}.")
        else:
            print(f"Error: Position {position} is out of range.")

    def to_dict(self):
        """Returns a dictionary representation of the protein for JSON formatting."""
        return {
            "protein": {
                "id": self.protein_id,
                "sequence": self.get_sequence()
            }
        }

    def __str__(self):
        """String representation of the protein."""
        return f"Protein ID: {', '.join(self.protein_id)}\nSequence: {self.get_sequence()}"


class ProteinSequenceReader:
    """Reads a protein sequence from a JSON file and allows in-place editing."""

    def __init__(self, file_path):
        self.file_path = file_path
        self.proteins = []  # List to store Protein objects

    def load_data(self):
        """Loads the protein sequence data from the JSON file."""
        try:
            with open(self.file_path, "r") as file:
                data = json.load(file)

            sequences = data.get("sequences", [])
            for entry in sequences:
                protein_data = entry.get("protein", {})
                protein_id = protein_data.get("id", [])
                sequence = protein_data.get("sequence", "")

                if protein_id and sequence:
                    protein = Protein(protein_id, sequence)
                    self.proteins.append(protein)

        except FileNotFoundError:
            print(f"Error: The file '{self.file_path}' was not found.")
        except json.JSONDecodeError:
            print("Error: Failed to decode JSON. Please check the file format.")

    def get_proteins(self):
        """Returns the list of proteins."""
        return self.proteins

    def edit_residue_in_protein(self, protein_id, positions, new_amino_acids):
        """
        Edits multiple residues in a given protein using a 4-mer amino acid sequence.

        Parameters:
        - protein_id (str): The ID of the protein to modify.
        - positions (list of int): List of residue positions to be modified.
        - new_amino_acids (str): Four-letter string representing new amino acids.

        Returns:
        - modifications (list of tuples): Tracked modifications (for filename generation).
        """
        modifications = []
        for protein in self.proteins:
            if protein_id in protein.protein_id:  # Match any listed ID
                for pos, new_aa in zip(positions, new_amino_acids):
                    protein.edit_residue(pos, new_aa)
                    modifications.append((pos, new_aa))  # Track modifications
                return modifications
        print(f"Error: Protein with ID '{protein_id}' not found.")
        return []

    def get_residue_codes(self, residue_numbers):
        """Prints the one-letter codes for given residue numbers in all stored proteins."""
        for protein in self.proteins:
            residues = protein.get_residues(residue_numbers)
            print(f"Protein {', '.join(protein.protein_id)}: Residues at {residue_numbers} -> {''.join(residues)}")

    def save_modified_json(self, original_filename="test.json", modifications=[], modification_id=""):
        """
        Saves the modified protein data to a new JSON file with a descriptive filename.

        Parameters:
        - original_filename (str): Name of the original file.
        - modifications (list of tuples): List of modifications as (position, new_amino_acid).
        - modification_id (str): A unique identifier for the modification batch (e.g., "AAPR").
        """
        modified_data = {
            "name": "Modified_Protein_Data",
            "sequences": [protein.to_dict() for protein in self.proteins],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1
        }

        # Generate a filename that includes modifications
        modification_tag = "_".join([f"{pos}{aa}" for pos, aa in modifications])
        base_name, ext = os.path.splitext(original_filename)
        new_filename = f"{base_name}_{modification_id}{ext}"

        # Save the modified JSON to a file
        with open(new_filename, "w") as file:
            json.dump(modified_data, file, indent=2)

        print(f"Modified JSON saved to: {new_filename}")


class AminoAcidSet:
    """Stores the 20 standard amino acid one-letter codes."""

    def __init__(self):
        self.amino_acids = set("ACDEFGHIKLMNPQRSTVWY")  # Set for fast lookup

    def get_acids(self):
        """Returns the set of amino acids."""
        return self.amino_acids


class AminoAcidCombinations:
    """Generates and prints all 4-letter combinations of amino acids with repetition."""

    def __init__(self, amino_acid_set):
        self.amino_acid_list = sorted(list(amino_acid_set))  # Sort for consistent ordering

    def generate_combinations(self, length=4):
        """Generates all possible combinations of a given length."""
        return product(self.amino_acid_list, repeat=length)

    def print_combinations(self, length=4):
        """Prints each generated combination and counts them."""
        count = 0
        for combination in self.generate_combinations(length):
            print("".join(combination))  # Convert tuple to string and print
            count += 1
        print(f"\nTotal number of combinations found: {count}")  # Print final count


# Init
if __name__ == "__main__":
    reader = ProteinSequenceReader("test.json")
    reader.load_data()

    # Generate 4-letter amino acid combinations
    aa_set = AminoAcidSet()
    comb_gen = AminoAcidCombinations(aa_set.get_acids())

    # Define positions to modify
    positions_to_modify = [70, 71, 313, 314]

    # Iterate over each combination and create a modified JSON file
    for combination in comb_gen.generate_combinations(length=4):
        modification_sequence = "".join(combination)  # Convert tuple to string
        modification_id = f"{positions_to_modify[0]}{modification_sequence[0]}_" \
                          f"{positions_to_modify[1]}{modification_sequence[1]}_" \
                          f"{positions_to_modify[2]}{modification_sequence[2]}_" \
                          f"{positions_to_modify[3]}{modification_sequence[3]}"

        # Apply mutations to a copy of the reader
        modifications = reader.edit_residue_in_protein("A", positions_to_modify, modification_sequence)

        # Save the modified sequence as a new JSON file
        if modifications:
            reader.save_modified_json("test.json", modifications, modification_id)
