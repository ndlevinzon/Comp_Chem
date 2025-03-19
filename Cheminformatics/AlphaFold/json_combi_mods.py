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

    def edit_residue_in_protein(self, protein_id, position, new_amino_acid, modifications):
        for protein in self.proteins:
            if protein_id in protein.protein_id:  # Match any listed ID
                protein.edit_residue(position, new_amino_acid)
                modifications.append((position, new_amino_acid))  # Track modifications
                return
        print(f"Error: Protein with ID '{protein_id}' not found.")

    def get_residue_codes(self, residue_numbers):
        """Prints the one-letter codes for given residue numbers in all stored proteins."""
        for protein in self.proteins:
            residues = protein.get_residues(residue_numbers)
            print(f"Protein {', '.join(protein.protein_id)}: Residues at {residue_numbers} -> {''.join(residues)}")

    def save_modified_json(self, original_filename="test.json", modifications=[]):
        """
        Saves the modified protein data to a new JSON file with a descriptive filename
        showing the amino acid changes and positions.

        Parameters:
        - original_filename (str): Name of the original file.
        - modifications (list of tuples): List of modifications as (position, new_amino_acid).
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
        new_filename = f"{base_name}_{modification_tag}{ext}"

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
    reader = ProteinSequenceReader("test.json")  # Load JSON file
    reader.load_data()

    # Print all loaded proteins
    for protein in reader.get_proteins():
        print(protein)

    aa_set = AminoAcidSet()  # Create the amino acid set
    comb_gen = AminoAcidCombinations(aa_set.get_acids())  # Pass the amino acid set

    comb_gen.print_combinations()  # Print all 4-letter combinations

    # Retrieve residues at specific positions
    residue_positions = [70, 71, 313, 314]
    reader.get_residue_codes(residue_positions)

    modifications = []

    # Example modifications
    reader.edit_residue_in_protein("A", 70, "A", modifications)
    reader.edit_residue_in_protein("A", 71, "Q", modifications)
    reader.edit_residue_in_protein("A", 313, "Q", modifications)
    reader.edit_residue_in_protein("A", 314, "Q", modifications)

    # Save modified JSON with descriptive filename
    reader.save_modified_json("./test.json", modifications)
