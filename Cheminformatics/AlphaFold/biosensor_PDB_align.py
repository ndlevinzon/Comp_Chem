from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
import numpy as np
import os

def align_ligand_and_mg_to_alphafold(af_pdb_path, domain_pdb_path, ligand_resname, mg_resname="MG", domain_chain_id="A", output_ligand_path="aligned_ligand_and_mg.pdb"):
    """
    Aligns a ligand and bound magnesium ion from a domain structure to an AlphaFold structure by aligning the domain.

    Parameters:
    - af_pdb_path (str): Path to AlphaFold full protein PDB
    - domain_pdb_path (str): Path to domain PDB containing ligand and Mg
    - ligand_resname (str): Residue name of the ligand (e.g., "ATP")
    - mg_resname (str): Residue name of the magnesium ion (default "MG")
    - domain_chain_id (str): Chain ID of the domain in both structures (default "A")
    - output_ligand_path (str): Output path for the aligned ligand+Mg PDB
    """
    parser = PDBParser(QUIET=True)

    # Load both structures
    af_structure = parser.get_structure("AF", af_pdb_path)
    domain_structure = parser.get_structure("DOMAIN", domain_pdb_path)

    def get_ca_atoms(structure, chain_id):
        return [res["CA"] for res in structure[0][chain_id] if "CA" in res and res.id[0] == " "]

    # Extract CA atoms
    try:
        af_atoms = get_ca_atoms(af_structure, domain_chain_id)
        domain_atoms = get_ca_atoms(domain_structure, domain_chain_id)
    except KeyError as e:
        raise ValueError(f"Chain {domain_chain_id} not found in one of the structures.") from e

    if len(af_atoms) != len(domain_atoms):
        min_len = min(len(af_atoms), len(domain_atoms))
        af_atoms = af_atoms[:min_len]
        domain_atoms = domain_atoms[:min_len]

    # Superimpose domain onto AlphaFold
    super_imposer = Superimposer()
    super_imposer.set_atoms(domain_atoms, af_atoms)
    super_imposer.apply(domain_structure.get_atoms())

    # Check for ligand and Mg presence
    ligand_residues = [res for res in domain_structure.get_residues() if res.get_resname() == ligand_resname]
    mg_residues = [res for res in domain_structure.get_residues() if res.get_resname() == mg_resname]

    if not ligand_residues:
        raise ValueError(f"Ligand '{ligand_resname}' not found in the domain structure.")
    if not mg_residues:
        raise ValueError(f"Magnesium ion '{mg_resname}' not found in the domain structure.")

    class LigandAndMgSelect(Select):
        def accept_residue(self, residue):
            return residue.get_resname() in {ligand_resname, mg_resname}

    # Write transformed ligand and magnesium
    io = PDBIO()
    io.set_structure(domain_structure)
    io.save(output_ligand_path, LigandAndMgSelect())

    print(f"Ligand and Mg aligned and written to: {output_ligand_path}")

# Example usage:
# align_ligand_and_mg_to_alphafold("alphafold.pdb", "domain_ligand.pdb", "ATP", "MG", "A", "aligned_ligand_and_mg.pdb")
