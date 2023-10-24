# FTIR Project

setup_pdbs2.py: a script that I wrote to modify pdb input files to allow leap to read them.
converting H2O to D2O
pdbNucleicAcidParse.py: formats source pdb into expected .list outputs

1. Generate nucleic acid PDB using Avogadro
2. Convert using pdb4amber.py utility
3. Run pdbNucleicAcidParse.py to generate ".list" files
4. Build using Amber tLeap, convert H2O to D2O
5. Minimize and equilibrate
6. Run MD
