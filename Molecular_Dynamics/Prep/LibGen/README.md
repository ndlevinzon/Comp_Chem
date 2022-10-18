# LibGen.pl
Script that takes a directory filled with ligand .PDB files, runs Gaussian16 geometry/charge optimization, and generates AMBER-compatible ligand parameters for use in Molecular Dynamics simulations

# Usage
Open LibGen.pl in a text editor and change the following lines in order to match your UNIX environment
```
# Load Modules
require "$ENV{\"LMOD_PKG\"}/init/perl";
module("load openbabel/2.4.1");                                 # OpenBabel
module("load gaussian16/SSE4.C01");                             # GAUSSIAN
module("load gcc/8.5.0 intel-oneapi-mpi/2021.4.0 amber/20.20"); # AMBER
