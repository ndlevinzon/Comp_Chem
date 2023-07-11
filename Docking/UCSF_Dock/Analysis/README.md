# Combining the results of all subdirectories
   
* in the subdirectory that contains all the individual directories for each chunk of the library, run <tt>$mud/combine.py</tt>. Then generate a file containing the top 500 molecules using '<tt>$mud/topdock.py -o top500.pdb</tt>', which you can read into ViewDOCK in chimera as a DOCK 4, 5, or 6 style file.
* to create an <tt>.eel1</tt> file containing the top 500 molecules just run <tt>$mud/topdock.py -e</tt>. If one wants to create an <tt>.eel1</tt> file for a different subset of the molecules, first create the list of molecule names plus their energies (on one line) and then feed it to <tt>getxpdb.pl name_energy.list < FF.test.eel1 > subset_name.eel1</tt>.

# Getting individual atom contributions with scoreopt_so

## First you need an <tt>.eel1</tt> file to be scored
### For the xtal-lig.mol2 in its crystallographic pose
New way that outputs your.eel1 starting from your.pdb directly
*run '<tt>$mud/to_eel1.csh your.pdb</tt>'. 

If that fails, use the old way to convert an input <tt>[http://www.tripos.com/index.php?family=modules,SimplePage,,,&page=sup_mol2&s=0 .mol2]</tt> file into an <tt>.eel1</tt> file
* run <tt>amsol</tt>  as described [[Preparing_the_ligand#Running amsol|here]] to calculate atomic solvation energies.
* run '<tt>file2file.py -s path/to/amsol.solv path/to/amsol.nmol2 ligand.eel1</tt>'.

###For molecules that have already been docked

* run '<tt>$mud/topdock.py -e -o top500.eel1' to generate an .eel1 containing the top 500 docked molecules.
* or unzip the dock output '<tt>gunzip -c test.eel1.gz > test.eel1</tt>'
* or to create an <tt>.eel1</tt> file for a different subset of the molecules, first create the list of molecule names plus their energies (on one line) and then feed it to '<tt>getxpdb.pl name_energy.list < FF.test.eel1 > subset_name.eel1</tt>'.

## Overall molecular score compiled from all scoreopt_so options

For default grids
* run <tt>'$mud/doscoreopt.csh your.eel1 ../path/to/grids'</tt>
Or for custom grids, used below to run SEV-based desolvation grids
* run <tt>'$mud/doscoreopt.csh your.eel1 ../path/to/grids rec+sph.phi chem solvmap_sev'</tt> 
The summary for the whole molecule is output to your.eel1.scores in combine.scores format 

## Atomic contributions to the coulombic energy

In your.eel1.delphi from the wrapper 
* in every ATOM line, columns 9, 10 and 11 are the partial charge, the electrostatic field and the energy in kT (i.e., 9 &times; 10) of the atom, respectively.  
* the DelPhi electrostatic score is the sum over the entries in column 11 times 0.5924 (conversion from kT to kcal/mol) and can be compared to the elect column in OUTDOCK.
Or to generate this data yourself
* start <tt>scoreopt_so</tt>  and choose option '2' in the first menu.  
* enter the name of the DelPhi potential file, presumably <tt>grids/rec+sph.phi</tt>.  
* enter the name of the ligand file, i.e., <tt>ligand.eel1</tt> or <tt>top500.eel1</tt>.  
* enter the name of the output file, e.g. <tt>ligand.delphi</tt> .

## Atomic contributions to the van der Waals energy

In your.eel1.vdw from the wrapper  
* be adequately [http://www.merriam-webster.com/dictionary/scared scared].  
* the van der Waals interaction energy is calculated as  <math>{vdW}_{(r)}=\frac{A}{r^{12}}-\frac{B}{r^6}=a-b</math>. In every ATOM line, columns 9, 10 and 11 are <math>a</math>, <math>b</math> and <math>a-b</math>,    respectively.
* DO NOT use the interaction energy, as we only use the vdw component now. Instead, use the vdwsum to compare with the vdW column in OUTDOCK.
Or to generate this data yourself
* start <tt>scoreopt_so</tt>  and choose option '3' in the first menu.  
* enter the prefix name of grids for ff scoring as a full path,    i.e., <tt>grids/chem</tt> .  
* enter the name of the van der Waals parameter file, presumably    <tt>grids/vdw.parms.amb.mindock</tt> .  
* answer the question about interpolation with 'yes'.  
* enter a sufficiently large number as maximal van der Waals    energy, e.g. 10000.  
* enter the name of the ligand file, i.e., <tt>ligand.eel1</tt> .  
* enter the name of the output file, e.g. <tt>ligand.vdw</tt> .  


## Atomic contributions to the desolvation
   
In your.eel1.solv from the wrapper  
* in every ATOM line, columns 9, 10, and 11 are the total atomic solvation energy (polar + apolar), percentage desolvation, and atomic desolvation energy (i.e. - 9 &times; 10) of the atom, respectively.
* the total desolvation is the sum over the entries in column 11 and can be compared to the sum of the polsol and apolsol columns in OUTDOCK.
Or to generate this data yourself
* start <tt>scoreopt_so</tt>  and choose option '4' in the first menu.  
* enter the name of the grid for partial desolvation, presumably <tt>grids/solvmap</tt> or <tt>grids/solvmap_sev</tt>.  
* enter the name of the ligand file, i.e., <tt>ligand.eel1</tt> .  
* enter the name of the output file, e.g. <tt>ligand.solv</tt> .

## Obtaining the net charge of a docked molecule
   
* take the output <tt>.eel1</tt> file and run <tt>molcharge_pdb.pl < output.eel1</tt>. This will output the sequential number of the molecule, the [http://zinc.docking.org/ ZINC] identifier, the total charge and the number of atoms for every molecule in the file.
