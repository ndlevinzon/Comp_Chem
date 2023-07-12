# UCSF DOCK Tutorial
## Preparing your Protein

Items which are prefixed with 'AH' are relevant for docking HEIs to amidohydrolases and can safely be ignored for most metal-free proteins. 

### HEIs

High-energy intermediate databases in mol2 format can be found at [http://hei.docking.org].

The files available are:

* rmkh.ohln.tar.gz: KEGG, more than one ring, leaving group neutral
* rmkh.ohlp.tar.gz: KEGG, more than one ring, leaving group protonated
* rmpkh.ohln.tar.gz: KEGG, more than one ring rings, leaving group neutral
* rnkh.ohln.tar.gz: KEGG, no rings, leaving group neutral
* rnkh.ohlp.tar.gz: KEGG, no rings, leaving group protonated
* rnpkh.ohln.tar.gz: KEGG, no rings, leaving group neutral
* rokh.ohln.tar.gz: KEGG, one ring, leaving group neutral
* rokh.ohlp.tar.gz: KEGG, one ring, leaving group protonated
* ropkh.ohln.tar.gz: KEGG, one ring, leaving group neutral

##### Recommended data structure

Generate 4 subdirectories:
* 1_SOMENAME2SDF,
* 2_OMEGA,
* 3_AMSOL, and
* 4_MOL2DB
SMILES files must be of the format (somesmiles somename).
Copy the database preparation scripts:
```
a1*.py-c3*.py to 1_SOMENAME2SDF.
```

##### HEI generation
Broadly speaking, the generation procedure involves three steps:
1. Conversion of the input .mol2 or .sdf files to SMILES and then isomeric SMILES.
2. Conversion ("reaction") of the appropriate group(s) in each molecule to form the HEI.
3. Generation of multiple protonation states, 3D-structures and partial charges for each HEI, resulting in .db files that can be fed to UCSF DOCK.

The scripts for each step are prefixed with the letter 'a' (step 1), 'b' (step 2), and 'c' (first part of step 3), respectively. Within each letter, the scripts are enumerated consecutively.
* every script takes a list of SMILES as input and outputs a list of SMILES (except for the b7 and b8 scripts, which output .sdf files), prefixed with the sequential number of the script.
* the scripts a3-a4 have to be run in sequence.
* scripts b1 to b8 all take the output of a4 as input. Each of these scripts describes a different reaction and each reaction will only happen when    the appropriate reacting groups are encountered in a molecule.
* each b script will generate an LN (neutral leaving group) and an LP (protonated leaving group) file.
* it is a VERY GOOD idea to keep the LN and LP separate throughout the entire procedure, especially when running the c scripts. This will make things easier lateron.
* c3_sdf2mol2_mysql_names.py has to be run on four files: the LN and LP files coming out of c2_ionizer_min.py and the sdf files resulting from the b7 and b8 scripts.
* 
##### a Scripts
```
 a1_create_sdf_from_fold.py folder(unpacked from KEGG-website)
 a2_corina.py molfile.smi
 a2_corina.py can be started with either .smi, .ism or .sdf files.

 a3.1_sdf2ism_filter.py molfile.sdf 
 a3.2_size_cutoff_filter.py molfile.smi
 a4_rm_doubles.py molfile.smi
```
##### b Scripts
```
 b1_rxn_carbonyl.py molfile.smi
 b1_rxn_lactone.pk.py molfile.smi
 b2.1_rxn_aromatic_cleav.py molfile.smi
 b2.2_rxn_aromatic_cleav.py molfile.smi
 b3_rxn_amidines.py molfile.smi
 b4.1_rxn_amidine_aromatic.py molfile.smi
 b4.2_rxn_amidine_aromatic.py molfile.smi
 b4.3_rxn_amidine_aromatic.py molfile.smi
 b4.4_rxn_amidine_aromatic.py molfile.smi
 b5_rxn_imin.py molfile.smi
 b6.1_rxn_imin_aromatic.py molfile.smi
 b6.2_rxn_imin_aromatic.py molfile.smi
 b6.3_rxn_imin_aromatic.py molfile.smi
 b7.1_parts_split_1.py molfile.smi
 b7.2_parts_split_2.py file-identifier(e.g. _mol_32007)
 b7.3_parts_connect_1.py file-identifier(e.g. _mol_32007)
 b7.4_parts_connect_2.py file-identifier(e.g. _mol_32007)
 b7.5_ionizer.py file-identifier(e.g. _mol_32007)
 b8.1_thio_parts_split_1.py molfile.smi
 b8.2_thio_parts_split_2.py file-identifier(e.g. _mol_32007)
 b8.3_thio_parts_connect_1.py file-identifier(e.g. _mol_32007)
 b8.4_thio_parts_connect_2.py file-identifier(e.g. _mol_32007)
 b8.5_thio_ionizer.py file-identifier(e.g. _mol_32007)
 b9_remove_doubles.py start-file-pattern end-file-pattern
```
##### c Scripts
```
 c1_corina.py start-file-pattern end-file-pattern
 c2_ionizer_min.py start-file-pattern end-file-pattern
 c3_sdf2mol2_mysql_names.py sdf-file(from corina+ionizer) suffix(for Folders after ring)
 c3_sdf2mol2_mysql_names_remove.py filename_containing_mol2_filenames(zipped does not hurt)
```

* commandline:
```
mrm_3_limit.py MOL_RAID DB_VERSION MOLS_SUBDIR JOB_ID OMEGA_PATH AMSOL_PATH CHECK WRITE_BROKEN
```
* The individual arguments and the <tt>pwd</tt> will be connected to form the path to the mol2-files generated in step 2:
```
/raid[MOL_RAID]/people/kolb/DB[DB_VERSION]/[MOLS_SUBDIR]/MOLS/[obtained from pwd: penultimate dir]/[obtained from pwd: last dir]
```
* CHECK gives the frequency of the check whether a molecule has already been processed or not: '0' &rarr; no check; '1' &rarr; check at the beginning of every job; '2' &rarr; check before processing each molecule.
* in case the script stops after just one molecule, do the following:
** check that the file .labels[JOB_ID].txt exists.
** create a file .dbnums[JOB_ID].txt and write something like "101 0" to it. The first number will be the starting number for the enumeration of the .db files, while the second is the current number of molecules already in that .db file.
** delete everything but the header from the <tt>.db</tt> file.
** start <tt>mrm_3_limit.py</tt> again.

##### Inserting the newly generated molecules into a mysql database

This step is essential to preserve knowledge about the correspondence between the original database name of a molecule, its HEI form, protonation states and conformations and the final name given by mol2db (of the form A00000000 [one letter + eight digits]).
* required files:
** mysql_insert_db6.py
* this also requires you to generate a mysql database of the proper format, best done with mysql_create_table_db5.pk.py
* commandline:
 ```
mysql_insert_db6.py MOLS_SUBDIR MYSQL_DB MOL2_SUBDIR DB_SUBDIR MYSQL_TABLE DB_VERSION MOL_RAID TAG
```
* the individual arguments will be connected to form the path to the .mol2 and files as follows:
```
/raid[MOL_RAID]/people/kolb/DB[DB_VERSION]/[MOLS_SUBDIR]/MOLS/MOL2_SUBDIR
```
* the .db files are expected in ./DB_SUBDIR
* TAG is optional and is the name with which the molecule names start.

### Modifying the PDB file
 
* prepare rec.pdb  by removing all lines that do not commence with 'ATOM', all columns to the right of the z-coordinate and the TER statements.  
* treat all selenomethionines (MSE) as methionines (MET) by replacing the selenium atom (SE&curren;) with sulphur (&curren;SD). Be careful about the correct alignment!  
* atom enumeration does not matter, so don't bother to renumber after any of the following steps. Unique numbers are a good idea, presumably.  
* select the protonation states of HIS residues to be either &delta;- (rename residue to HID), &epsilon;- (rename residue to HIE) or doubly protonated (rename residue to HIP). HIS on the surface should be HIP. HIS residues coordinating the metal ions should have their protons pointing away from the ions. Base your decision on the immediate environment of the HIS residue: are there potential hydrogen bonds that can be formed?; are there charged residues close by?; would a certain protonation lead to clashes with other residues?; etc.
* AH: the carboxylated LYS of subtype I is CYK, but this is not tolerated by <tt>startdockblaster5</tt> , so store and delete the 3 surplus atoms and call the residue LYS.  
* AH: the more buried metal ion is ZB (charge 1.4), the other one ZA (charge 1.3). Atom names are right-aligned!

### Running startdockblaster5
   
* generate the file xtal-lig.pdb, which should only contain atoms of the molecular modeling force field. Prepare it in the same way as above: remove all columns to the right of the z-coordinate and the TER statements. Change HETATM to ATOM.
* generate the files <tt>.only_spheres</tt> and &ndash; in case you would like the matching spheres to be based on the heavy atoms in <tt>xtal-lig.pdb</tt> &ndash; <tt>.useligsph</tt> and write `on' to the latter. Be careful not to add blank lines at the end, this will not be understood by <tt>makespheres2.pl</tt> . In any case, the entry in <tt>.useligsph</tt>  will be ignored by <tt>makespheres1.pl</tt> .  
* on sgehead (or, as of [[dock67]], on any machine), run <tt>startdockblaster5</tt>  to set up the data structure and copy all relevant files. It is a good idea to use csh and to <tt>source .login</tt> beforehand.  
* if <tt>startdockblaster5</tt>  doesn't finish for any obvious reason and with no clear error message, or <tt>rec.crg</tt> has very funny hydrogen placements, make sure that you have no non-printing characters in <tt>rec.pdb</tt> or <tt>xtal-lig.pdb</tt>. Do that by running your file through <tt>pc2unix rec.pdb</tt>. Check that your file is clean by looking at it with <tt>od -c rec.pdb | less </tt>. The only character with a backslash should be \n &mdash; you should see no \t, \r, etc. If this doesn't solve the problem, your best bet is to re-prepare <tt>rec.pdb</tt> and <tt>xtal-lig.pdb</tt> from scratch &mdash; it is likely that there are some blanks or hidden characters that are causing the problems.  
* Take any WARNING messages emitted seriously, and continue only if you know why each one is there. Furthermore, verify that <tt>rec.crg</tt> still contains ''all'' atoms.  
* if you do not want to do anything special with the protein, like tarting some residues or modifying the spheres, go directly to chapter [[Running DOCK|3]].

### Removing and modifying files
   
* go to <tt>./grids</tt>  and remove the surplus files from this directory (some would cause error messages from the subsequent programs):<br><tt>rm -f PDBPARM chem.* rec+sph.phi solvmap_sev tart.txt OUT*</tt>
* modify <tt>rec.crg</tt>: 
** AH: CYK: put the three missing atoms, delete the surplus hydrogens specific for LYS and rename the carboxylated lysine residue CYK.  
** remove all TER statements that might have been added.  
** AH: set the atom names of the metal ions to ZA and ZB and the residue name to ZN.  
** take care of disulfide bonds. Remove the thiol hydrogens (if they have been added) and change the residue name from CYS to CYX.  
* look at the <tt>box</tt> and maybe move it, so that the ligands won't stick out. Modify the 'center' and 'coordinates' statement in the preamble.  
* all residues and atoms have to be listed in <tt>prot.table.ambcrg.ambH</tt> and <tt>vdw.parms.amb.mindock</tt>, respectively &rArr; do not tart any residues in this file!

### Running Chemgrid
   
* run <tt>chemgrid</tt> and check <tt>OUTPARM</tt> for the correct van der Waals parameters of all residues.  
* grep for <tt>0.000</tt>  in <tt>PDBPARM</tt>: if any atom has this value in the 3<sup>rd</sup> and 4<sup>th</sup> column, it has not been recognized by <tt>chemgrid</tt>  (because it is not listed in <tt>prot.table.ambcrg.ambH</tt>) and is thus ''ignored'' in the van der Waals-maps. There will be no other errors, the docking will finish showing some "bumping" ligands which have extremely favorable energies (&le; -200).
* Another sign of a problem with atomic radii are any 'WARNING's issued in OUTPARM
* if one has to run <tt>chemgrid</tt>  again, first remove <tt>PDBPARM OUTPARM OUTCHEM</tt> and <tt>chem.*</tt>.

### Tarting the protein
   
* cp <tt>rec.crg</tt> to <tt>rec+sph.crg</tt> and continue with the latter file.
* tarted residues can be found in <tt>$DOCK_BASE/scripts/grids</tt>, they are the files with the extension <tt>prot2</tt>.
* add the relevant resides to the bottom of your <tt>prot.table.ambcrg.ambH</tt> file, being very precise to match the current formatting
* generate the new <tt>amb.crg.oxt</tt> from the edited <tt>prot.table.ambcrg.ambH</tt> using:<br><tt>$mud/prot2crg.py < prot.table.ambcrg.ambH > amb.crg.oxt</tt>
* AH: select the appropriate version of <tt>amb.crg.oxt</tt> depending on the subtype. Files are called <tt>amb.crg.oxt.N</tt>, where <tt>N</tt> can be <tt>I, III</tt> or <tt>VI</tt>.  
* AH: edit the residues in the binding site (i.e., all residues complexing the metal ions in the binding site), so that their names conform to the names of the modified residues in <tt>amb.crg.oxt.N</tt>  
* optionally tart the residues that are in contact with a crystallographic ligand, if any.  
* AH: check that ZA and ZB, respectively (left-aligned in the atom column), have corresponding entries in <tt>amb.crg.oxt.N</tt>  and <tt>vdw.siz</tt>.

### Modifying the Delphi spheres
   
* load <tt>match1.sph.pdb</tt> (i.e., the DelPhi spheres).  
* delete the spheres that are too close to the solvent.  
* (AH:) add spheres so that there is one sphere ''between'' the metals, several spheres ''around'' the metals and some spheres close to polar residues.  
* a good number for DelPhi spheres is 120.  
* append the spheres to the end of <tt>rec.crg</tt> to make <tt>rec+sph.crg</tt> and put a TER statement after each sphere. Don't use tabs for whitespace, can cause problems with DelPhi!  

### Modifying the Matching spheres

* load <tt>match2.sph.pdb</tt> for sparse initial spheres or <tt>match3.sph.pdb</tt> denser spheres.
* If you selected <tt>.useligsph</tt> be careful not to move any spheres based on the ligand atoms.  
* (AH:) put at least one sphere between the metals and increase the sampling in the region around the metal ions by putting some spheres there.
* a good number for matching spheres is 50-60.  
* run <tt>pdbtosph matchN.sph.pdb mysph.sph</tt> to generate the files that will be read by [[DOCK]].   
* if color matching is desired, run <tt>colorspheres.pl sph/match2.sph</tt> in the parent directory of the docking run (i.e., <tt>..</tt> to <tt>sph</tt> ) to put some color on your spheres.  
* run <tt>cat $mud/header.sph match2.sph</tt> .

### Running Delphi

With the Makefile all someone needs to do is go up one directory and run the make command by typing:
```
 > cd ../
 > make grids/rec+sph.phi
```
And the Make file will do the work.

### Running newsolv.sev
   
* if you changed rec.crg or the box above, you need to run newsolv.sev   
* check that all atoms are present in <tt>rec.crg</tt> and run <tt>newsolv.sev</tt> .

## Preparing your Ligand

### Automatic way, starting from <tt>SMILES</tt>

This way, you will make use of John's automatic scripts for database preparation and actually upload new molecules to a special section of [http://zinc.docking.org/ ZINC].
   
* it is advisable to create a special subdirectory, since many new    files will be generated.  
* the file containing the [http://www.daylight.com/smiles/ SMILES] strings should contain a string    followed by an identifier on each line.  
* OPTIONAL: run <tt>convert.py --i=yourname.smi --o=yourname.ism</tt> . This will  convert your SMILES to ''isomeric'' SMILES.
* run <tt>dbgen.csh yourname.smi</tt>.  
** Note that the dbgen.csh does not work for more that 1000 molecules.
** Brake up your molecules into chunks of 1000 and run dbgen on each chunk.
** Clean up your directory afterwards. dbgen.csh generates a lot of files that you do not need if it ran correctly.  
* you should obtain a file <tt>somename.db.gz</tt> .

#### Optional
To increase the number of molecules that are written out for the database generation, copy the file $DOCK_BASE/data/omega.parm into the directory that dbgen.csh is going to be run in.
* At the end of the omega.parm file you will see a section called "Torsion Driving Parameters", here you will find three variables that can be changed.
** SetMaxConfs(600)   #set to higher numbers ie. 1000
** SetRMSThreshold(0.80)  #set to lower numbers ie. 0.50
** SetEnergyWindow(12.5)  #can be changed but this can often generate broken molecules
* WARNING this should only be done if generating conformations for a small set of compounds!!!

### Manual way

#### Isolating the ligand as <tt>.mol2</tt> file
   
* extract the ligand structure from the <tt>.pdb</tt>  file.  
* assign hydrogens.  
* assign all atom ([http://www.tripos.com/mol2/atom_types.html Sybyl/TAFF]) and bond types.  
* save it as <tt>ligandname.mol2</tt>  file.  

#### Running <tt>Omega</tt>
   
* run [http://www.eyesopen.com/products/applications/omega.html OMEGA], but don't ask me how to do that yet.

#### Running <tt>Amsol</tt>

* find more information about amsol [http://comp.chem.umn.edu/amsol/ on its homepage].  
* <tt>mkdir ./amsol2</tt>   
* Use file2file.py to get the right formal charge to feed to AMSOL. It is also important to change the name, otherwise the original <tt>.mol2</tt> file will be overwritten!
<tt>file2file.py -g ligandname.mol2 ./amsol2/someothername.mol2</tt>      
* edit <tt>./amsol2/someothername.mol2</tt> :     
** delete all lines prior to <tt>@<TRIPOS>MOLECULE</tt>   
** change line 2 (molecule name) to something of the format <tt>ABCD12345678</tt> (four capital letters followed by eight numbers).  
** line 3 should be <tt>n<sub>atoms</sub> n<sub>bonds</sub> 0 0 0</tt>
** the <tt>@<TRIPOS>MOLECULE</tt>  section must consist of exactly '''5''' lines (adjust by adding/deleting blanks).  
** remove all sections after the <tt>@<TRIPOS>BOND</tt>  section.
** delete the blank lines between the <tt>ATOM</tt>  and <tt>BOND</tt>     sections, if there are any.    
** run <tt>RunAMSOL3.csh WAIT</tt>   
* the output <tt>someothername.solv</tt>  file will contain the following:
{| style="text-align: center; border:1px solid #aaa; margin: 1em 1em 1em 0; background: #f9f9f9; border-collapse: collapse;" cellpadding="5" cellspacing="0" 
|+ '''AMSOL output'''
|-
! style="border:1px #aaa solid; padding: 0.2em;" | line #1
| style="border:1px #aaa solid; padding: 0.2em;" | molname || style="border:1px #aaa solid; padding: 0.2em;" | <math>n_{atoms}</math> || style="border:1px #aaa solid; padding: 0.2em;" | charge || style="border:1px #aaa solid; padding: 0.2em;" | pol_solv || style="border:1px #aaa solid; padding: 0.2em;" | ? || style="border:1px #aaa solid; padding: 0.2em;" | apol_solv || style="border:1px #aaa solid; padding: 0.2em;" | total_solv
|-
! style="border:1px #aaa solid; padding: 0.2em;" | other lines 
| style="border:1px #aaa solid; padding: 0.2em;" | charge || style="border:1px #aaa solid; padding: 0.2em;" | pol_solv || style="border:1px #aaa solid; padding: 0.2em;" | ? || style="border:1px #aaa solid; padding: 0.2em;" | apol_solv || style="border:1px #aaa solid; padding: 0.2em;" | total_solv
|-
| style="border:1px #aaa solid; padding: 0.2em;" | ''(per_atom)''
|}


*furthermore, there will be <tt>someothername.nmol2</tt>  file    which contains the correct partial charges.

#### Running <tt>Mol2DB</tt>
   
* edit <tt>someothername.nmol2</tt>  so that the <tt>@<TRIPOS>MOLECULE</tt> section consists of exactly '''6''' lines.  
* edit the <tt>inhier</tt>  file so that the 'mol2_file',    'db_file' and 'solvation_table' entries are correct.  
* run <tt>mol2db inhier</tt>   
* add the preamble at the top of the file.  
* <tt>gzip</tt>  the resulting file so that it can be used by <tt>DOCK</tt> .

## Understanding MakeDOCK

MakeDOCK automates the process of sphere and grid generation for a target.

### MakeDOCK Features

* If you have the MakeDOCK Makefile, you simply type 'make' to generate all the spheres and grids. 
* All required files and structures will be copied and created for you by MakeDOCK.  
* Automatically displays non-fatal WARNINGS that occurred during grid generation, which often lead to incorrect docking results.
* Full dependency resolution means you can change any input or parameter file and type 'make' to generate all files that depend on that change (i.e. tarting or editing spheres). 

### Input Files for MakeDOCK

Create a new directory containing the following input files to MakeDOCK: 
* rec.pdb: Prepared receptor file
* xtal-lig.pdb: Ligand specification file

These are required to specify the target for docking. For help preparing these files, see Receptor Preparation and Ligand Preparation.

### How to Use MakeDOCK

In the directory containing rec.pdb and xtal-lig.pdb, do the following to run MakeDOCK.
```
 setenv DOCK_BASE ~mysinger/xyz/dockenv
 source $DOCK_BASE/etc/login
 cp $DOCK_BASE/scripts/Makefile3 Makefile
 make
```

#### Making Changes

Any changes to any file used by MakeDOCK will automatically run all the commands necessary to generate all spheres and grids that depend on that change. For example, tart up a residue in grids/rec+sph.crg and type 'make' to regenerate delphi grids.

#### MakeDOCK Extras

* Once spheres and grids are generated, `make clean` will delete all input and log files, leaving you will the essential output files
* 'make distclean' will remove all files created by MakeDOCK

##### MakeDOCK as a Learning Tool

* Look inside the Makefile. The Makefile begins with variable definition and is then divided into target sections. The target sections have been specifically ordered to roughly follow how sphere and grid generation proceeds. Thus the Makefile can be used to learn exactly how spheres and grids are generated for DOCK.

* Variables are set using 'VARNAME = contents'. They are later referenced using $VARNAME. The general syntax for a target section is 'target: dependencies'. Changes to any of the dependencies cause that target to be made. Targets are made using the rules (lines) in that section. The command 'make' actually runs the 'all:' target. You can specify a different target using 'make target'. 

### Receptor MakeDOCK

To prepare a typical pdb structure for use as the receptor in MakeDOCK, do the following:

* Remove solvent (i.e. crystallographic waters)
* Remove other non-essential hetero-atoms (everything but perhaps a co-factor essential for ligand binding)
* Remove the ligand itself
* Delete all hydrogens
* Save as rec_original.pdb
* Run 'paranoia.csh' to clean up the rec_original.pdb file into the final rec.pdb file 

### Ligand MakeDOCK

To prepare the ligand specification, reload the original pdb structure and perform the following steps:

* Delete everything BUT the ligand atoms
* Save as xtal-lig.pdb

