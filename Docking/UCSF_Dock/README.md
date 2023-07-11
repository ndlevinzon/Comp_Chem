# UCSF DOCK Tutorial
## Preparing your Protein

Items which are prefixed with 'AH' are relevant for docking [[HEI]]s to amidohydrolases and can safely be ignored for most metal-free proteins. 

### Modifying the PDB file
 
* prepare <tt>rec.pdb</tt>  by removing all lines that do not commence with 'ATOM', all columns to the right of the z-coordinate and the TER statements.  
* treat all selenomethionines (MSE) as methionines (MET) by replacing the selenium atom (SE&curren;) with sulphur (&curren;SD). Be careful about the correct alignment!  
* atom enumeration does not matter, so don't bother to renumber after any of the following steps. Unique numbers are a good idea, presumably.  
* select the protonation states of HIS residues to be either &delta;- (rename residue to HID), &epsilon;- (rename residue to HIE) or doubly protonated (rename residue to HIP). HIS on the surface should be HIP. HIS residues coordinating the metal ions should have their protons pointing away from the ions. Base your decision on the immediate environment of the HIS residue: are there potential hydrogen bonds that can be formed?; are there charged residues close by?; would a certain protonation lead to clashes with other residues?; etc.
* AH: the carboxylated LYS of subtype I is CYK, but this is not tolerated by <tt>startdockblaster5</tt> , so store and delete the 3 surplus atoms and call the residue LYS.  
* AH: the more buried metal ion is ZB (charge 1.4), the other one ZA (charge 1.3). Atom names are right-aligned!

### Running startdockblaster5
   
*generate the file <tt>xtal-lig.pdb</tt> , which should only contain atoms of the molecular modeling force field. Prepare it in the same way as above: remove all columns to the right of the z-coordinate and the TER statements. Change HETATM to ATOM.
*generate the files <tt>.only_spheres</tt> and &ndash; in case you would like the matching spheres to be based on the heavy atoms in <tt>xtal-lig.pdb</tt> &ndash; <tt>.useligsph</tt> and write `on' to the latter. Be careful not to add blank lines at the end, this will not be understood by <tt>makespheres2.pl</tt> . In any case, the entry in <tt>.useligsph</tt>  will be ignored by <tt>makespheres1.pl</tt> .  
*on sgehead (or, as of [[dock67]], on any machine), run <tt>startdockblaster5</tt>  to set up the data structure and copy all relevant files. It is a good idea to use csh and to <tt>source .login</tt> beforehand.  
*if <tt>startdockblaster5</tt>  doesn't finish for any obvious reason and with no clear error message, or <tt>rec.crg</tt> has very funny hydrogen placements, make sure that you have no non-printing characters in <tt>rec.pdb</tt> or <tt>xtal-lig.pdb</tt>. Do that by running your file through <tt>pc2unix rec.pdb</tt>. Check that your file is clean by looking at it with <tt>od -c rec.pdb | less </tt>. The only character with a backslash should be \n &mdash; you should see no \t, \r, etc. If this doesn't solve the problem, your best bet is to re-prepare <tt>rec.pdb</tt> and <tt>xtal-lig.pdb</tt> from scratch &mdash; it is likely that there are some blanks or hidden characters that are causing the problems.  
*Take any WARNING messages emitted seriously, and continue only if you know why each one is there. Furthermore, verify that <tt>rec.crg</tt> still contains ''all'' atoms.  
*if you do not want to do anything special with the protein, like tarting some residues or modifying the spheres, go directly to chapter [[Running DOCK|3]].

### Removing and modifying files
   
* go to <tt>./grids</tt>  and remove the surplus files from this directory (some would cause error messages from the subsequent programs):<br><tt>rm -f PDBPARM chem.* rec+sph.phi solvmap_sev tart.txt OUT*</tt>
* modify <tt>rec.crg</tt>: 
* *AH: CYK: put the three missing atoms, delete the surplus hydrogens specific for LYS and rename the carboxylated lysine residue CYK.  
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
   
*load <tt>match1.sph.pdb</tt> (i.e., the DelPhi spheres).  
*delete the spheres that are too close to the solvent.  
*(AH:) add spheres so that there is one sphere ''between'' the metals, several spheres ''around'' the metals and some spheres close to polar residues.  
*a good number for DelPhi spheres is 120.  
*append the spheres to the end of <tt>rec.crg</tt> to make <tt>rec+sph.crg</tt> and put a TER statement after each sphere. Don't use tabs for whitespace, can cause problems with DelPhi!  

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
 > cd ../
 > make grids/rec+sph.phi

And the Make file will do the work.

### Running newsolv.sev
   
* if you changed rec.crg or the box above, you need to run newsolv.sev   
* check that all atoms are present in <tt>rec.crg</tt> and run <tt>newsolv.sev</tt> .

## Preparing your Ligand

### Automatic way, starting from SMILES

This way, you will make use of John's automatic scripts for database
preparation and actually upload new molecules to a special section of
[http://zinc.docking.org/ ZINC].
   
* it is advisable to create a special subdirectory, since many new    files will be generated.  
* the file containing the [http://www.daylight.com/smiles/ SMILES] strings should contain a string    followed by an identifier on each line.  
* OPTIONAL: run <tt>convert.py --i=yourname.smi --o=yourname.ism</tt> . This will  convert your SMILES to ''isomeric'' SMILES.
* run <tt>dbgen.csh yourname.smi</tt>.  
** Note that the dbgen.csh does not work for more that 1000 molecules.
** Brake up your molecules into chunks of 1000 and run dbgen on each chunk.
** Clean up your directory afterwards. dbgen.csh generates a lot of files that you do not need if it ran correctly.  
* you should obtain a file <tt>somename.db.gz</tt> .

#### optional
To increase the number of molecules that are written out for the database generation, copy the file $DOCK_BASE/data/omega.parm into the directory that dbgen.csh is going to be run in.
*At the end of the omega.parm file you will see a section called "Torsion Driving Parameters", here you will find three variables that can be changed.
***SetMaxConfs(600)   #set to higher numbers ie. 1000
***SetRMSThreshold(0.80)  #set to lower numbers ie. 0.50
***SetEnergyWindow(12.5)  #can be changed but this can often generate broken molecules
*WARNING this should only be done if generating conformations for a small set of compounds!!!

### Manual way

===Isolating the ligand as <tt>.mol2</tt>  file===
   
*extract the ligand structure from the <tt>.pdb</tt>  file.  
*assign hydrogens.  
*assign all atom ([http://www.tripos.com/mol2/atom_types.html Sybyl/TAFF]) and bond types.  
*save it as <tt>ligandname.mol2</tt>  file.  

#### Running <tt>omega</tt>
   
* run [http://www.eyesopen.com/products/applications/omega.html OMEGA], but don't ask me how to do that yet.

#### Running amsol

*find more information about amsol [http://comp.chem.umn.edu/amsol/ on its homepage].  
*<tt>mkdir ./amsol2</tt>   
*Use file2file.py to get the right formal charge to feed to AMSOL. It is also important to change the name, otherwise the original <tt>.mol2</tt> file will be overwritten!
<tt>file2file.py -g ligandname.mol2 ./amsol2/someothername.mol2</tt>      
*edit <tt>./amsol2/someothername.mol2</tt> :     
*
*delete all lines prior to <tt>@<TRIPOS>MOLECULE</tt>   
*
*change line 2 (molecule name) to something of the format <tt>ABCD12345678</tt> (four capital letters followed by eight numbers).  
*
*line 3 should be <tt>n<sub>atoms</sub> n<sub>bonds</sub> 0 0 0</tt>
*
*the <tt>@<TRIPOS>MOLECULE</tt>  section must consist of exactly '''5''' lines (adjust by adding/deleting blanks).  
*
*remove all sections after the <tt>@<TRIPOS>BOND</tt>  section.
*
*delete the blank lines between the <tt>ATOM</tt>  and <tt>BOND</tt>     sections, if there are any.    
*run <tt>RunAMSOL3.csh WAIT</tt>   
*the output <tt>someothername.solv</tt>  file will contain the following:
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

#### Running <tt>mol2db</tt>
   
*edit <tt>someothername.nmol2</tt>  so that the <tt>@<TRIPOS>MOLECULE</tt> section consists of exactly '''6''' lines.  
*edit the <tt>inhier</tt>  file so that the 'mol2_file',    'db_file' and 'solvation_table' entries are correct.  
*run <tt>mol2db inhier</tt>   
*add the preamble at the top of the file.  
*<tt>gzip</tt>  the resulting file so that it can be used by <tt>DOCK</tt> .

