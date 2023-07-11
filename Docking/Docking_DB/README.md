==Automated Docking Database Tools==

#[[#Automatic Database Generation]]: You want to generate DOCK 3.5 hierarchy databases for your molecules
#[[#Automatic Decoy Generation]]: You want to generate DUD style decoys property matched to your active ligands

==Automatic Database Generation==

If you just want to make a docking database for some molecules, then this is the section for you.

===Simple Database Generation===
*For automated database generation on small input files (say < 5000 molecules). Making a new directory is a good idea because these scripts generate a LOT of output files.
 mkdir new_dir_name
 cd new_dir_name
 dbgen.csh INPUT
 < or >
 dbgen.csh INPUT [PROTONATION]

*Setup:
Scripts in this section are automatically put in your path by the DOCK login scripts. If they are not, then inside a csh first do the following:
 setenv DOCK_BASE /raid1/soft/dockenv
 source $DOCK_BASE/etc/login

*Options:
You can use dbgen.csh alone to get help. INPUT is a file containing the input molecules. Either use a .smi file containing lines of smiles strings and ids or some other file type easily converted to smiles (i.e. multi .mol2 or .sdf). The optional PROTONATION argument can be used to generate databases containing extended protonation states. The available protonation types are as follows:
#ref - only the reference protonation
#mid - reference plus middle protonation [default]
#lo  - reference, middle, and lo protonation
#hi  - reference, middle, and hi protonation
#all - all protonation ranges

*Caveats:
dbgen.csh is most useful when you want to test out a dockable database without worrying about ZINC. If you like the molecules and decide to add them to ZINC, it should be easy using the output of dbgen.csh. Please contact me (Michael Mysinger) if you want to do this at any time, as it should be easy but is untested at the moment. If you want to add the molecules to ZINC from the beginning then you can use the XML-RPC interface of DOCKBlaster like so:
 xmlclient.py upload my.smi   # uploads ligands to server
 xmlclient.py qup ID          # later on, get docking database back
where ID is the job id returned by the upload command.

===Complex Database Generation===

*If you want to do automated database generation on a large scale (> 5000 molecules), then look here. First, <i>NOTA BENE</i> that this process is demanding and has been known to fill all space on the file servers, slam them into submission, or overload the entire SGE cluster. For estimation purposes, assume the processes take ~40GB of disk per 100k molecules.
 ssh sgehead.compbio.ucsf.edu                 # go to cluster head node
 set dud=~mysinger/code/dud/trunk             # setup variable for path
 $dud/dbmake.py -i my.smi |& tee dbmake.log   # run the large scale database generation

Options:
Use -n to change the name of the output database. Use -p to change the protonation type to one of the protonation options listed for dbgen.csh above.

If instead you want to load the molecules directly into ZINC, then John's corresponding script is 'mas.csh SETNAME INFILE'.

==Automatic Decoy Generation==

If you just want DUD style computational decoys property matched specifically to your ligands, then look here. For the fastest turnaround given the molecules already have ZINC ids, use decoys.py. To regenerate everything from scratch like [[#Automatic Database Generation]] above, use dudgen.py.

*Generate DUD style decoys starting from known ZINC ids, pulling the pre-generated hierarchy databases directly from ZINC

 set dud=~mysinger/code/dud/trunk             # setup variable for path
 $dud/decoys.py -i INFILE                     # make DUD style decoys
It pulls ZINC IDs from the first column of stdin (or INFILE if -i is used), but the id column can be changed with -c COL. Use -n X to change the number of decoys generated per ligand from 35 to X. Use -s to output selected decoys but skip final database assembly.  

*Generate DUD style decoys starting from smiles. This script first builds the ligands, then parses them to get the ligand properties,  the selecting and building 

The ids in the smiles file must be either ZINC ids or numbers. When not using ZINC ids the ligands and decoys will have the same ids (P000001, P000002, ...), so you will need to edit the ligands database file and change the ids.

 ssh sgehead.compbio.ucsf.edu                 # go to cluster head node
 $dud/dudgen.py -i INFILE                     # Fully rebuild ligands and DUD style decoys
This last script could probably have many more command line options, so let me know 1) if you use it and 2) if you think of something you would like to change.

Note that this script has some specific naming requierments:
For a ligand name with CHEMBL1200666 results in the following error:

  Error: invalid literal for int() with base 10: 'BL120066'
  
  Likely caused by my rudimentary handling of molecule ids.
  They should look like TEMP12345678, P12345678, or 12345678
  Sorry...

It is recomended that you not give this script smiles with names that excede 9 characters sence
DOCK3.6 will truncate the molcule name if it is grater than 9 characters long.  

here is a awk script to rename the smiles file:
  
   awk 'BEGIN{count=1}{printf "%s TEMP%d\n",$1,count; count=count+1}' ligands.smi > ligands.mod.smi

To generate decoys and ligands for DOCK3.7, you could run the these scripts which generated smiles file for the decoys.  you can then run this smiles file with [[Ligand preparation]].

&rarr; Back to [[Tutorials]]

[[Category:Tutorials]]
