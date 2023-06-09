PyRED Website: https://upjv.q4md-forcefieldtools.org/REDServer-Development/ 

# Prepare PDB input file(s)
Make sure to name the files correctly (i.e. Mol_red$n.pdb where $n = the molecule number (1, 2, …)) 
You can upload several PDB files just make sure naming is correct. One PDB file can contain several conformations but make sure each PDB file corresponds to a single molecule. PDB’s with multiple conformations will produce a <.log> file corresponding to each optimized geometry.
Make sure charge/multiplicity is correct
If not in the PDB file(s), can edit through the addition of a <.config> file
Create PDB for PyRED input:
cp 5hnj.pdb Mol_red1.pdb
# Generate “Project.config” and/or “System.config” input files. 
THESE FILES ARE NOT NEEDED IF THE MULTIPLICITY = 1 AND CHARGE = 0
If this file is absent, default tasks are executed.
System.config : information related to tasks performed by PyRED and allows for overwriting default tasks (geometry optimization, charge fitting, etc.…).
Project.config : information about the molecules involved in the project and allows for modifying of default information (charge, multiplicity, etc.…). 
PyRED will optimize your structure’s geometry unless you tell it not to. If you’ve already optimized in Gaussian, you can save time by turning this function off (Example: OPT_Calc = OFF1).
You will also need to add in your structure’s multiplicity/charge (Example: MOLECULE1 -TOTCHARGE 2+).
You can do A LOT with this file and there is really helpful documentation on PyRED’s Mini-How-To: https://upjv.q4md-forcefieldtools.org/Tutorial/Mini-HowTo-InputFiles.pdf
System.config input created in text editor of choice:
OPT_Calc = OFF1
MOLECULE1 -TOTCHARGE 2+

As more PyRED jobs are performed on a system, you may need to define forcefield parameters using a “frcmod.user” file. This of course is not needed if you are just generating <.mol2> files for AMBER purposes. 
Compress all files for submission into a single “archive” folder (i.e. <.tar> folder) for PyRED submission.
Compress to a .tar folder with all the necessary files using the following the command line format of: tar -cvzf <Archive.tgz> <files-you-want-to-compress>
tar -cvzf Archive.tgz Mol_red1.pdb Project.config System.config
Register as a new user on the PyRED server – makes it so you will get an email when your job starts/ends and you will get a personalized link to download your project files.
Submit your job through the PyRED online server through a “private account” – the account you just created.

Here is where you name your project. If you are wanting to use PyRED to optimize geometry, here is where you select what QM interface to use. There are multiple Gaussian versions available amongst a few other programs.

Here is a summary of your project information – a chance to double check everything is correct – and where you upload your archive folder with your PDB and <.config> files.

A link to be able to download your project files once the job is complete is provided, but this will also be in the email notifying job completion.
Go to the link provided in the email sent when your job completes and login to view your project files. 
In the emails, you may get a general sense if your job was completed successfully or not but peruse through the files to double check.
You have 20 days to download your project files before they are removed from the server.
Download all the project files to your local system and have fun with your new <.mol2> file!
Your project folder will be named after the job number ID and you will have to unzip it to access all your new files.
Your <.mol2> file will be in Data-R.E.D.Server directory within the project folder
To unzip your project file, follow this command line format: tar -xzvf <folder-you-want-to-unzip.tar.gz>
tar -xvf P7726.tar.gz
cd P7726/Data-R.E.D.Server/Mol_m1
more Mol-sm_m1-c1.mol2

GENERATING FRCMOD + LIBRARY FILES
1. <.frcmod>
    	parmchk2 -i m1A.mol2 -f mol2 -o m1A.frcmod
2. <.lib>
    	tleap
    	source leaprc.RNA.OL3
    	source leaprc.gaff
    	loadamberparams m1A.frcmod
    	m1A = loadmol2 m1A.mol2
    	check m1A
    	saveoff m1A m1A.lib
    	quit (edited)

