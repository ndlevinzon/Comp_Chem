# Analog Utilities and the Gradient Descent Pipeline
## Background:
The basis for these programs originates from the notion of “Chemical Space Travel.” In 2007, work published by van Deursen and Reymond [1] reported a “spaceship” program that travels from a starting molecule A to a target molecule B through a continuum of structural mutations, and thereby charts unexplored chemical space. The authors describe chemical space as a structural continuum to enable movement in an unexplored chemical space and the discovery of new structures. Rather than referring to proximity in property space, the authors defined a finite set of “Nearest Neighbor Mutations” related through a single structural mutation: Atom Type Exchange, Atom Inversion, Atom Removal, Atom Addition, Bond Saturation, Bond Unsaturation, Bond Rearrangement, and Aromatic Ring Addition. This description organizes chemical space as a graph in which nodes represent molecules and edges represent mutations. Theoretically, one can go from Molecule A to Molecule B in a finite time by simply applying the correct series of Nearest Neighbor Mutations to Molecule A sequentially. The original Analogs.py and Analog_Methods.py programs were built on top of RDKit and designed to perform this Nearest Neighbor Mutation analog procedure for molecules specified in a .SMI file.

Work Cited:

[1] van Deursen, R. and Reymond, J.-L. (2007), Chemical Space Travel. ChemMedChem, 2: 636-640. https://doi.org/10.1002/cmdc.200700021

# The Analog Pipeline
Generating analogs from an input .SMI can be performed using either Analog_Generator_Dev.py (requiring users to change input files within the code) or Analogs_Run.py and Analogs_Methods.py in concert (requiring users to input files on the command line). Both files must be in the same directory if the latter option is used.

## Analogs_Run.py:
The Python environment must be sourced for the Python code to access the correct dependencies. As of 7/31/23, you can source my Python Environment:
```
>>source ./mnt/nfs/home/nlevinzon/bashrc.env
>>conda activate venv
```
Usage: 
```
>>python3 Analogs_Run.py -i input.smi
```
Analogs_Run.py wraps Analog_Methods.py and serves as the location where most user adjustments can be made. The most important user adjustment is the analog_methods list in main(). Here, all of the Nearest Neighbor Mutations are specified:
 ```
Trim_Extremities: Trim Parent Molecule Extremity Atoms One At A Time if M.W. > 500 Da
BO_Stepup: Recursively Increases Bond Order Until Maximum Conjugation
BO_Stepdown: Decreases Bond Order of Non-Single Bonds
Ring_Breaker: Enumerates Rings In Parent And Opens Rings
Ring_Maker: Enumerates Terminal -CH3, Finds Paths Of Length (4, 5) And Forms Rings With SP3 Carbons On Path
Walks: Performs Walks On Parent Molecule
Scans: Performs Scans On Parent Molecule
```
* Note: A ‘Scan’ replaces hydrogens in a X-H bond with R, forming X-R. A ‘Walk’ replaces a heavy atom in a X-H bond with R, forming R-H.
* Note: Some mutations (like some ‘Walks’ and ‘Scans’) produce more conservative changes when generating analogs, while other mutations (like ‘Ring_Breaker’) have the potential to generate analogs very different from the starting compound.
To control this, the analogging methods performed can be specified simply by commenting lines (using the ‘#’ character) containing unwanted operations in the analog_methods list. For example, if I only wanted the generator to produce analogs from halogen scans (F, Cl, Br, and I), my analog_methods list would look like: 	
```
analog_methods = [
    	# [trim_extremities, "trim"],
    	# [BO_stepup, "increase-bond-order"],
    	# [BO_stepdown, "decrease-bond-order"],
    	# [ring_breaker, "ring-opening"],
    	# [ring_maker, "ring-closure"],
    	# [nitrogen_walks, "n-walk"],
    	# [oxygen_walks, "o-walk"],
    	# [sulfur_walks, "s-walk"],
    	# [heterocycle_walks, "heterocycle-walk"],
    	# [methyl_scanning, "ch3-scan"],
    	# [amine_scanning, "nh2-scan"],
    	# [hydroxyl_scanning, "oh-scan"],
    	# [thiol_scanning, "sh-scan"],
    	[fluorine_scanning, "f-scan"],
    	[chlorine_scanning, "cl-scan"],
    	[bromine_scanning, "br-scan"],
    	[iodine_scanning, "i-scan"],
    	# [bioisosters_scanning, "bioisoster-scan"],
	]
```

Once the analog_methods list has been adjusted, the code can be run on the command line. For the code to run, do the following:
* Both Analogs.py and Analog_Methods.py need to be in the same directory. As of 7/31/23, you can source the code here: /mnt/nfs/home/nlevinzon/analog_gen/new
* The code must be run in a location supporting Python3. As of 7/31/23, I am using Gimel2
* A Python Environment containing RDKit must be sourced

Once this has been completed, you can run the analog generator on the command line. The code behaves generally by going into the specified .SMI file, generating analogs for each line, generating stereoisomers for each analog (if applicable), and writing to a new .SMI file. For easier record keeping, the code will read the first molecule in the .SMI and generate all analogs for that first molecule before moving onto the second molecule. Analogs can be easily seen in the output .SMI file because their generated ZINC IDs will take the format “PARENT_ID_analog_number.” For example, if my .SMI input file contained the ZINC ID “ZINC184991516,” then all analogs generated from that compound would have the ZINC ID “ZINC184991516_analog0001” The “analog_number” counts up to enumerate how many analogs each parent compound produced. 
Currently, the code run on the command line in Gimel2 produces ~250 molecules/second.

## Analog_Methods.py
Analog_Methods.py contains the actual operations specified as Nearest Neighbor Mutation. Here you can adjust the specifics of each method currently implemented and add new ones. When developing new methods, the following framework should be used:
```
def new_method(smiles):
	"""Describe your new method"""

	# Convert the SMILES code to a molecule object
	mol = Chem.MolFromSmiles(smiles)
	if not mol:
    		return []

	analogs = []

	# Specify new operations
	# …
	
	# Remove duplicates and the initial SMILES from the analogs list
	analogs = [mol for mol in analogs if mol != mol and Chem.MolToSmiles(mol) != smiles]

	# Return analogs as RWMol objects
	return analogs
```
## Preparing for Docking: 3D Building Pipeline and DB2 Generation
First, you will need the correct environment:
```
source /nfs/soft/dock/versions/dock38/pipeline_3D_ligands/env.(sh|csh)
```

This environment will set up most of the required variables for you, as well as adds the submission scripts to your PATH, which means submission can be as simple as:

bash:
```
export INPUT_FILE=$HOME/myligands.smi
export OUTPUT_DEST=$HOME/myoutput
submit-all-jobs-slurm.bash
```
csh:
```
setenv INPUT_FILE $HOME/myligands.smi
setenv OUTPUT_DEST $HOME/myoutput
submit-all-jobs-slurm.bash
```

After building, we must make, extract, and generate a file named "job_input_list" from a series of .DB2 files:
```
>> make_tarballs.bash </OUT Directory> <One Level Above /OUT Directory>
>> tar xf *.db2.tar.gz
>> find $PWD -type f -name "*.db2*" > my_db2_list
>> split -a 3 --lines=5000 my_db2_list db2_chunk.
>> for db2_chunk in db2_chunk.*; do
>>	tar -czf $db2_chunk.db2.tgz --files-from $db2_chunk
>>	done
>> find $PWD -type f -name "db2_chunk.*.db2.tgz" > job_input_list
```
Use job_input_list, which contains paths to the correct .DB2 archives, as the input for the docking script

## Extracting and Analyzing OUTDOCK
The first thing you will need to do is format the OUTDOCK into a .CSV file containing the following standard headers:
```
|Dir_Source|mol#|id_num|flexiblecode|matched|nscored|time|hac|setnum|matnum|rank|charge|elect|gist|vdW|psol|asol|tStrain|mStrain|rec_d|r_hyd|Total|
```
To do this, I first use process_dock_output.py to generate a dir_list file. Then I use a Python environment running Python2.7 in order to run extract_all_blazing_fast.py on my dir_list:
```
>> process_dock_output.py
>> python2 extract_all_blazing_fast.py ./dir_list
```
From here, I usually import the data into Microsoft Excel to generate the correctly-formatted .CSV (but I am sure there is a more elegant way to do this).
Finally, run ligand_analog_curves and 





## SMI_Filter.py
This application is built on top of ChemAxon's CXCALC to take in a .SMI input and calculate each entry's pKa/b and Lipinski RO5
Usage: 
```
>>python3 smi_filter.py -i input.smi
```
Before running this code, ensure that the lines interacting with the CLI point to a directory containing ChemAxon. Also, several temporary files will be generated to parse the .SMI and CXCALC outputs (ensure that your working directory has space and RW permissions)

## fix_truncated_zinc.py
When building from an .SMI in the UCSF DOCK 3D building pipeline, sometimes the ZINC IDs get truncated when printing the OUTDOCK. These impartial ZINC IDs make grouping ligands by family difficult after docking is complete. To solve these partial ZINC_IDs, a "fuzzy sort" is implemented to replace the incomplete ZINC IDs from the OUTDOCK using the output of Analogs.py as a key.
Usage: 
```
>>python3 fix_broken_zinc.py -i OUTDOCK.csv -k KEY.smi
```
