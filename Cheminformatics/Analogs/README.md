# Analog Utilities and the Gradient Descent Pipeline
## Background:
The basis for these programs originates from the notion of “Chemical Space Travel.” In 2007, work published by van Deursen and Reymond [1] reported a “spaceship” program which travels from a starting molecule A to a target molecule B through a continuum of structural mutations, and thereby charts unexplored chemical space. To enable movement in an unexplored chemical space and the discovery of new structures, the authors describe chemical space as a structural continuum. Rather than referring to proximity in property space, the authors defined a finite set of “Nearest Neighbor Mutations” related through a single structural mutation : Atom Type Exchange, Atom Inversion, Atom Removal, Atom Addition, Bond Saturation, Bond Unsaturation, Bond Rearrangement, and Aromatic Ring Addition. This description organizes chemical space as a graph in which nodes represent molecules and edges represent mutations. In theory, one can go from any Molecule A to any Molecule B in a finite time by simply applying the correct series Nearest Neighbor Mutations to Molecule A sequentially. The original Analogs.py and Analog_Methods.py programs were built on top of RDKit and designed to perform this Nearest Neighbor Mutation analog procedure for molecules specified in a .SMI file.

Work Cited:

[1] van Deursen, R. and Reymond, J.-L. (2007), Chemical Space Travel. ChemMedChem, 2: 636-640. https://doi.org/10.1002/cmdc.200700021
## Analogs.py:
Usage: 
```
$python3 Analogs.py -i input.smi -o output
```
Analogs.py wraps Analog_Methods.py and serves as the location to which most user adjustments can be made. The most important user adjustment is the analog_methods list in main(). Here, all of the Nearest Neighbor Mutations are specified:
 ```
Trim_Extremities: Trim Parent Molecule Extremity Atoms One At A Time if M.W. > 500 Da
BO_Stepup: Recursively Increases Bond Order Until Maximum Conjugation
BO_Stepdown: Decreases Bond Order of Non-Single Bonds
Ring_Breaker: Enumerates Rings In Parent And Opens Rings
Ring_Maker: Enumerates Terminal -CH3, Finds Paths Of Length (4, 5) And Forms Rings With SP3 Carbons On Path
Walks: Performs Walks On Parent Molecule
Scans: Performs Scans On Parent Molecule
```
Note: A ‘Scan’ replaces hydrogens in a X-H bond with R, forming X-R. A ‘Walk’ replaces a heavy atom in a X-H bond with R, forming R-H.
Note: Some mutations (like some ‘Walks’ and ‘Scans’) produce more conservative changes when generating analogs, while other mutations (like ‘Ring_Breaker’) have the potential to generate analogs very different from the starting compound. To control for this, the analogging methods performed can be specified simply by commenting lines (using the ‘#’ character) containing unwanted operations in the analog_methods list. For example, if I only wanted the generator to produce analogs from halogen scans (F, Cl, Br, and I), my analog_methods list would look like: 	
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
Once the analog_methods list has been adjusted, the code can be run on the command line. In order for the code to run, do the following:
* Both Analogs.py and Analog_Methods.py need to be in the same directory
**As of 7/31/23, you can source the code here: /mnt/nfs/home/nlevinzon/analog_gen/new
* The code must be run in a location supporting Python3
**As of 7/31/23, I am using Gimel2
* A Python Environment containing RDKit must be sourced
**As of 7/31/23, you can source my Python Environment:
```
$source ./mnt/nfs/home/nlevinzon/bashrc.env
$conda activate venv
```
Once this has been completed, you can run the analog generator on the command line. The code behaves generally by going into the specified .SMI file, generating analogs for each line, generating stereoisomers for each analog (if applicable), and writing to a new .SMI file. For easier record keeping, the code will read the first molecule in the .SMI and generate all analogs for that first molecule before moving onto the second molecule. Analogs can be easily seen in the output .SMI file because their generated ZINC IDs will take the format “PARENT_ID_analog_number.” For example, if my .SMI input file contained the ZINC ID “ZINC184991516,” then all analogs generated from that compound would have the ZINC ID “ZINC184991516_analog0001” The “analog_number” counts up to enumerate how many analogs each parent compound produced. 
Currently, the code run on the command line in Gimel2 produces ~250 molecules/second.

## Analog_Methods.py
Analog_Methods.py contains the actual operations specified as Nearest Neighbor Mutation. Here you can adjust the specifics of each method currently implemented, as well add new ones. When developing new methods, the following paradigm should be used:
```
def new_method(smiles):
	"""Describe your new method"""

	# Convert the SMILES code to a molecule object
	mol = Chem.MolFromSmiles(smiles)
	if not mol:
    		raise ValueError(f"Invalid SMILES: {smiles}")

	# Specify new operations
	# …
	
	# Remove duplicates and the initial SMILES from the analogs list
	analogs = [mol for mol in analogs if mol != mol and Chem.MolToSmiles(mol) != smiles]

	# Return analogs as RWMol objects
	return analogs
```
## SMI_Filter.py
