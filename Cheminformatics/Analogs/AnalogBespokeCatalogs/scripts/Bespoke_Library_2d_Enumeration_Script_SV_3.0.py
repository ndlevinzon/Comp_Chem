# Written by Seth Vigneron, Shoichet Lab UCSF, 2021
# Takes input files of smiles for each reagent Building Block (BB) type
# Loops through every combination of BBs (one from each list)
# Runs all the Reaction SMARTs that should be used for each specific combination of BBs as defined in the .json input file
# Builds the 2D structure of the products to be formed
# Checks to see if the product passes disqualification checks based on the disqualification critera defined in the .json input file
# If the product passes, writes it to the output file and then moves onto the next set of BBs
import sys
import json
import itertools
import Bespoke_Disqualification_Criterea
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

sys.path.append(".")

#####################################################################################################################
#####################################################################################################################

# Class Objects for Building Blocks, Product Molecules, and Reactions
#####################################################################################################################

class BuildingBlock:

	def __init__(self, smiles, bbid, purch, rxn_tag, mol):
		self.smiles = smiles
		self.bbid = bbid
		self.purch = purch
		self.rxn_tag = rxn_tag
		self.mol = mol

#######################################

class Product_Identifiers:

	def __init__(self, mol, smiles, purch, tranche, BB_IDs, BB_smiles, code_descriptor):
		self.mol = mol
		self.smiles = smiles
		self.purch = purch
		self.tranche = tranche
		self.BB_IDs = BB_IDs
		self.BB_smiles = BB_smiles
		self.code_descriptor = code_descriptor

	def BB_component_IDs(self):
		components_ID = ''
		for bb in self.BB_IDs:
			components_ID += '_'
			components_ID += bb
		return components_ID 

	def BB_component_Smiles(self):
		components_smiles = ''
		for bb in self.BB_smiles:
			components_smiles += bb
			components_smiles += '\t'
		return components_smiles

#######################################	

class Reaction_Identifiers:

	def __init__(self, rxn, Product_Code_Descriptor):
		self.rxn = rxn
		self.Product_Code_Descriptor = Product_Code_Descriptor

#######################################

class Restrictions_Identifiers:

	def __init__(self, mol, BB_combo):
		self.mol = mol
		self.BB_combo = BB_combo

#####################################################################################################################
#####################################################################################################################

# Functions
#####################################################################################################################
# Store each line of the BB input files as a BB class object
def extract_bbs(file_name):

	reagent_bb_objects = []

	reagent_file = open(file_name, 'r')
	reagent_file_lines = reagent_file.readlines()

	for bb_line in reagent_file_lines:
		bb_char_list = bb_line.strip().split('\t')
		bb_smiles = bb_char_list[0]
		bb_ID = get_BB_ID(bb_char_list[1])
		bb_purch = bb_char_list[2]
		bb_rxn_tag = bb_char_list[3]
		bb_mol = Chem.AddHs(Chem.MolFromSmiles(bb_smiles))

		reagent_bb_objects.append(BuildingBlock(bb_smiles, bb_ID, bb_purch, bb_rxn_tag, bb_mol))
	
	reagent_file.close()
	return reagent_bb_objects

#######################################
# Remove the ZINC00... from the ZincIDs to save space
def get_BB_ID(zincID):
	count = 7
	while zincID[count] == '0':
		count += 1
	return zincID[count:]

#######################################
# Get the built product molecules from a given combination of BBs and write them to the output file
def Build_Products_and_Outputs(BB_combo):
		
	All_Products = Determine_and_Run_Reactions(BB_combo)

	for product in All_Products:
		Product_code = "%s_%s%s_%s" % (product.purch, Lib_ID, product.BB_component_IDs(), product.code_descriptor)
		Output_line = "%s\t%s\t%s\t%s\n" % (Product_code, product.smiles, product.tranche, product.BB_component_Smiles())
		Product_library.write(Output_line)

#######################################
# Check if the reactions for a given combination of BBs have already been run, skip this BB combo if so 
def BB_combo_chars_and_diqualification_check(BB_combo):
	disqualification = False # Default is the current BB combo is not disqualified
	
	BB_ID_string = ''
	for bb in BB_combo:
		BB_ID_string += bb.bbid

	return BB_ID_string, disqualification

#######################################
# Calculate the product molecule's purchasability tranche based on the purchasability of the BBs
def get_purchasability(BB_list):
	total_bbs = len(BB_list)
	max_purch_count = 50 * total_bbs
	purch_count = 0
	for bb_purch in BB_list:
		purch_count += int(bb_purch)

	purch_percent = purch_count / max_purch_count
	if purch_percent == 1:
		product_purchasability = 1
	elif purch_percent >= 0.9:
		product_purchasability = 2
	elif purch_percent >= 0.8:
		product_purchasability = 3
	elif purch_percent >= 0.7:
		product_purchasability = 4
	else:
		product_purchasability = 5

	return product_purchasability

#######################################
# Take the given combination of BBs and return a list of product molecules
def Determine_and_Run_Reactions(BB_combo):
	
	BB_combo_rxn_tag = ''
	products_list = []
	
	product_BB_mol_list = []
	product_BB_smi_list = []
	product_BB_ids_list = []
	product_BB_purch_list = []

	# Get the rxn tag for this combination of BBs
	for BB in BB_combo:
		BB_combo_rxn_tag += str(BB.rxn_tag)
		product_BB_mol_list.append(BB.mol)
		product_BB_smi_list.append(BB.smiles)
		product_BB_purch_list.append(BB.purch)
		product_BB_ids_list.append(BB.bbid)

	# Run each reaction associated with the BB combo rxn tag
	for reaction in rxn_tag_to_rxnClass_dict[BB_combo_rxn_tag]:
		rxn = reaction.rxn
		product_code_descriptor = reaction.Product_Code_Descriptor

		product_mol = rxn.RunReactants(tuple(product_BB_mol_list))[0][0]	# Run Rxn on reactants
		
		# Check if product molecule passes all the given disqualification criterea
		product_char_object = Restrictions_Identifiers(product_mol, BB_combo)
		product_diqualification = Product_Char_and_Disqualification_Check(product_char_object)
		if product_diqualification:
			continue
		
		# If product molecule passes disqualification criterea, collect its identifiers, save it to a list of products built from this set of BBs
		product_tranche = get_tranche(product_mol)	
		product_smiles = Chem.MolToSmiles(Chem.RemoveHs(product_mol))
		product_purch = get_purchasability(product_BB_purch_list)
		products_list.append(Product_Identifiers(product_mol, product_smiles, product_purch, product_tranche, product_BB_ids_list, product_BB_smi_list, product_code_descriptor))

	return products_list	# return the list of all product molecules built from the provided BB combo 

#######################################
# Check if Product is Disqualified Based on Critera Provided in the json file input
def Product_Char_and_Disqualification_Check(product_object):
	disqualification = False  # Default is the current product is not disqualified

	# For each restriction listed in the json file, run the associated restriction class function to check if the product molecule passess 
	for product_restriction in Bespoke_Lib_Map["Product_Restrictions_Map"].keys():
		restriction_function_call_string = product_restriction + "_check"
		current_restriction_function = getattr(Bespoke_Disqualification_Criterea, restriction_function_call_string)
		disqualification = current_restriction_function(product_object)
		if disqualification:
			return disqualification
			
	return disqualification # If product molecule passes all the criterea, return the disqualification = False

#######################################
# Get the product molecules tranche based on its Heavy Atom Count and LogP
def get_tranche(mol):
	Product_HAC = mol.GetNumHeavyAtoms()
	logP = Chem.Crippen.MolLogP(Chem.RemoveHs(mol))
	tranche_logP, sign = scale_logp_value(logP)

	if Product_HAC > 99:
		Product_HAC = 99

	tranche_code = "H%02d%s%03d" % (Product_HAC, sign, tranche_logP)

	return tranche_code

#######################################
# Round LogP to Tranche Format
def scale_logp_value(logp):
	if logp < -9.0:
		logp = -9.0
	elif logp > 9.0:
		logp = 9.0
	if logp < 0.0 or logp >= 5.0:
		logp = 100*int(logp)
	else:
		logp = 10*int(10*logp)

	sign = 'M' if logp < 0.0 else 'P'

	tranche_logP = abs(logp)

	return tranche_logP, sign

#######################################
# Create Reaction Class Objects for every rxn SMARTs listed in the json file
def create_rxnClasses(json_Map):
	for rxn_type in Bespoke_Lib_Map["Reaction_Map"]:
		reactions_list = []

		reaction_tag = rxn_type["Reaction_Tag"]
		
		for i in range(0, len(rxn_type["Reaction_SMARTs"])):
			rxnSMARTs = rxn_type["Reaction_SMARTs"][i]
			rxn = AllChem.ReactionFromSmarts(rxnSMARTs)
			prod_code_descriptor = rxn_type["Product_Code_Descriptor"][i]
			
			reactions_list.append(Reaction_Identifiers(rxn, prod_code_descriptor))

		rxn_tag_to_rxnClass_dict[reaction_tag] = reactions_list # Populate a dictionary of all the reaction class objects associated with each rxn tag

#####################################################################################################################
#####################################################################################################################


#####################################################################################################################

# Get the Bespoke Library Map Json File
json_file = open(sys.argv[1],)
Bespoke_Lib_Map = json.load(json_file)


rxn_tag_to_rxnClass_dict = {}
create_rxnClasses(Bespoke_Lib_Map) # Create Reaction Class Objects for the provided rxn SMARTs

Lib_ID = Bespoke_Lib_Map["Library_ID"] # Get the Bespoke Library ID

# For each of the disqualification criterea listed in json file, create a class variable that equals the specified cut off value
for restriction_type in Bespoke_Lib_Map["Product_Restrictions_Map"].keys():
	setattr(Restrictions_Identifiers, restriction_type, Bespoke_Lib_Map["Product_Restrictions_Map"][restriction_type])

# Get the BB Input Files
Reagent_input_files = sys.argv[2:-1]

# Get the name of the output libraries 
Product_library_file_name = sys.argv[-1]
Product_library = open(Product_library_file_name, 'w')

all_BBs_objects_list = []
for reagent_file in Reagent_input_files:
	all_BBs_objects_list.append(extract_bbs(reagent_file))

# Get every combination of a single BB from each BB input file and begin running the 2D enumeration
for BB_combo_objects in itertools.product(*all_BBs_objects_list):
	Build_Products_and_Outputs(BB_combo_objects)

Product_library.close()

