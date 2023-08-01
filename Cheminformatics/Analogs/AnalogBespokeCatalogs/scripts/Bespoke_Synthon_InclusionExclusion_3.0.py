# Written By Seth Vigneron, Shoichet Lab UCSF, 2022
import json
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

#################################################

class Synthon_Characteristics:

	def __init__(self, name, maxHAC, allowSymm, inSMARTS, exSMARTS, rxnTag):
		self.name = name
		self.maxHAC = maxHAC
		self.allowSymm = allowSymm
		self.inSMARTS = inSMARTS
		self.exSMARTS = exSMARTS
		self.rxnTag = rxnTag



#################################################
# Determine if the current compound has the correct Reactive Center (RC) inclusion pattern
def inclusion_rules_check(synthon_rules, compound_mol):
#If the BB has a match to the inclusion Smarts, how many? if only 1; return true, elif if more than 1 and synthon allowed multiple symmetric synthons; return true, else return false
	for inclusion_SMARTS in synthon_rules.inSMARTS:
		inclusion_pattern_mol = Chem.MolFromSmarts(inclusion_SMARTS)
		if compound_mol.HasSubstructMatch(inclusion_pattern_mol, useChirality=True): 
			num_substructs = compound_mol.GetSubstructMatches(inclusion_pattern_mol,uniquify =True) 
			if len(num_substructs) == 1:
				return True
 
			elif synthon_rules.allowSymm:
				return allowed_reactive_symm_check(reagent_name, compound_mol)

			else:
				return False		
		else:
			continue

	return False
#################################################
# Determine if the current compound sidechains contain a specificed exclusion pattern
def exclusion_rules_check(synthon_rules, compound_mol):
#if the BB has a match to any of the exclusion SMARTS, return the text version of that SMARTS, else return a blank string
	for exclusion_pattern_mol in synthon_rules.exSMARTS:
		if compound_mol.HasSubstructMatch(exclusion_pattern_mol, useChirality=True):
			global pattern_text_lookup_dict
			matched_SMARTS = SMARTS_to_text_lookup_dict[exclusion_pattern_mol]
			#print(Chem.MolToSmiles(compound_mol),Chem.MolToSmiles(exclusion_pattern_mol))
			return matched_SMARTS

	return ''
	
#################################################
# Determine if the current compound contains, and is allowed to contain, multiple symmetric RC
def allowed_reactive_symm_check(compound_mol):
	symm = Chem.RemoveHs(compound_mol).GetSubstructMatches((Chem.RemoveHs(compound_mol)), useChirality =True, uniquify =False)
	return len(symm) == 2 or len(symm) == 4

#################################################

# Determine if the current compound pass the inclusion/exclusion rules defined in the rxn rules file
def filter_building_block(compound_line, synthon_selection_rules):
	# Collect info on the compound
	# Arthor BB TSV lines in the form: count \t smiles \s zincID \t \t purch_catagory 
	BB_count, BB_smi_ID, blank, BB_purchasability_phrase = compound_line.strip().split('\t')
	BB_smi, BB_ID  = BB_smi_ID.split(' ')
	
	# Test if BB is too expensive
	purch_code = int(purch_phrase_dict[BB_purchasability_phrase])
	if purch_code < min_purch_code:
		BB_purch_fail_line = '%s\t%s\tpurchasability\t%s\t%s\t%s\n' % (BB_smi, BB_ID, purch_code, min_purch_code, synthon_selection_rules.name ) 
		excluded_compounds_file.write(BB_purch_fail_line)	
		return

	# Convert BB to mol object
	BB_mol = Chem.AddHs(Chem.MolFromSmiles(BB_smi))
	
	# Test if BB is too large
	BB_HAC = BB_mol.GetNumHeavyAtoms()
	if BB_HAC > synthon_selection_rules.maxHAC:
		BB_HAC_HAC_fail_line = '%s\t%s\tHeavy_Atom_Count\t%s\t%s\n' % (BB_smi, BB_ID, BB_HAC, synthon_selection_rules.name)
		excluded_compounds_file.write(BB_HAC_HAC_fail_line)
		return

	# Test if BB has correct synthon, then if has any excluded moieties
	if inclusion_rules_check(synthon_selection_rules, BB_mol): # BB contains the inclusion pattern? Else, write to excluded_compounds with 
		matched_exlusion_pattern = exclusion_rules_check(synthon_selection_rules, BB_mol) # BB contains an exclusion pattern? If not, it passes! Else, write to excluded_compounds with matched exSMARTS
		if not matched_exlusion_pattern:
			BB_success_output_string = '%s\t%s\t%s\t%s\n' % (BB_smi, BB_ID, purch_code, synthon_selection_rules.rxnTag) # If pass all checks then write to the output file
			output_file.write(BB_success_output_string)
	
		else: # If BB matched to an exSMARTS, write to excluded_compounds with matched exSMARTS
			BB_exSMARTS_fail_line = '%s\t%s\texclusion\t%s\t%s\n' % (BB_smi, BB_ID, matched_exlusion_pattern, synthon_selection_rules.name) 
			excluded_compounds_file.write(BB_exSMARTS_fail_line)
	else: # If BB doesn't contain the inclusion SMARTS pattern, write to excluded_compounds
		BB_inSMARTS_fail_line =	'%s\t%s\tinclusion\t%s\n' % (BB_smi, BB_ID, synthon_selection_rules.name)
		excluded_compounds_file.write(BB_inSMARTS_fail_line)

#################################################
#def test_RSF(RSF, reagent_counts):
	
	

#################################################
# Dictionary associating Arthor BB purch catagories with their purch code value
purch_phrase_dict = {'Premier'  : 50 ,
                     'In-Stock' : 40 ,
                     'Agent'    : 30 ,
                     'On-Demand': 20 ,
                     'For-Sale' : 10 ,
                     'Annotated': 1  ,
                     'BB-All-For-Sale-19Q4-26M' : 40,
                     'BB-Stock-19Q4-1.6M' : 40,
                     'In-Stock-19Q4-13.8M' : 40,
					 'BB-50-21Q4-1.5M' : 50,
                     'BB-40-21Q4-590K' : 40,
                     'BB-30-21Q4-3K' : 30,
					'BB-10-21Q4-1.2M' : 10,
					 'BB-50-22Q1' : 50,
					 'BB-40-22Q1' : 40,
					 'BB-30-22Q1' : 30

}

#################################################
#bad_RSF, error_message = test_RSF(inclusion_exclusion_rules_lines, num_reagents)
#if bad_RSF:
#	print(error_message)
#	exit()


#####################################################################################################################

# Read in json RSF_Map File
json_file = open(sys.argv[1],)
RSF_Map = json.load(json_file)
# test_RSF(RSF_Map) # Check to ensure .json file was in correct format


Lib_ID = RSF_Map["Library_ID"]
synthon_name_list = RSF_Map["All_Synthons"]
min_purch_code = int(RSF_Map["Min_Purch_Code"])
SMARTS_to_text_lookup_dict = {}

# Get unfiltered BB files
BB_input_files = sys.argv[2:]

if len(BB_input_files) != len(RSF_Map["Synthon_Map"]):
	print("Error, different number of Synthon types than BB input files")
	exit()

# Iterate over all the provided unfiltered BB input files
for i in range(0, len(BB_input_files)):
	synthon_smi_input_file_name = sys.argv[i+2]
	synthon_smi_input_file = open(synthon_smi_input_file_name, 'r')
	synthon_BB_smi_lines = synthon_smi_input_file.readlines()
	synthon_name = synthon_name_list[i]
	synthon_chars_dict = RSF_Map["Synthon_Map"][i]

	if synthon_name != synthon_chars_dict["Synthon_Type"]:
		print("Error, input RSF synthons out of order")
		exit()
	
	# Collect info about synthon type, save as a class object
	synthon_max_HAC = int(synthon_chars_dict["Synthon_max_HAC"])
	synthon_allow_symm_RC = synthon_chars_dict["Synthon_Allowed_Symm"]
	synthon_inclusion_SMARTS = synthon_chars_dict["Synthon_SMARTS"]
	synthon_exclusion_SMARTS = synthon_chars_dict["Synthon_Exclusion_SMARTS"]
	synthon_output_file_name = synthon_chars_dict["Synton_Output_File_Name"]
	synthon_rxnTag = synthon_chars_dict["Synthon_Reaction_Tag"]

	# Associate the SMARTS text with the SMARTS mol obj to show which SMARTS called for BB exclusion
	synthon_exclusion_SMARTS_mol = []
	for exSMARTS in synthon_exclusion_SMARTS:
		exMol = Chem.MolFromSmarts(exSMARTS)
		SMARTS_to_text_lookup_dict[exMol] = exSMARTS
		synthon_exclusion_SMARTS_mol.append(exMol)

	# Create Synthon Rules class obj that holds common info/rules for the synthon type
	synthon_selection_rules = Synthon_Characteristics(synthon_name, synthon_max_HAC, synthon_allow_symm_RC, synthon_inclusion_SMARTS, synthon_exclusion_SMARTS_mol, synthon_rxnTag)

	# Open files for excluded BBs and successful BBs
	excluded_compounds_file_name = 'excluded_%s' % (synthon_output_file_name)
	excluded_compounds_file = open(excluded_compounds_file_name, 'w')
	output_file = open(synthon_output_file_name, 'w')

	# Iterate over all BBs in the unfiltered BB input file
	for building_block_line in synthon_BB_smi_lines:
		filter_building_block(building_block_line, synthon_selection_rules)

	# Close files
	excluded_compounds_file.close()
	synthon_smi_input_file.close()
	output_file.close()
