# Written by Seth Vigneron, Shoichet Lab UCSF, 2022

# Takes input file of reactions and adds to an existing JSON file their associated SMARTS


import json
#import pickle # I believe pickle should not be used for security reasons as it would allow people to put in exicutable code
import rdkit
from rdkit.Chem import AllChem
import sys


# Load in rxn and bb input files
reactions_input_filename = sys.argv[1]
reactions_input_file = open(reactions_input_filename, 'r')
reactions_input_lines = reactions_input_file.readlines()
reactions_input_file.close

# Load in existing reactions JSON file
json_file_name = sys.argv[2]
json_file = open(json_file_name,)
Analoging_Reactions_Map = json.load(json_file) # Load in the Master Map dictionary
appended_json_file_name = sys.argv[3]


rxn_IDs_list = Analoging_Reactions_Map["Reaction_IDs"] # Have a list of all breakdown rxn IDs


for line in reactions_input_lines:
	splitline = line.strip().split("\t")

	rxn_name = splitline[0]
	rxn_ID = splitline[1]
	rxn_Description = splitline[2]
	rxn_build_SMARTS = splitline[3]
	rxn_synthonIDs = splitline[4]
	#rxn_build_SMARTS_obj = AllChem.ReactionFromSmarts(rxn_build_SMARTS)

	# Convert forward reaction SMARTS to breakdown reaction SMARTS
	rxn_build_SMARTS_split = rxn_build_SMARTS.split(">>")
	numReactants = len(rxn_build_SMARTS_split[0].split('.'))
	rxn_breakdown_SMARTS = "%s>>%s" % (rxn_build_SMARTS_split[1], rxn_build_SMARTS_split[0])
	#rxn_rxn_breakdown_SMARTS_obj = AllChem.ReactionFromSmarts(rxn_breakdown_SMARTS)

	rxn_IDs_list.append(rxn_ID)

	rxn_dict = {
		"Name" : rxn_name,
		"ID" : rxn_ID,
		"Description" : rxn_Description,
		"SynthonIDs" : rxn_synthonIDs,
		"Build_SMARTS" : rxn_build_SMARTS,
		"Breakdown_SMARTS" : rxn_breakdown_SMARTS,
		"NumReactants" : numReactants,
		}

	Analoging_Reactions_Map["Reactions"][rxn_ID] = rxn_dict
	

Analoging_Reactions_Map["Reaction_IDs"] = rxn_IDs_list



outfile_name = appended_json_file_name + ".json"
with open(outfile_name, 'w') as outfile:
	json.dump(Analoging_Reactions_Map, outfile)
outfile.close()

