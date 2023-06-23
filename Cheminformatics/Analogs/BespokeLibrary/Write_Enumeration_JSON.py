# Writen By Seth Vigneron, Shoichet Lab, 2021
# Takes a text file that outlines Reaction SMARTs and Product Disqualification Criterea
# and generates a .json file that can be read by the Bespoke 2D Enumeration Script


import json
import sys

#######################################
class Reaction_Characteristics:
	
	def __init__(self, SMARTs, RxnTag, ProdIDDescriptor):
		self.rxnSMARTs = SMARTs
		self.rxnTag = RxnTag
		self.descriptor = ProdIDDescriptor

#######################################
# Get the value specified a tab after the substring
def substring_index(lines_list, substring):
	line_index = [index for index, s in enumerate(lines_list) if substring in s][0]
	line = lines_list[line_index].strip().split('\t')
	return line[1]


#######################################

bespoke_lib_map = {} # Create the Master Map dictionary

reactions_input_file = open(sys.argv[1], 'r')
rxn_lines = reactions_input_file.readlines()

#######################################

# Get the bespoke library ID
LibID = substring_index(rxn_lines, "Bespoke Library ID:")
bespoke_lib_map["Library_ID"] = LibID

#######################################
# Collect all the Product Disqualification Criterea listed in the text file input
bespoke_lib_map["Product_Restrictions_Map"] = []

start_product_restrictions_line_Index = rxn_lines.index("Product Restrictions:\n") +1
end_product_restrictions_line_Index = rxn_lines.index("end\n", start_product_restrictions_line_Index)

# Create a dictionary of the Restriction Type (ie Max LogP) and its associate cutoff value
restrictions_type_dict = {}
for i in range(start_product_restrictions_line_Index, end_product_restrictions_line_Index):
	restriction_line = rxn_lines[i].strip().split('\t')
	restriction_type = restriction_line[0]
	restriction_value = restriction_line[1]
	restrictions_type_dict[restriction_type] = restriction_value

bespoke_lib_map["Product_Restrictions_Map"] = restrictions_type_dict # Add the restriction dict to the Master Map


#######################################
# Collect all the Reaction Smarts and their associated identifiers
bespoke_lib_map["Reaction_Map"] = []

start_reactions_line_Index = rxn_lines.index("Reactions:\n") +1
end_reactions_line_Index = rxn_lines.index("end\n", start_reactions_line_Index)

rxnTag_list = []
rxnTag_to_Reaction_dict = {}

# Loop through the lines in the input text file that specify the reaction SMARTs and its identifiers
for i in range(start_reactions_line_Index, end_reactions_line_Index):
	current_rxn_line = rxn_lines[i].strip().split('\t')
	current_rxn_SMARTs = current_rxn_line[0]
	current_rxn_tag = current_rxn_line[1]
	current_prodID_descriptor = current_rxn_line[2]
	
	# Create class object for each reaction listed in input text file that consists of the reaction and its identifiers
	current_reaction = Reaction_Characteristics(current_rxn_SMARTs, current_rxn_tag, current_prodID_descriptor)

	# Create a dictionary of each rxn tag and all the reaction class objects that are associated with that tag
	if current_rxn_tag not in rxnTag_list:
		rxnTag_list.append(current_rxn_tag)
	try:
		reactions_list = rxnTag_to_Reaction_dict[current_rxn_tag]
		reactions_list.append(current_reaction)
		rxnTag_to_Reaction_dict[current_rxn_tag] = reactions_list
	except KeyError:
		rxnTag_to_Reaction_dict[current_rxn_tag] = [current_reaction]

# For each rxn tag, append the Reaction Map section of the Master Map with lists of the reactions to run and what their identifiers are
for Tag in rxnTag_list:
	list_of_reactions = rxnTag_to_Reaction_dict[Tag]
	
	rxn_SMARTs_list = []
	Product_Code_Descriptor_list = []
	for reaction in list_of_reactions:
		rxn_SMARTs_list.append(reaction.rxnSMARTs)
		Product_Code_Descriptor_list.append(reaction.descriptor)
	
	bespoke_lib_map["Reaction_Map"].append({
		"Reaction_Tag" : Tag,
		"Reaction_SMARTs" : rxn_SMARTs_list,
		"Product_Code_Descriptor" : Product_Code_Descriptor_list
		 })

# Write .json file
outfile_name = '%s.json' % (sys.argv[1][:-4])
with open(outfile_name, 'w') as outfile:
	json.dump(bespoke_lib_map, outfile)

