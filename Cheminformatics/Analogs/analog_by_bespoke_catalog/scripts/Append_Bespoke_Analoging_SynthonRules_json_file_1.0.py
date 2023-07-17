# Written by Seth Vigneron, Shoichet Lab UCSF, 2022

# Takes input file of building block synthon types with their inSMARTS and exSMARTS, 
# Appends them to an existing json file to be used as the input for Analog_Breakdownrxns script


import json
import rdkit
from rdkit import Chem as C
from rdkit.Chem import AllChem
import sys

#######################################

class Synthon_Chars: # Class to hold all the info for each synthon ID (both CID and NID are stored as the same obj)
	
	def __init__(self, SynthClass, SynthID, inSMARTS, exSMARTS):
		self.SynthClass = SynthClass
		self.SynthID = SynthID
		self.inSMARTS = inSMARTS
		self.exSMARTS = exSMARTS

######################################
def extract_synthon_rules(in_line, ex_line): # Takes 2 lines, one for inSMARTS and one for exSMARTS

	global failed_synthon_IDs # Global list to keep track of failed input lines

	# Split the input lines seperatly
	in_splitline = in_line.strip().split("\t")
	ex_splitline = ex_line.strip().split("\t")

	# Collect info of the inSMARTS line
	in_synthon_class_fullname = in_splitline[0]
	in_synthon_class_ID_long = in_splitline[1]
	if not "C" in in_synthon_class_ID_long: # 'ShouldContain' inSMARTS IDs should start with C
		failed_synthon_IDs.append(in_synthon_class_ID_long)
		return
	in_synthon_class_ID = in_synthon_class_ID_long[1:]
	in_synthon_class_rules = in_splitline[2:]

	# Get the name of the Synthon class and subclass for these 2 lines; SynthonClass__SynthonSubClass__ShouldContainAtLeastOne (ex. Acid__Aliphatic_Acid__ShouldContainAtLeastOne__-)
	synthon_class_name = "__".join([in_synthon_class_fullname.split('__')[i] for i in [0,1]])

	# Collect info of the exSMARTS
	ex_synthon_class_fullname = ex_splitline[0]
	ex_synthon_class_ID_long = ex_splitline[1]
	if not "N" in ex_synthon_class_ID_long: # 'shouldNotContain' exSMARTS IDs should start with N
		failed_synthon_IDs.append(ex_synthon_class_ID_long)
		return
	ex_synthon_class_ID = ex_synthon_class_ID_long[1:]
	ex_synthon_class_rules = ex_splitline[2:]

	# Make sure both lines are representing the same Synthon ID
	if in_synthon_class_ID == ex_synthon_class_ID:
		# If so, store as a Synthon_Chars obj and return it
		return Synthon_Chars(synthon_class_name, in_synthon_class_ID, in_synthon_class_rules, ex_synthon_class_rules)		

	else:
		# If not, add both lines IDs to the failed IDs list 
		failed_synthon_IDs.append(in_synthon_class_ID)
		failed_synthon_IDs.append(ex_synthon_class_ID)
		return

######################################


# Load in synthon rules input file
synthon_rules_input_filename = sys.argv[1]
synthon_rules_input_file = open(synthon_rules_input_filename, 'r')
synthon_rules_input_lines = synthon_rules_input_file.readlines()
synthon_rules_input_file.close()

# Load in existing synthon rules JSON file
json_file_name = sys.argv[2]
json_file = open(json_file_name,)
Synthon_Map = json.load(json_file)
appended_json_file_name = sys.argv[3]


# Store Synthon IDs for all incorrect line in synthon input file
failed_synthon_IDs = []

# Synthon Input File in the form of:
# SynthonAClass__SynthonASubclass__ShouldContainAtLeastOne__- \t C<A> \t InclusionSMARTS_1 \t InclusionSMARTS_2 \t ...
# SynthonAClass__SynthonASubclass__shouldNotContain \t N<A> \t ExclusionSMARTS_1 \t ExclusionSMARTS_2 \t ...
# SynthonBClass__SynthonASubclass__ShouldContainAtLeastOne__- \t C<B> \t InclusionSMARTS_1 \t InclusionSMARTS_2 \t ...
# SynthonBClass__SynthonASubclass__shouldNotContain \t N<B> \t ExclusionSMARTS_1 \t ExclusionSMARTS_2 \t ...
# ...
# So loop over every other line, getting the line and line +1 to make each synthon class obj
synthon_rules = [extract_synthon_rules(synthon_rules_input_lines[i], synthon_rules_input_lines[i+1]) for i in range(0, len(synthon_rules_input_lines), 2) ]
# synthon_rules is list of Synthon_Chars objs

for synthon in synthon_rules: # For each synthon class, add a dictionary to the Synthon_Map dict. Get synthon info by Synthon_Map[ID]["whatever"]
	try:
		 Synthon_Map[synthon.SynthID] = {
			"Name" : synthon.SynthClass,
			"inSMARTS" : synthon.inSMARTS,
			"exSMARTS" : synthon.exSMARTS,
			}

	except AttributeError: # If a line in the input file is bad it will return a None in synthon_rules, giving an AttributeError None.SynthID
		continue

# If there were any bad input file lines, print those synthon IDs to an output file
if len(failed_synthon_IDs) != 0:
	errorfile_name = synthon_rules_input_filename + "_FAILED_SYNTHON_IDs"
	errorfile = open(errorfile_name, 'w')
	for ID in failed_synthon_IDs:
		errorstring = "%s\n" % (ID)
		errorfile.write(errorstring)
		
# Write out the appended Synthon_Map to a new JSON file
with open(appended_json_file_name, 'w') as outfile:
	json.dump(Synthon_Map, outfile)
outfile.close()

