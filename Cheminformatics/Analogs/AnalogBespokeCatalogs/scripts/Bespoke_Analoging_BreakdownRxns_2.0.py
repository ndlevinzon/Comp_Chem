# Written by Seth Vigneron, Shoichet Lab UCSF, 2022
# Takes an input file of molecules and an optional portions of those compounds to preserve, and returns reactions that could be used for analoging

import json
import rdkit
from rdkit import Chem as C
from rdkit.Chem import AllChem
import sys

#####################################################################################################################
#####################################################################################################################

# Class Objects for Base Molecules to Analog
#####################################################################################################################
class BaseMolecule:
	
	def __init__(self, smiles, cid, mol):
		self.smi = smiles
		self.cid = cid
		self.mol = mol


#######################################
class Reaction:
	
	def __init__(self, Name, ID, Synthons, breakdownSMARTS, breakdownRxn, buildSMARTS, buildRxn, numReactants):
		self.Name = Name
		self.ID = ID
		self.Synthons = Synthons
		self.breakdownSMARTS = breakdownSMARTS
		self.breakdownRxn = breakdownRxn
		self.buildSMARTS = buildSMARTS
		self.buildRxn = buildRxn
		self.numReactants = numReactants



#####################################################################################################################
#####################################################################################################################

# Functions
#####################################################################################################################

def Extract_Rxns(rxnIDs_list):
	rxns_list = []
	for rxn_ID in rxnIDs_list:
		name = Analoging_Rxn_Map["Reactions"][rxn_ID]["Name"]
		synthonIDs = Analoging_Rxn_Map["Reactions"][rxn_ID]["SynthonIDs"]
		synthonIDs_list = synthonIDs.split(',')
		breakdown_SMARTS = Analoging_Rxn_Map["Reactions"][rxn_ID]["Breakdown_SMARTS"]
		breakdown_rxn = AllChem.ReactionFromSmarts(breakdown_SMARTS)
		build_SMARTS = Analoging_Rxn_Map["Reactions"][rxn_ID]["Build_SMARTS"]
		build_rxn = AllChem.ReactionFromSmarts(build_SMARTS)
		numReactants = Analoging_Rxn_Map["Reactions"][rxn_ID]["NumReactants"]
	
		rxns_list.append(Reaction(name, rxn_ID, synthonIDs_list, breakdown_SMARTS, breakdown_rxn, build_SMARTS, build_rxn, numReactants))
		
	return rxns_list


# Generate BaseMolecule Class objects for each input compound to be analoged
def Extract_Base_Mols(file_lines):
	mol_list = []
	for line in file_lines:
		stripline = line.strip().split()
		smiles = stripline[0]
		cid = stripline[1]
		mol = C.MolFromSmiles(smiles)
		
		mol_list.append(BaseMolecule(smiles, cid, mol))

	return mol_list
# Check to to see if all the mols of the BBs are big enough to be considered an actual BB (ie not just NH3)
def check_bbs_HAC_min(bb_mols_list):
        not_too_small = True
        for bb_mol in bb_mols_list:
                if bb_mol.GetNumHeavyAtoms() <= 3:
                        not_too_small = False
        return not_too_small
#######################################
# Check if a molecule can be broken down each of the analoging reactions and return the synthons
def Run_Analgoing_Assess(base_mol):
	analog_dict = {} # Start a dict to add rxnID:[[bb1, synthonID_bb1], [bb2, synthonID_bb2]] for each rxn
	utilized_rxns_list = []
	for analog_rxn in Analoging_Rxns: # Loop through all the reactions
		rxn = analog_rxn.breakdownRxn
		num_reactants = int(analog_rxn.numReactants)
		#try: # Try using each of the analoging breakdown reactions, if they are successfull, label the base molecule with that reaction ID
		BB_products = rxn.RunReactants((base_mol.mol,))
		try:
			BBs = BB_products[0]
			if len(BBs) == num_reactants and check_bbs_HAC_min(BBs):
				rxnID = analog_rxn.ID
				rxnSynthonIDs = analog_rxn.Synthons
				BB_mols = []
				for i in range(0, num_reactants):
					#print(BB_mols)
					#print(num_reactants)
					BB_out = [BBs[i], rxnSynthonIDs[i]]	
					BB_mols.append(BB_out)
				
				analog_dict[rxnID] = BB_mols # Potentially have an if statement that checks to see if the BBs pass the rxn rules before adding it to the dict
				utilized_rxns_list.append(analog_rxn)
		except IndexError:
			continue

	return analog_dict, utilized_rxns_list
		




#####################################################################################################################
#####################################################################################################################

# Main
#####################################################################################################################


# load input files. Input file for compounds to analog should be in the form: cid smiles smiles_to_preserve
analoging_reactions_json = sys.argv[1]
base_molecules_file_name = sys.argv[2]
base_molecules_file = open(base_molecules_file_name, 'r')
base_molecules_file_lines = base_molecules_file.readlines()
base_molecules_file.close()


# Read in json Analoging_Rxn_Map_json
json_file = open(analoging_reactions_json,)
Analoging_Rxn_Map = json.load(json_file)
Analoging_Rxn_IDs = Analoging_Rxn_Map["Reaction_IDs"]
# Json map should contain 2 keys:
#       Analoging_Rxn_Map["Reaction_IDs"] gives all rxn IDs of analoging reactons stored in the Json map
#		Analoging_Rxn_Map["Reactions"] gives a dictionary where the keys are the rxn IDs

# Create a list of Reaction Class objects for all analoging reactions in the Analoging Rxn Map
Analoging_Rxns = Extract_Rxns(Analoging_Rxn_IDs)

# Get all the mol objects for the compounds to analog
base_mol_obj_list = Extract_Base_Mols(base_molecules_file_lines)

for molecule in base_mol_obj_list:
	analog_dict, utilized_rxns = Run_Analgoing_Assess(molecule)

	outfile_name = '%s_analoging' % (base_molecules_file_name)
	outfile = open(outfile_name, 'w')

	# Write out resulsts of rxns that can be used to analog the given base compound
	# BBs of the broken down base compound listed in form BB1_smi,BB2_smi,...,BBn_smi \t BB1_synthonID,BB2_synthonID,...,BBn_synthonID
	for rxn in utilized_rxns:
		bb_list = analog_dict[rxn.ID]
		bb_outstring = ''
		bb_smiles_list = []
		bb_synthon_ID_list = []
		for bb in bb_list:
			bb_smi = C.MolToSmiles(bb[0])
			bb_synthon_ID = bb[1]
			outstring = '%s\t%s\t%s\t%s\t%s\n' % (molecule.smi, rxn.ID, rxn.Name, bb_smi, bb_synthon_ID)
			outfile.write(outstring)

	outfile.close()
	
		
	

