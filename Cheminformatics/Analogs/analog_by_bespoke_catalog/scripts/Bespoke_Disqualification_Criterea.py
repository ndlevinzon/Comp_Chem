# Written by Seth Vigneron, Shoichet Lab UCSF, 2021 

# Define disqualification functions that take return true if the product compound is disqualified base on the cut off value provided in the json input file
# Note: The name of the function must be the same as the name of the restriction type in the json file + '_check'
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def Max_Product_HAC_check(self):
	mol_HAC = self.mol.GetNumHeavyAtoms()
	if mol_HAC > int(self.Max_Product_HAC):
		return True

def Min_Product_HAC_check(self):
	mol_HAC = self.mol.GetNumHeavyAtoms()
	if mol_HAC < int(self.Min_Product_HAC):
		return True

def Max_Product_LogP_check(self):
	mol_logP = Chem.Crippen.MolLogP(Chem.RemoveHs(self.mol))
	if float(mol_logP) > float(self.Max_Product_LogP):
		return True
	
def Max_Product_RotatableBonds_check(self):
	mol_rotbonds = rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(self.mol))
	if mol_rotbonds > int(self.Max_Product_RotatableBonds):
		return True

def Min_Product_RotatableBonds_check(self):
	mol_rotbonds = rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(self.mol))
	if mol_rotbonds < int(self.Min_Product_RotatableBonds):
		return True

def Max_Purchasability_check(self):
	purch_list = []
	for bb in self.BB_combo:
		purch_list.append(int(bb.purch))
	max_purch = 50 * len(purch_list)
	prod_purch_count = sum(purch_list)
	purch_percent = prod_purch_count / max_purch
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

	if product_purchasability > int(self.Max_Purchasability):
		return True

def Max_MolWeight_check(self):
	mol_MolWeight = Descriptors.MolWt(Chem.RemoveHs(self.mol))
	if mol_MolWeight > float(self.Max_MolWeight):
		return True

def Min_MolWeight_check(self):
	mol_MolWeight = Descriptors.MolWt(Chem.RemoveHs(self.mol))
	if mol_MolWeight < float(self.Min_MolWeight):
		return True

def Max_HBond_Donors_check(self):
	mol_HBondDonors = Chem.Lipinski.NumHDonors(Chem.RemoveHs(self.mol))
	if mol_HBondDonors > int(self.Max_HBond_Donors):
		return True

def Min_HBond_Donors_check(self):
	mol_HBondDonors = Chem.Lipinski.NumHDonors(Chem.RemoveHs(self.mol))
	if mol_HBondDonors < int(self.Min_HBond_Donors):
		return True

def Max_HBond_Acceptors_check(self):
	mol_HBondAcceptors = Chem.Lipinski.NumHAcceptors(Chem.RemoveHs(self.mol))
	if mol_HBondAcceptors > int(self.Max_HBond_Acceptors):
		return True

def Min_HBond_Acceptors_check(self):
	mol_HBondAcceptors = Chem.Lipinski.NumHAcceptors(Chem.RemoveHs(self.mol))
	if mol_HBondAcceptors < int(self.Min_HBond_Acceptors):
		return True

def Max_Unspecified_Stereocenters_check(self):
	mol_UnspecStereocenters = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(Chem.RemoveHs(self.mol))
	if mol_UnspecStereocenters > int(self.Max_Unspecified_Stereocenters):
		return True

def Max_TPSA_check(self):
	mol_TPSA = rdMolDescriptors.CalcTPSA(Chem.RemoveHs(self.mol))
	if mol_TPSA > float(self.Max_TPSA):
		return True

def Min_TPSA_check(self):
	mol_TPSA = rdMolDescriptors.CalcTPSA(Chem.RemoveHs(self.mol))
	if mol_TPSA < float(self.Min_TPSA):
		return True

def Max_Fraction_sp3_C_check(self):
	mol_FracSp3C = Chem.Lipinski.FractionCSP3(Chem.RemoveHs(self.mol))
	if mol_FracSp3C > float(self.Max_Fraction_sp3_C):
		return True

def Min_Fraction_sp3_C_check(self):
	mol_FracSp3C = Chem.Lipinski.FractionCSP3(Chem.RemoveHs(self.mol))
	if mol_FracSp3C < float(self.Min_Fraction_sp3_C):
		return True

def Max_NHOH_Count_check(self):
	mol_NHOH_Count = Chem.Lipinski.NHOHCount(Chem.RemoveHs(self.mol))
	if mol_NHOH_Count > int(self.Max_NHOH_Count):
		return True

def Min_NHOH_Count_check(self):
	mol_NHOH_Count = Chem.Lipinski.NHOHCount(Chem.RemoveHs(self.mol))
	if mol_NHOH_Count < int(self.Min_NHOH_Count):
		return True

def Max_NO_Count_check(self):
	mol_NO_Count = Chem.Lipinski.NOCount(Chem.RemoveHs(self.mol))
	if mol_NO_Count > int(self.Max_NO_Count):
		return True

def Min_NO_Count_check(self):
	mol_NO_Count = Chem.Lipinski.NOCount(Chem.RemoveHs(self.mol))
	if mol_NO_Count < int(self.Min_NO_Count):
		return True

def Max_Alkyl_Halide_check(self):
	mol_alkyl_halides = Chem.Fragments.fr_alkyl_halide(Chem.RemoveHs(self.mol))
	if mol_alkyl_halides > int(self.Max_Alkyl_Halide):
		return True

def Min_Alkyl_Halide_check(self):
	mol_alkyl_halides = Chem.Fragments.fr_alkyl_halide(Chem.RemoveHs(self.mol))
	if mol_alkyl_halides < int(self.Min_Alkyl_Halide):
		return True

# Add more disqualification functions here
