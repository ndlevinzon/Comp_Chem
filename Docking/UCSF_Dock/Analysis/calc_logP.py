from __future__ import print_function, absolute_import
import os, sys
from rdkit import Chem as C
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD


def get_stuff(smiles):

        mol = C.MolFromSmiles(smiles)
        hac = D.HeavyAtomCount(mol)
        if len(str(mol)) == 0:
                return(0, 0, 0, 0, 0, 0)

        else:
                mw = CD.CalcExactMolWt(mol)
                logp = C.Crippen.MolLogP(mol)
                #rotB = D.NumRotatableBonds(mol)
                #HBA = CD.CalcNumHBA(mol)
                #HBD = CD.CalcNumHBD(mol)
                #q = C.GetFormalCharge(mol)

                return(mw, logp, hac)


def main():

	smi_file = open(sys.argv[1], 'r')
	read_smi = smi_file.readlines()
	smi_file.close()

	output = open("mw_logp.txt", 'w')
	for line in read_smi:
		splt = line.strip().split()
		smiles = splt[0]
		zinc = splt[1]
		mw, logp, hac = get_stuff(smiles)
		output.write(smiles+"\t"+zinc+"\t"+str(mw)+"\t"+str(logp)+"\n")

	output.close()
			


main()
