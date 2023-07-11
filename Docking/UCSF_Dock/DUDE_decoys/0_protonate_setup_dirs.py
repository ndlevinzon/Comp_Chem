import os, sys

def protonate_ligands(path, prot_dict):

        prot_dir = path+"test_ligand_protonation/"
        if not os.path.isdir(prot_dir):
                os.system("mkdir "+prot_dir)
        os.chdir(prot_dir)

        to_prot_file = prot_dir+"test_ligand_ids.smi"
        output = open(to_prot_file, 'w')
        for zinc_id in prot_dict:
                smiles = prot_dict[zinc_id][0]
                output.write(smiles+" "+zinc_id+"\n")
        output.close()

        os.system("source ~tbalius/.cshrc_dbgen_corina")
        os.system("bash /mnt/nfs/export/rstein/DUDE_Z/new_DUDE_pipeline_4_2_2018/generate_protomers_from_smiles.sh -H 7.4 "+to_prot_file)

        built_dir = prot_dir+"/working/protonate/"
        built_file = built_dir+"test_ligand_ids-protomers.ism"

        new_prot_dict = {}
	repeat_list = []
        if os.path.isfile(built_file):
                open_built = open(built_file, 'r')
                read_built = open_built.readlines()
                open_built.close()

                for line in read_built:
                        splitline = line.strip().split()
                        smiles = splitline[0]
                        lig_ID = splitline[-2]
			lig_ID_count = repeat_list.count(lig_ID)
			dict_entry = lig_ID+"_"+str(lig_ID_count)
			repeat_list.append(lig_ID)
			
                        if dict_entry not in new_prot_dict:
                                new_prot_dict[dict_entry] = [smiles]
                        else:
                                new_prot_dict[dict_entry].append(smiles)

        os.chdir(path)
        os.system("rm -rf "+prot_dir)

        return(new_prot_dict)

def gather_ligands(ligand_file):

        open_lig = open(ligand_file,'r')
        read_lig = open_lig.readlines()
        open_lig.close()

        lig_dict = {}
        for line in read_lig:
                splitline = line.strip().split()
                lig_ID = splitline[-1]
                smiles = splitline[0]
                lig_dict[lig_ID] = [smiles]

        return(lig_dict)

def create_lig_dirs(smiles_dir, prot_dict, infile):

        lig_count = 0
        for prot in prot_dict:
                smiles = prot_dict[prot][0]
                lig_count += 1

                print(smiles, prot)
                lig_dir = smiles_dir+"/ligand_"+str(lig_count)+"/"
                if not os.path.isdir(lig_dir):
                        os.system("mkdir "+lig_dir)

                        os.chdir(lig_dir)
			os.system("cp "+infile+" .")
                        lig_set_name = "lig_set_"+str(prot)+".smi"
                        output = open(lig_set_name, 'w')
                        output.write(smiles+" "+prot+"\n")
                        output.close()

                        os.system("cp "+lig_set_name+" ligand_"+str(lig_count)+".smi")

                        os.chdir(smiles_dir)

def read_infile(infile):

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()

	for line in read_in:
		splitline = line.strip().split()
		if len(splitline) > 0:
			if splitline[0].lower() == "protonate":
				tf = splitline[-1]

	return(tf)

def main():

	
	pwd = os.getcwd()+"/"
	if len(sys.argv) != 3:
		print("Syntax: python 0000_protonate_setup_dirs.py smiles_file.smi directory_name")
		sys.exit()

	infile = pwd+"decoy_generation.in"
	protonate_ligs = False
	if not os.path.isfile(infile):
		print("decoy_generation.in does not exist")
		sys.exit()
	else:
		tf = read_infile(infile)
		if tf.lower()[0] == "y":
			protonate_ligs = True
		
	smiles_file = pwd+sys.argv[1]
	smiles_dir = pwd+sys.argv[2]+"/"

	if not os.path.isdir(smiles_dir):
		os.system("mkdir "+smiles_dir)

	os.chdir(smiles_dir)
	os.system("cp "+infile+" .")
		
	lig_dict = gather_ligands(smiles_file)
	if protonate_ligs == True:
		prot_dict = protonate_ligands(smiles_dir, lig_dict) 
		create_lig_dirs(smiles_dir, prot_dict, infile)
	else:
		create_lig_dirs(smiles_dir, lig_dict, infile)


main()
