import sys, os


def main():
	pwd = os.getcwd() + '/'
	
	ligands = [name for name in os.listdir(pwd) if os.path.isdir(pwd + name)]
	sdi_file = open('ligand_database_index', 'a')
	
	for l in ligands:
		
		protomers = [x for x in os.listdir(pwd + l) if x.endswith('db2.gz')]
		for p in protomers:

			sdi_file.write(pwd + l + '/'+ p + '\n')

	sdi_file.close()
	







main()
