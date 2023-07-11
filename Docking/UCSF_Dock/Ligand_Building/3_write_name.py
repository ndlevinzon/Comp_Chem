import os, sys


def main():

	dec_dict = {}
	dec_list = [name for name in os.listdir(".") if (os.path.isfile(name) and name.endswith(".db2.gz"))]
	
	output = open("decoys.name", 'w')
	for d in dec_list:
		if "_" in d:
			name = d.split("_")[0]
			output.write(name+"\n")
		else:
			name = d.split(".")[0]
			output.write(name+"\n")
	output.close()

main()
