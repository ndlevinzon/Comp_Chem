import os,sys,traceback
import subprocess

## Script for retrieving poses after docking
## Written by Reed Stein with assistance from T.E. Balius, 2019
## Modified by Reed Stein to make 3x faster, Dec, 2019

def get_zinc_names_and_chunks_from_extract_all_sort_uniq(filename,N):
	print("running function: get_zinc_names_and_chunks_from_extract_all_sort_uniq")
	fh = open(filename,'r')
	lig_dict = {}
	chunk_dict = {}
	count = 0
	for line in fh:
		splitline = line.split()
		zincid = splitline[2]
		chunk  = splitline[0]
		outindex = splitline[22] # 23rd element contains the output index
		if not (chunk in chunk_dict):
			chunk_dict[chunk] = [outindex]
		else:
			if not (outindex in chunk_dict[chunk]):
				chunk_dict[chunk].append(outindex)
		lig_dict[zincid] = [500] ## make large energy to start with
		#print chunk, zincid 
		count = count + 1
		if count > N: 
			break
	return(lig_dict, chunk_dict)

def process_gz(filename, lig_dict, chunk):
	print("running function: process_gz")
        #this function collects the name, grids, and total energy of the best scoring pose of each ligand

	lig_list = []
	ligandFile = subprocess.Popen(['gzip', '-cdfq', filename], stdout=subprocess.PIPE)

	better_energy = False
	comp_found = False
	for line in ligandFile.stdout:
		splt = line.strip().split()
		lig_list.append(line)

		if len(splt) > 2:
			if splt[0] == "##########" and splt[1] == "Name:" and better_energy == True:
		#		print(name, splt[2], comp_found, better_energy, temp_tote, new_tote)
				del lig_dict[name]
				lig_dict[name] = [new_tote, lig_list[0:-1]]
				lig_list = []
				better_energy = False
				comp_found = False

			if splt[0] == "##########" and splt[1] == "Name:" and better_energy == False:
				name = splt[2]
				if name in lig_dict:
					comp_found = True
					#temp_tote = lig_dict[name][0]
					lig_list = []
					lig_list = [line]
				else:
					comp_found = False # make sure to redefine
					lig_list = []

			if splt[0] == "##########" and splt[1] == "Total" and splt[2] == "Energy:" and comp_found == True:
				temp_tote = lig_dict[name][0]
				new_tote = float(splt[3])
				if new_tote < temp_tote:
					better_energy = True
				#else:
				#	comp_found = False

	#### only when molecule is the last in the test.mol2.gz file
	if comp_found == True and better_energy == True:
		del lig_dict[name]
		lig_dict[name] = [new_tote, lig_list]
		
	ligandFile.terminate()
	
	return(lig_dict)


def main():

	if (len(sys.argv) != 6):
		print("Give this script 4 inputs:")
		print("(1) path to where docking is located. uses '' if not needed. ")
		print("(2) path to where the extract all file is. ")
		print("(3) number of molecules (poses) to get. ")
		print("(4) file name for poses eg. 'poses.mol2'")
		print("(5) file name for input eg. 'test.mol2.gz'")
		exit()

	docking_dir = sys.argv[1]

        #extractname = sys.argv[2]
	extractfile = sys.argv[2]
	number_of_poses = int(sys.argv[3])
	posesfilename   = sys.argv[4]
	gz_file_name    = sys.argv[5]
	print("docking_dir: "+docking_dir)
	#print "extractname: "+extractname
	#extractfile = docking_dir+extractname
	print("extract file path: "+extractfile)
	print("number_of_poses: "+str(number_of_poses))        
	#os.chdir(lig_dir)
	
	#if os.path.isfile(lig_dir+"poses.mol2"):
	#if os.path.isfile(docking_dir+"poses.mol2"):
	#if os.path.isfile("poses.mol2"):
	if os.path.isfile(posesfilename):
		print(posesfilename+" already exists. Quitting.")
		sys.exit()
	
	#extractfile = docking_dir+"extract_all.sort.uniq.txt"
	if not os.path.isfile(extractfile):
		print("there needs to be an extract_all.sort.uniq.txt. ")
		exit()

	lig_dict, chunk_dict = get_zinc_names_and_chunks_from_extract_all_sort_uniq(extractfile,number_of_poses)


	for chunk in chunk_dict:
		#print(chunk) 
		#gz_file = docking_dir+chunk+"/"+gz_file_name
		#lig_dict = process_gz(gz_file, lig_dict, chunk)
		try:
			for outindex in chunk_dict[chunk]:
				#print chunk, chunk_dict[chunk]
                	        # jji adding element 23 here
				gz_file = docking_dir+chunk+"/"+gz_file_name + outindex
				print gz_file
				lig_dict = process_gz(gz_file, lig_dict, chunk)
		except:
			traceback.print_exc()
			print chunk+" is broken"


	output = open(posesfilename, 'w')
	for l in lig_dict:
		if len(lig_dict[l]) < 2:
			print(l+" may be above your energy cutoff. This is typically -10 kcal/mol. Check your OUTDOCK/extract_all.sort.uniq.txt file to confirm.")
		else:
			for i in range(len(lig_dict[l][1])):
				output.write(lig_dict[l][1][i])

	output.close()

	#write_out_poses(lig_dict, docking_dir, posesfilename)
	
main()

