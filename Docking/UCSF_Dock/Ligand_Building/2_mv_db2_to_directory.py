import sys, os

#This file will take the output of the dbgen script and will copy all the db2 files into a common directory 


def main():
	if len(sys.argv) != 2:
		print('This script needs: ')
		print('1: Directory path to copy all db2 files to')
		
		sys.exit()
	
	pwd = os.getcwd() + '/'
	copy_dir = sys.argv[1]
	current_dir = [name for name in os.listdir(pwd)]# if name.startswith('ZINC')]
	print(current_dir)
	
	for d in current_dir:
		print(d)

		os.chdir(d)
		print('Working on this mol: ')
		npwd = os.getcwd() + '/'
		
		
		rename_list = [name for name in os.listdir(npwd) if name.endswith('db2.gz')]
		print(rename_list)
		
	#	zinc_id = [name for name in os.listdir(npwd) if name.endswith('.solv')]# and name.startswith('ZINC')]
	#	z = zinc_id[0][:-5]
		print('Renaming protomers for this mol: ' + d)
		
		
		for n in rename_list:
			new_name = d + '_'+  n
			print('Moving ' + new_name)
		

			os.system('cp ' + n + ' ' + new_name)
		
			os.system('cp ' + new_name + ' ' + copy_dir)
		os.chdir(pwd)
main()
