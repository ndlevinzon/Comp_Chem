import sys
import os
import arthor
import json
import pandas as pd


def run_sub_query(string, limnum, db_path, outfile):

	db = arthor.SubDb('{}.smi.atdb'.format(db_path), smifile='{}.smi'.format(db_path))
	try:
		qry = arthor.Query(string, "Sma:QI")
	except Exception as err:
		qry = arthor.Query(string, "Sma:QI")
	if limnum is None:
		rs = db.search(qry)
	else:
		rs = db.search(qry, limit=limnum)
	print(rs)
	#with open("outbytes", 'wb') as f:
	#	f.write(rs)
	print(len(rs))
	#print(rs.decode('UTF-8'))
	for i in rs:
		outfile.write("%s\n" % (i.decode('utf-8')))
	
if __name__ == '__main__':
	infile = open(sys.argv[1], "r")
	lines = infile.readlines()
	query_string = lines[0].strip('\n')
	limit_no=None
	if (sys.argv[2]) != 'All':
		limit_no = int(sys.argv[2])
	dest_dir = os.getcwd()
	db_str = sys.argv[3]
	db_file = open('/nfs/home/khtang/work/Projects/arthor_batch_search/arthor_databases.txt', 'r')
	db_list = db_file.readlines()
	db_path = ''
	for line in db_list:
		line.strip('\n')
		x = line.split(':')
		if db_str == x[0]:
			db_path = x[1].strip('\n')
	if len(sys.argv) == 5:
		dest_dir = sys.argv[4] 
	#for line in lines:
	#	smiles, cid = line.strip('\n').split(' ')
	#	print(smiles, cid)
	outfile = open("{}/result_query.txt".format(dest_dir), 'w')
	run_sub_query(query_string, limit_no, db_path, outfile)
