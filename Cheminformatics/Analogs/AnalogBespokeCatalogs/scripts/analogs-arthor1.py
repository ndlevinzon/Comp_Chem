from ast import arg
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import urllib.parse
import requests
import os
from rdkit import DataStructs

# Imports from Khanh's script
import sys
import os
import arthor
import json
import pandas as pd



#.smi file of a BB
file=sys.argv[1] #'.sdf'

#BB SMARTS
# TO DO: pass as an argument
smarts="[NH2;+0][CX4;!$(C[OH,SH])][a!$(c1ccccc1)]"

def run_sub_query(string, limnum, db_path):
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
    bbs_dict = {}
    for i in rs:
        l_split = i.split()
        bbs_dict[l_split[1].decode('UTF-8')] = l_split[0].decode('UTF-8')
    print(bbs_dict)
    return(bbs_dict)
    # for i in rs:
    #         outfile.write("%s\n" % (i.decode('utf-8')))



def arthor_query(smarts):
    db="zinc22-all"
    query_string = smarts
    limit_no = 20000
    db_str = db
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
    bbs_dict=run_sub_query(query_string, limit_no, db_path)
    return(bbs_dict)

 

def read_arthor(filename):
    bbs_dict = {}
    with open(filename, 'r') as f:
        for line in f:
                # print(line)
            l_split = line.split()
            bbs_dict[l_split[1]] = l_split[0]

    print(bbs_dict)
    return(bbs_dict)


def calc_sim(bbs_dict, input, tc):
    bb_mol = Chem.MolFromSmiles(input)
    bb_fp = AllChem.GetMorganFingerprintAsBitVect(bb_mol,2,1024)
    bbs_dict_out = {}
    for keys in bbs_dict:
        mol = Chem.MolFromSmiles(bbs_dict[keys])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,1024)
        curr_tc = DataStructs.TanimotoSimilarity(fp, bb_fp)
        print(keys, curr_tc)
        if curr_tc >= tc:
            print(keys)
            bbs_dict_out[keys] = bbs_dict[keys]


        # bbs_dict_mol[keys] = [bbs_dict[keys], mol, AllChem.GetMorganFingerprintAsBitVect(mol,2,1024)]
    print("=================final BBs")
    print(bbs_dict_out)
    return(bbs_dict_out)

def write_out(bbs_dict_out,outfile):
    with open(outfile, 'w') as f:
        for keys in bbs_dict_out:
            write_str = bbs_dict_out[keys] +  " " + keys +  "\n"
            f.write(write_str)



def main(args):
    bbs_dict=arthor_query(args.smarts)
    # bbs_dict = read_arthor(args.arthor)
    bbs_dict_out = calc_sim(bbs_dict, args.input, float(args.tc))
    write_out(bbs_dict_out, args.output)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Query Arthor for a SMARTS and pick BBs similar to the defined BB with TC <= tc_cutoff",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input BB SMILES.")
    parser.add_argument("-s", "--smarts", type=str, help="Inclusion SMARTS.")
    parser.add_argument("-t", "--tc", type=str, help="Cutoff Tanimoto coefficient.")
    parser.add_argument("-o", "--output", type=str, help="Output files suffix name.")
#    parser.add_argument("-a", "--arthor", type=str, help="Arthor responce")
    args = parser.parse_args()
    main(args)


