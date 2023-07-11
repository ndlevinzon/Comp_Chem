#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:23:23 2020

@author: stefan
"""

import sys, os
import argparse
import time


start=time.time()
def main():

#pwd = os.getcwd()+"/"

    parser = argparse.ArgumentParser(description='''This script calculates the RMSD between docked poses and reference structrues (e.g. xtal-lig). You need to generate reference files first. Remove hydrogen atoms (openbabel) and save as NAME-ref.mol2''')
	
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

    parser.add_argument('-i', '--input', type=str, dest="inlist", required=True, help="List of molecule names to check")
    parser.add_argument('-p', '--poses', type=str, dest="poses", required=True, help="Full path to poses.mol2")
    parser.add_argument('-o', '--output', type=str, dest="outlist", default='RMSD-results.txt', help="Output file")

    args = parser.parse_args()

    open_in = open(args.inlist, 'r')
    read_in = open_in.readlines()
    open_in.close()

    output=open(args.outlist, 'w')
 #   in_dict={}
    for line in read_in:
        splt=line.strip().split()
        ID=splt[0]
        os.system("echo "+ID+" > ID")
        print("python lc_blazing_fast_collect_mol2.py ID "+args.poses+" "+ID+".mol2")
        print("/nfs/soft/openbabel/current/bin/obabel -imol2 "+ID+".mol2"+" -d -omol2 -O "+ID+"-noH.mol2")
        print("csh calculate_lig_rmsd.csh "+ID+"-noH.mol2 "+ID+"-ref.mol2 > RMSD.dat")
        print("grep HA_RMSDh RMSD.dat | awk '{print $3}' > val.dat")
        val=float(open("val.dat", 'r').readlines()[0].split()[0])
        output.write(ID+" "+val+"\n")
        print("rm RMSD.dat val.dat")
           
if __name__ == '__main__':
	main()

end=time.time()
print(end - start)
