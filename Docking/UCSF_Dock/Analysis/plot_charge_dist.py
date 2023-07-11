#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 13:46:00 2019

@author: stefan
"""

import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp
from optparse import OptionParser
import math

def plot_charges(datalist, output):
    plt.figure()
    plt.xlabel('Charge')
    plt.ylabel('Energy in DOCK units')
    plt.xlim([-2.5,2.5])
    plt.ylim([-70,-30])
    for x in range(len(datalist)):
        data=pd.read_csv(datalist[x], header=None)
        if x==0:
            plt.plot(data[1], data[2], 'x', markersize=3.5)
        elif x%2==0:
            plt.plot(data[1]+0.075*(math.ceil(x/2.0)), data[2], 'x', markersize=3.5)
        elif x%2==1:
            plt.plot(data[1]-0.075*math.ceil(x/2.0), data[2], 'x', markersize=3.5)
            
    outputname=output+".png"
    plt.legend(datalist, loc='lower right', fontsize='small')
    plt.title("Charge distribution")
    plt.savefig(outputname, dpi=600)

def main(argv):
    description = "Plot charge distributions"
    usage = "%prog [options]"
    parser = OptionParser(usage=usage, description=description)
    parser.set_defaults(inputfiles=[], output="Charge-distributions")
    parser.add_option("-i", "--inputfiles", action="append", help="Input. You can input multiple files")
    parser.add_option("-o", "--output", help="Output")
    
    options, args = parser.parse_args(args=argv[1:])
    
    plot_charges(datalist=options.inputfiles, output=options.output)
    return 0

if __name__ == "__main__":
    plt.matplotlib.use('Agg')
    sys.exit(main(sys.argv))
