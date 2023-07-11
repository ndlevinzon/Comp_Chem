#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 17:18:30 2019

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
import pylab


def plot_roc(datalist, output):
    #set equal aspect
    plt.figure().set_size_inches(8., 8.)
    plt.rc('font', size=15)
    plt.xlabel('% Decoys docked')
    plt.ylabel('% Ligands docked')
    plt.ylim([0,100])
    plt.xlim([0.1,100])
    seq=np.linspace(0,100, 1000)
    plt.plot(seq, seq, 'k--')
    for i in range(len(datalist)):
        data=pd.read_csv(datalist[i], header=None, sep='\t') #read in data from roc_own.txt
        #Transform data into vectors for plotting
        x=data[1:][0]
        xdat=list(map(float, list(x)))
        xvec=np.array(xdat)
        y=data[1:][1]
        yvec=np.array(list(y))
        #plot % decoys docked on logarithmic scale
        plt.semilogx(xvec, yvec, label=datalist[i]+": "+str(data[3][0]))
    
    outputname=output+".png"
    plt.legend(loc='best')
    plt.title("LogAUC")
    plt.savefig(outputname, dpi=600)


def main(argv):
    description = "Plot Multiple LogAUC curves"
    usage = "%prog [options]"
    parser = OptionParser(usage=usage, description=description)
    parser.set_defaults(inputfiles=[], output="LogAUC")
    parser.add_option("-i", "--inputfiles", action="append", help="Input (roc_own.txt). You can input multiple files")
    parser.add_option("-o", "--output", help="Output")
    
    options, args = parser.parse_args(args=argv[1:])
    
    plot_roc(datalist=options.inputfiles, output=options.output)
    return 0

if __name__ == "__main__":
    plt.matplotlib.use('Agg')
    sys.exit(main(sys.argv))
