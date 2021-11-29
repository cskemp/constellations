import networkx as nx
import matplotlib.pyplot as plt
import sys, os, glob
import json
sys.path.append(os.path.abspath(os.path.join('functions')))
import readgroups
import plotmodelsoln_single
import starhelpers as sh
import pandas as pd

# make plot for each culture
infiles  = ['../data/na_examples.txt']
#infiles  = ['../output/data/western_stellarium_processed.txt']


starfile = '../output/data/stars.txt'
brightthresh = 5.5
brightthresh = 5.92  # leopard with Zet2CrB
brightthresh = 6.11  # leopard with HR5827
stargraph = sh.makestargraph(starfile, [brightthresh])
stargraph = nx.create_empty_copy(stargraph)

invbayerd = sh.invbayerd(starfile)
invhd = sh.invhdict(starfile)

for infile in infiles:
    print(infile)
    #hc = readgroups.group(infile, invhd)
    #he = sh.mstedges(hc, stargraph)
    [hc, he] = readgroups.readbayer(infile, stargraph, invbayerd)

    prefix = os.path.splitext(os.path.basename(infile))[0]

    # make single plot with entire sky
    figfile = "../output/figures/single_" + prefix
    plotmodelsoln_single.plotmodelsoln(hc, he, stargraph, prefix ='', magthresh=brightthresh, edgecol=0, showlabs=0, show=1)
    plotmodelsoln_single.plotfullwithexamples(hc, he, stargraph, flip=0, prefix=figfile + '_egs', samecol=0,  magthresh=brightthresh, show=0)
    #plt.close()

