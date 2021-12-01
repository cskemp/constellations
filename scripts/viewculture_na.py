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
infiles  = ['../data/na_councilofchiefs.txt']
infiles  = ['../data/na_bearden.txt']
infiles  = ['../data/na_kaitusiyuman.txt']

infiles  = ['../data/na_councilofchiefs.txt',
            '../data/na_bearden.txt',
            '../data/na_kaitusiyuman.txt']

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
    thisstargraph = nx.create_empty_copy(stargraph)
    [hc, he] = readgroups.readbayer(infile, thisstargraph, invbayerd)

    # if we want to plot all stars
    # thisstargraph = nx.create_empty_copy(stargraph)
    # to plot stars in first asterism only
    thisstargraph = nx.create_empty_copy(stargraph.subgraph(hc[1]))

    prefix = os.path.splitext(os.path.basename(infile))[0]

    # make single plot with entire sky
    figfile = "../output/figures/single_" + prefix
    plotmodelsoln_single.plotmodelsoln(hc, he, thisstargraph, prefix ='', magthresh=brightthresh, edgecol=0, showlabs=0, show=0)
    plotmodelsoln_single.plotna(hc, he, thisstargraph, flip=0, prefix=figfile + '_egs', samecol=0,  magthresh=brightthresh, show=1)
    #plt.close()

