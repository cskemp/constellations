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
infiles  = list(glob.iglob(os.path.abspath(os.path.join('..', 'output/data/*_processed.txt'))))

starfile = '../output/data/stars.txt'
brightthresh = 4.5
stargraph = sh.makestargraph(starfile, [3.5, 4.0, 4.5])

invhd = sh.invhdict(starfile)

conncnts = []
conncntsfile = '../output/results/conncounts.csv'

for infile in infiles:
    print(infile)
    hc = readgroups.group(infile, invhd)
    he = sh.mstedges(hc, stargraph)

    bc = sh.brightclusters(hc, stargraph, brightthresh)
    ccnt, ncnt = sh.printclusters(bc, stargraph)
    conncnts.append( (os.path.basename(infile), ccnt, ncnt))

    prefix = os.path.splitext(os.path.basename(infile))[0]

    # make single plot with entire sky
    figfile = "../output/figures/single_" + prefix
    plotmodelsoln_single.plotmodelsoln(hc, he, stargraph, prefix ='', magthresh=brightthresh, edgecol=0, showlabs=0, show=0)
    plotmodelsoln_single.plotfullwithexamples(hc, he, stargraph, flip=0, prefix=figfile + '_egs', samecol=0,  magthresh=brightthresh, show=0)
    plt.close()

d = pd.DataFrame(conncnts)
d.to_csv (conncntsfile, index=False)
