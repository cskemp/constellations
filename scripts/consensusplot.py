import functions.starhelpers as sh
import functions.readgroups as readgroups
import functions.plotmodelsoln_single as plotmodelsoln_single
import json
import numpy as np
import networkx as nx
import itertools

starfile = '../output/data/stars.txt'
brightthresh = [3.5, 4.0, 4.5]
stargraph = sh.makestargraph(starfile, brightthresh)
invhd = sh.invhdict(starfile)

with open("../output/data/hrnames.txt", "r") as f:
    hrnames = dict(json.load(f))

# human data
hdfile_bright = '../output/data/allhuman_bright_ourids.txt'
with open(hdfile_bright, 'r') as f:
    hd_bright = json.load(f)
clist_bright = [[]] + list(itertools.chain.from_iterable(hd_bright[1]))
hc = sh.llist2lset(clist_bright)

he = sh.mstedges(hc, stargraph)

sh.edge_increment(stargraph, he)

edgefile = '../output/results/allhumanedges_' + str(brightthresh) + '.txt'
sh.printedgecounts(stargraph, [], edgefile)

figfile = "../output/figures/consensus_edgethick"
magthresh = np.max(brightthresh)

# drop edges that appear fewer than 4 times

thresh_sg = sh.dropsmalledges(stargraph, 'ecount', 4, alledge=1)
cs = list(nx.connected_components(thresh_sg))
cs = [ c for c in cs if len(c) > 1]
cs = [[]] + cs
es = sh.allmodeledges(cs, thresh_sg)

plotmodelsoln_single.plotmodelsoln(cs, es, thresh_sg, show = 0, prefix='', samecol=0, edgecount=1, magthresh = magthresh)
plotmodelsoln_single.plotfullwithexamples(cs, es, thresh_sg, flip = 0, prefix=figfile+'_egs', samecol=0, edgecount=1, magthresh = magthresh)
