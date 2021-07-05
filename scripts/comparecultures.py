import json
import numpy as np
import functions.starhelpers as sh

starfile = '../output/data/stars.txt'
brightthresh = [3.5, 4.0, 4.5]
stargraph = sh.makestargraph(starfile, brightthresh)

allhumanfile = '../output/data/allhuman_bright_ourids.txt'

with open(allhumanfile, "rb") as hf:
    hall = json.load(hf)

scorefn = sh.scorecluster_wrtsystem
weights = [1,0,-1]
modelfile = '../output/results/model_n_F10.pickle'

cdiffs = [ [] for name in hall[0]]

hsys = [sh.llist2lset(csys) for csys in hall[1]]

allis = np.arange(len(hall[0]))
for i, cname in enumerate(hall[0]):
    print(cname)
    minusi= np.delete(allis, i)
    sysi = [hsys[i]]
    rest = [hsys[mi] for mi in minusi]
    cdiffs[i] = sh.compare_humantomodelsys(rest, sysi, scorefn, weights)[0]

outfile = '../output/results/culturesvsotherscores.csv'
sh.printperculturescores(outfile, hall[0], cdiffs)

