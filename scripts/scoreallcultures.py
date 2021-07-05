import json
import numpy as np
import functions.starhelpers as sh


starfile = '../output/data/stars.txt'
brightthresh = [3.5, 4.0, 4.5]
stargraph = sh.makestargraph(starfile, brightthresh)

allhumanfile = '../output/data/allhuman_bright_ourids.txt'

with open(allhumanfile, "rb") as hf:
    hall = json.load(hf)

magthresh = np.max(brightthresh)

cscores = [ [] for name in hall[0]]
scorefn = sh.scorecluster_wrtsystem
weights = [1,0,-1]
modelfile = '../output/results/model_n_F10.pickle'

hsys = [sh.llist2lset(csys) for csys in hall[1]]
mscores = sh.consensus_modelscores(modelfile, hsys, scorefn, weights)

outfile = '../output/results/perculturescores.csv'
sh.printperculturescores(outfile, hall[0], mscores)

