import json
import functions.starhelpers as sh
import numpy as np
import pandas as pd
import itertools

def scorecommongroups(includemixed = 0):

    scorefn = sh.scorecluster_wrtsystem
    weights = [1,0,-1]

    # all human data
    hdfile = '../output/data/allhuman_bright_ourids.txt'   

    with open(hdfile, 'r') as f:
            hd = json.load(f)

    cnames, csystems = hd[0], hd[1]

    csys_sets = [sh.llist2lset(sys) for sys in csystems]
    cscores = sh.consensusscores(csys_sets, scorefn, weights, np.mean)
    cscoresfull = sh.consensusscores(csys_sets, scorefn, weights, lambda x: x)

    modelfile = '../output/results/model_n_F10.pickle'
    mscores = sh.consensus_modelscores(modelfile, csys_sets, scorefn, weights)

    starfile = '../output/data/stars.txt'
    with open(starfile, "rb") as f:
        stars = json.load(f)

    # convert list to tuple so we can easily remove duplicates
    clist = list(map(lambda x: tuple(sorted(x)), map(set, itertools.chain.from_iterable(csystems))))
    scorelist = list(itertools.chain.from_iterable(cscores))
    mscorelist = list(itertools.chain.from_iterable(mscores))
    fullscorelist = list(map(tuple, itertools.chain.from_iterable(cscoresfull)))

    cscores = list(zip(scorelist, mscorelist, clist, fullscorelist))
    # compute statistics reported in paper
    paperpercent = list(zip(scorelist, mscorelist, clist))
    pc = pd.DataFrame(paperpercent, columns = ["humanscore", "modelscore", "concept"])
    pc = pc.sort_values(by=['humanscore', 'concept'], kind ='stable')
    pc['clen'] = pc['concept'].map(len)
    pc.to_csv('../output/results/paperpercent_df.csv', index_label="index")

    cscores = list(set(cscores))
    sortedscores = sorted(cscores, reverse=True)

    comps = list(zip(*sortedscores))
    clist, cscores, mscores, fullscores = comps[2], comps[0], comps[1], comps[3]
    mmscores = tuple([0 for csc in cscores])

    df = pd.DataFrame(list(map(list, fullscores)), columns=cnames)
    df['cluster'] = list(map(str, clist))
    df.to_csv('../output/results/commongroups_df.csv', index_label="index")

    if includemixed:  # replace scores with mixedmodel scores -- need to have run scripts/mixedmodel.R first
        mmscores = pd.read_csv('../output/results/commongroups_inc.csv', names=["index", "mmscore"])
        mmscores = mmscores.mmscore.values
    else:
        mmscores = [0 for c in cscores]

    cscores_mm = list(zip(cscores, mmscores, mscores, clist))
    sortedscores_mm = sorted(cscores_mm, reverse=True)
    comps = list(zip(*sortedscores_mm))
    clist, cscores, mmscores, mscores = comps[3], comps[0], comps[1], comps[2]
    # sh.printcscores(stars, clist, cscores, mmscores, mscores, outfile='../output/results/TMPTEST')

    # sort by mmscores to check which groups have highest adjusted scores
    # cscores_mm = list(zip(mmscores, cscores, mscores, clist))
    # sortedscores_mm = sorted(cscores_mm, reverse=True)
    # comps = list(zip(*sortedscores_mm))
    # clist, cscores, mmscores, mscores = comps[3], comps[1], comps[0], comps[2]
    # sh.printcscores(stars, clist, cscores, mmscores, mscores, outfile='../output/results/TMPTEST')

    clset = sh.llist2lset(clist)
    finalc, finalscore, finalmmscore, finalmscore = [], [], [], []
    for ci, c in enumerate(clset):
        prevmatch = scorefn(c, clset[0:ci], weights=weights)
        if prevmatch < 0.5 and cscores[ci] > 0.01:
            finalc.append(clist[ci])
            finalscore.append(cscores[ci])
            finalmscore.append(mscores[ci])
            finalmmscore.append(mmscores[ci])

    sh.printcscores(stars, finalc, finalscore, finalmmscore, finalmscore, outfile = '../output/results/tableS2_clusters.tex')
