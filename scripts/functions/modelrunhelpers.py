import pickle
import functions.starhelpers as sh
import functions.readgroups as readgroups
import numpy as np
import functions.parmaps as parmaps

pm, invpm = parmaps.parmaps()


def score_format(rawscore, allscores):
    #score = round(rawscore, 2)
    score = rawscore
    allscores = [round(s, 2) for s in allscores]
    #scorestr = str(score) + ' ' + str(allscores)
    scorestr = str(score)
    return score, allscores, scorestr

# run model for a given set of parameters. Score based on comparing training system with single model system
def fscores_sys(x, *args):
    modelf, starfile, optfile, weights, invhd, scoref, fixp, alledge = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]
    # model
    mcs, stargraph, actualpars = modelf(starfile, np.append(x, alledge), fp=fixp)

    # flat list of human groups
    # hcs = readgroups.group(optfile, invhd)
    # rawscore, _, allscores = sh.scoreclusters(mcs, hcs, weights=weights)

    # list of groups for each culture
    hss = readgroups.cultures(optfile)
    #rawscore, _, allscores = sh.scoreworld(mcs, hss, weights=weights)
    rawscore, _, allscores = scoref(mcs, hss)

    score, allscores, scorestr = score_format(rawscore, allscores)
    return mcs, stargraph, rawscore, actualpars, scorestr, allscores


# Alternative based on correlation between human counts, model weights (not tested recently)
def fscores_corr(x, *args):
    modelf, starfile, optfile, weights, invhd, _, fixp, alledge = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]
    # run model. need alledges = 1 to preserve all edges
    mcs, stargraph, actualpars = modelf(starfile, np.append(x,1), fp=fixp)
    bval = pm['bmap'](fixp['b'])
    humanfile = '../output/results/allhumanedges_' + str(bval) + '.txt'

    #newmc = [[]] + mcs
    #newme = sh.mstedges(newmc, stargraph)
    #edgefile = '../output/results/modeledges_' + str(bval) + '.txt'
    #sh.printedgecounts(stargraph, newme, edgefile, att='weight')

    #rawscore = sh.compare_modelhumanedges_hfocused(humanfile, stargraph, printflag=0)

    rawscore = sh.compare_modelhumanedges(humanfile, stargraph, printflag=0)
    score = round(rawscore, 4)
    allscores = [score]
    return mcs, stargraph, rawscore, actualpars, str(score) + ' ' + str(allscores), allscores

# Alternative based on ml fit to human counts
def fscores_ml(x, *args):
    modelf, starfile, optfile, weights, invhd, _, fixp, alledge = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]
    # run model. need alledges = 1 to preserve all edges
    mcs, stargraph, actualpars = modelf(starfile, np.append(x,1), fp=fixp)
    bval = pm['bmap'](fixp['b'])
    humanfile = '../output/results/allhumanedges_' + str(bval) + '.txt'
    rawscore = sh.compare_modelhumanedges_ml(humanfile, stargraph, printflag=0)
    score = round(rawscore, 4)
    allscores = [score]
    return mcs, stargraph, rawscore, actualpars, str(score) + ' ' + str(allscores), allscores


# wrapper that discards everything except model score
def fwrap(x, *args):
    if np.min(x) < 0:  # keep parameters positive
        return 1000000000
    fscores = args[-1]
    restargs = args[0:-1]
    _, _, score, actualpars, scorestr, _, = fscores(x, *restargs)
    print((round(score, 4), actualpars, scorestr))
    return -score

# save grid results
def savegridresults(outfile, outd):  #modelname, grids, npvecs, mcss, scores, allscoress):
    with open(outfile, "wb") as f:
        pickle.dump(outd, f)
