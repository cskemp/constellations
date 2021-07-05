import numpy as np
import itertools
import scipy.optimize as spo
import pickle
import json
import os
import re
import functions.plotmodelsoln_single as plotmodelsoln_single
import functions.starhelpers as sh
import functions.models as models
import functions.parmaps as parmaps
import functions.readgroups as readgroups
import functions.modelspecs as modelspecs
import functions.modelrunhelpers as modelrunhelpers


# main function for running models

def runmodel(modelnames = ['model_n'], scorefns = ['F10'], optimize=0, weights = [1,0,-1], showgraph = 0, showmodelgraph = 0, listclusters=0, printedges=0, gridinit=0, showgriditems=0, modelhumanprofile=0, fscoretype = 'sys'):

    # ----------------------------------------------------------------------
    # Files and paths
    starfile = '../output/data/stars.txt'
    optfile = '../output/data/allhuman_opt_bright_ourids.txt'  # use when optimizing parameters
    outdir = '../output/results/'
    # ----------------------------------------------------------------------


    # create mapping from hrid to our indices
    invhd = sh.invhdict(starfile)
    with open(starfile, "rb") as f:
        stars = json.load(f)

    fixp, pgrids, modelfs, parnamess, guesss, gridss = modelspecs.modelspecs(starfile, optfile)
    printnames = modelspecs.printnames
    if fscoretype == 'sys':
        fscores = modelrunhelpers.fscores_sys  
    elif fscoretype == 'corr':
        fscores = modelrunhelpers.fscores_corr  
    else: 
        raise Exception('unknown fscoretype')

    pm, invpm = parmaps.parmaps()

    scorefm = {}
    scorefm['adj_rand'] = lambda x, y: sh.scoreworld_arand(x, y, starfile=starfile, bval=pm['bmap'](fixp['b']))
    scorefm['F'] = lambda x, y: sh.scoreworld(x, y, weights=weights, beta=1)
    scorefm['F10'] = lambda x, y: sh.scoreworld(x, y, weights=weights, beta=10)

    for scoref in scorefns:
        sf = scorefm[scoref]
        for modelname in modelnames:
            bestmcss = None
            modelf = modelfs[modelname]
            parnames = parnamess[modelname]
            guess = guesss[modelname]

            grids = []
            for par in parnames:
                grids.append(pgrids[par])

            if optimize == 0:
                bestpars = guess
            elif optimize == 1:
                outfile = outdir + modelname + '_' + scoref + '_opt.pickle'
                if gridinit == 1:
                    gridfile = outdir + modelname + '_' + scoref + '.pickle'
                    with open(gridfile, 'rb') as f_gridfile:
                        gridresults = pickle.load(f_gridfile)
                    bestparind = np.argmax(gridresults['scores'])
                    guess = list(gridresults['pvals'][bestparind])

                if os.path.isfile(outfile):
                    with open(outfile, 'rb') as f_outfile:
                        outd = pickle.load(f_outfile)
                    bestpars = outd['bestpars']
                    bestscore = np.round(outd['fun'], 4)
                    print('best score:' + str(bestscore))
                else:
                    # bnds = [(0.00000001, 100000) for i in guess]  # Powell can't handle bounds
                    od = {'disp': True, 'xtol': 0.001, 'ftol': 0.001}
                    methname = 'Powell'
                    allargs = (modelf, starfile, optfile, weights, invhd, sf, fixp, 0, fscores),
                    result = spo.minimize(fun=modelrunhelpers.fwrap,
                                          args=(modelf, starfile, optfile, weights, invhd, sf, fixp, 0, fscores),
                                          x0=guess,
                                          method=methname,
                                          options=od)
                    bestpars = result['x']
                    outd = {'modelname': modelname, 'guess': guess, 'bestpars': bestpars, 'result': result,
                            'fun': result['fun'], 'args': str(allargs), 'method': methname, 'options': od}
                    modelrunhelpers.savegridresults(outfile, outd)
            elif optimize == 2:
                print(modelname)
                outd = {'modelname': modelname}
                outd['pvals'] = list(itertools.product(*grids))
                npvecs = len(outd['pvals'])
                mcs, stargraph, score, actualpars, scorestr, allscores = fscores(guess, modelf, starfile, optfile, weights,
                                                                                 invhd, sf, fixp, 1)
                origstargraph = stargraph
                outd['allscoress'] = np.zeros((npvecs, len(allscores)))
                outd['scores'] = np.zeros(npvecs)
                outd['mcss'] = [[] for i in outd['scores']]
                outfile = outdir + modelname + '_' + scoref + '.pickle'
                # hcs = readgroups.group(optfile, invhd)
                hss = readgroups.cultures(optfile, invhd)
                # if outfile already exists (possibly incomplete) then load it
                if os.path.isfile(outfile):
                    with open(outfile, 'rb') as f_outfile:
                        outd = pickle.load(f_outfile)
                if not outd['mcss'][-1]:  # only run through outd['pvals'] if mcss incomplete
                    for i in np.arange(npvecs):
                        if not outd['mcss'][i]:
                            ps = list(outd['pvals'][i])
                            if parnames == ['n']:  # don't need to recreate origstargraph -- just threshold differently
                                newg = origstargraph.copy()  # because running the model removes some edges in stargraph
                                topnc = models.topn_clusters(newg, pm['nmap'](ps[0]), alledge=0)
                                mcs = [[]] + topnc[0]
                                stargraph = topnc[1]
                                rawscore, doscore, allscores = sf(mcs, hss)
                                score, allscores, scorestr = modelrunhelpers.score_format(rawscore, allscores)
                            else:
                                mcs, stargraph, score, actualpars, scorestr, allscores = fscores(ps, modelf, starfile,
                                                                                                 optfile, weights, invhd,
                                                                                                 sf, fixp, 0)
                                print(actualpars)
                            print(str(ps) + ' ' + scorestr)
                            outd['scores'][i] = score
                            outd['allscoress'][i] = allscores
                            outd['mcss'][i] = mcs
                            if showgriditems:
                                # need to work on code for plotting edges in right colors
                                me = sh.allmodeledges(mcs, stargraph)
                                nval = pm['nmap'](ps[0])
                                plotfile = '../output/figures/videoframes/frame_' + str(i)
                                plotmodelsoln_single.plotmodelsoln(mcs, me, stargraph, '', title='n = ' + str(nval),
                                                                   samecol=1, edgecol=0, show=0)
                                plotmodelsoln_single.plotfullwithexamples(mcs, me, stargraph, flip=1,
                                                                          prefix=plotfile + '_egs',
                                                                          title='n = ' + str(nval), show=0)
                            if np.mod(i, 1000) == 0:
                                print(i)
                                # modelrunhelpers.savegridresults(outfile, outd)
                    modelrunhelpers.savegridresults(outfile, outd)

                bestscore = np.max(outd['scores'])
                print('best score:' + str(bestscore))
                bestparind = np.argmax(outd['scores'])
                bestpars = list(outd['pvals'][bestparind])
                bestmcss = outd['mcss'][bestparind]

            # Now work with best model solution
            parstr = ''
            for xi, xval in enumerate(list(bestpars)):
                pname = parnames[xi]
                parstr = parstr + pname + ': ' + str(round(pm[pname + 'map'](xval), 2)) + ' '
            print(modelname + ': ' + parstr)

            if showgraph:
                mcs, stargraph, score, actualpars, scorestr, _ = fscores(bestpars, modelf, starfile, optfile, weights,
                                                                         invhd, sf, fixp, 0)
                if bestmcss:
                    mcs = bestmcss  # use best soln from optimization sweep

                newmc = [[]] + mcs
                newme = []
                newme = sh.allmodeledges(newmc, stargraph)
                plottitle = modelname + ': ' + parstr + ': score=' + scorestr
                plottitle = ''
                plotfile = ''
                plotfile = '../output/figures/best_' + modelname + '_' + str(optimize)
                if (re.match("model_(knn|code)", modelname)):
                    sh.addecount(stargraph, const=1)
                else:
                    sh.addecount(stargraph)
                magthresh = pm['bmap']([np.max(fixp['b'])])[0]
                plotmodelsoln_single.plotmodelsoln(newmc, newme, stargraph, plotfile, title=plottitle, magthresh=magthresh,
                                                   edgecol=0, edgecount=1, showlabs=0, show=0)
                plotmodelsoln_single.plotfullwithexamples(newmc, newme, stargraph, flip=1, prefix=plotfile + '_egs',
                                                          edgecount=1, magthresh=magthresh)

                print(plottitle)
                if 0:  # filter out pairs
                    newmc_fp, newme_fp = sh.filterpairs(newmc, newme, stargraph)
                    plotfile = plotfile + '_filterpair'
                    plotmodelsoln_single.plotmodelsoln(newmc_fp, newme_fp, stargraph, plotfile, title=plottitle,
                                                       magthresh=magthresh, edgecol=0, edgecount=1, showlabs=0)

            if 0:  # plot mst
                newmc = [[]] + [set(itertools.chain(*stargraph.edges()))]
                sh.addinvweight(stargraph)
                newme = [[]] + [nx.minimum_spanning_tree(stargraph, weight='invweight')]
                plottitle = modelname + ': ' + parstr + ' mst'
                plotmodelsoln_single.plotmodelsoln(newmc, newme, stargraph, title=plottitle,
                                                   magthresh=pm['bmap'](fixp['b']), edgecol=1, showlabs=0)

            if printedges and optimize != 2:
                bval = pm['bmap'](fixp['b'])
                humanfile = '../output/results/allhumanedges_' + str(bval) + '.txt'
                plotfile = '../output/figures/humanvsmodel_' + modelname
                edgefile = '../output/results/modelhumanedges_' + modelname + str(bval) + '.csv'
                _, stargraphall1, _, _, _, _ = fscores(bestpars, modelf, starfile, optfile, weights, invhd, sf, fixp, 1)
                # need to have computed stargraph with alledges=1
                sh.compare_modelhumanedges(humanfile, stargraphall1, modelname=printnames(modelname), edgefile=edgefile,
                                           plotfile=plotfile, printflag = 0)
                # print model graph
                if showmodelgraph:
                    plotfile = '../output/figures/model_graph'
                    plotmodelsoln_single.plotmodelsoln([], [], stargraphall1, plotfile, title='',
                                                   magthresh=pm['bmap']([np.max(fixp['b'])])[0], samecol=1, edgecol=1,
                                                   showlabs=0, ethick = 1)

            if modelhumanprofile:
                prof = sh.profile_system(newmc, newme, stargraph)
                magthresh = pm['bmap'](fixp['b'])
                fpref = '../output/figures/' + modelname + str(magthresh) + '_profile'
                sh.plot_profile(prof, fpref)

            if listclusters:
                # print clusters with human scores
                conncounts = sh.printclusters(mcs, stargraph, humanfile='../output/data/allhuman_bright_ourids.txt', outfile = '../output/results/tableS3_clusters.tex')
