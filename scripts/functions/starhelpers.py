import numpy as np
import pandas as pd
import json
import functions.neighborhood as neighborhood
import networkx as nx
import itertools
from scipy.stats import rankdata
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import statsmodels.api as sm
from skimage.future import graph
import sklearn.neighbors
import pickle
import re
import sys 
import collections
import pandas as pd
from scipy.stats import halfgennorm, lognorm, expon

# sum two star magnitudes
def magsum(m1, m2):
    s = -2.5 * np.log10(10 ** (-m1 * 0.4) + 10 ** (-m2 * 0.4))
    return s


# thresh should be list
def makestargraph(starfile, thresh, outprefix=''):
    with open(starfile, "rb") as f:
        stars = json.load(f)
    suff = "_" + str(thresh)
    #nb = neighborhood.makeneighborhood(stars, [thresh])
    nb = neighborhood.makeneighborhood(stars, thresh)
    G, _ = makedistgraph(nb, stars)
    if outprefix != '':
        outfile = outprefix + suff + '.graphml'
        nx.write_graphml(G, outfile)
    return G

def makestargraph_closenb(starfile, thresh, maxd, outprefix=''):
    with open(starfile, "rb") as f:
        stars = json.load(f)
    suff = "_closenb_" + str(maxd) + '_' + str(thresh)
    nb = neighborhood.makecloseedges_nb(stars, thresh, maxd)
    G, _ = makedistgraph(nb, stars)
    if outprefix != '':
        outfile = outprefix + suff + '.graphml'
        nx.write_graphml(G, outfile)
    return G

def namedict(starfile):
    with open(starfile, "rb") as f:
        stars = json.load(f)
    named= {}
    for (i, (h,r,d,m,n)) in enumerate(stars):
        named[i] = n
    return named 


def invhdict(starfile):
    with open(starfile, "rb") as f:
        stars = json.load(f)
    invhd = {}
    for (i, (h,r,d,m,n)) in enumerate(stars):
        invhd[h] = i
    return invhd

def invbayerd(starfile):
    with open(starfile, "rb") as f:
        stars = json.load(f)
    invbayerd = {}
    for (i, (h,r,d,m,n)) in enumerate(stars):
        invbayerd[n] = i
    return invbayerd



def starinddict(stars):
    sd = {}
    for si, s in enumerate(stars):
        sd[s[0]] = si
    return(sd)


def name2starind(sg, name):
    nodeid = [e['i'] for n, e in sg.nodes(data=True) if e['name'] == name]
    return nodeid[0]


def mstedges(clusters, stargraph, ename = 'd'):
    msts = [[]]
    for c in clusters[1:]:
        msts.append( nx.minimum_spanning_tree(stargraph.subgraph(c), weight=ename ) )
    return msts

def allmodeledges(clusters, stargraph):
    msts = [[]]
    for c in clusters[1:]:
        msts.append( stargraph.subgraph(c) )
    return msts

# filter out all stars that don't make brightness cutoff
def brightclusters(cc, stargraph, brightthresh):
    bcs = []
    for c in cc:
        bc = [ci for ci in list(c) if stargraph.nodes[ci]['m']<= brightthresh]
        #if bc:
        if len(bc) > 1: # drop asterism unless 2 or more stars remain
            bcs.append(set(bc))
    return bcs


# print names for stars in clusters
def printclusters(cc, stargraph, humanfile = '', outfile = ''):
    cc = np.array(cc)
    scores = np.array([0. for c in cc])
    matchnames = [ '' for c in cc]

    handle = open(outfile, "w") if outfile else sys.stdout

    if humanfile:
        scorefn = scorecluster_wrtsystem
        weights = [1, 0, -1]
        
        with open(humanfile, "rb") as hf:
            hall = json.load(hf)
        hnames = np.array(hall[0])
        hsys = hall[1]

        for ci, c in enumerate(cc):
            hsysscores = np.array([])
            for hs in hsys:
                csc = scorefn(c, hs, weights=weights)
                hsysscores= np.append(hsysscores, csc)

            scores[ci] = np.max(hsysscores)
            mns = hnames[np.argwhere(hsysscores== np.amax(hsysscores)).flatten()]
            matchnames[ci] = ','.join(mns)

        sortind = scores.argsort()
        sortind = sortind[::-1]
        scores = scores[sortind]
        matchnames = np.array(matchnames)
        matchnames = matchnames[sortind]
        cc = cc[sortind]       

    conn_count, not_conn_count = 0, 0
    for ci, c in enumerate(cc):
        namelist = sortnamelist([stargraph.node[ci]['name'] for ci in c])
        subg = stargraph.subgraph(c)
        # move next line inside following if statement to get disconnected groups only
        handle.write( str(ci+1) + ' & ' + str(round(scores[ci],2)) + ' & ' + ', '.join(namelist) + ' & ' + matchnames[ci] + ' \\\\\n')
        if c and not nx.is_connected(subg):
            handle.write('    ** ^ is not connected\n')
            not_conn_count = not_conn_count + 1
        elif c:
            conn_count = conn_count + 1

    if handle is not sys.stdout:
        handle.close()

    return( conn_count, not_conn_count)

# print consensusscores
def printcscores(stars, cc, cscores, mmscores, mscores, outfile = '' ):

    handle = open(outfile, "w") if outfile else sys.stdout

    sd = starinddict(stars)
    for ci, c in enumerate(cc):
        #namelist = [stars[sd[ci]][4] for ci in c if ci in sd]
        namelist = sortnamelist([stars[ci][4] for ci in c])
        cscore = round(cscores[ci], 2)
        mmscore = round(mmscores[ci], 2)
        mscore = round(mscores[ci], 2)
        handle.write(str(ci+1) + ' & ' + str(cscore) + ' & '  + str(mmscore) + ' & '+ str(mscore) +' & ' +  ', '.join(namelist) + ' &  \\\\\n')

    if handle is not sys.stdout:
        handle.close()


# make networkx graph from sparse matrix, list of stars
def makedistgraph(nb, stars):
    G = nx.from_scipy_sparse_matrix(nb, edge_attribute='d')

    for i, (h, r, d, m, n) in enumerate(stars):
        G.add_node(i, i=i, m=m, ra=r, dec=d, h=h, name=n)

    maxd = 0;
    for v1, v2 in G.edges():
        dist = neighborhood.angsep((G.nodes[v1]['ra'], G.nodes[v1]['dec']), (G.nodes[v2]['ra'], G.nodes[v2]['dec']))
        if dist > maxd:
            maxd = dist
        G.add_edge(v1, v2, d=dist)

    return G, maxd

def dropsmalledges(sg, att, smallval, alledge=0):
    smalledges = [(u, v) for u, v, e in sg.edges(data=True) if e[att] < smallval]
    if alledge:
        newsg = sg.copy()
    else:
        newsg = sg
    newsg.remove_edges_from(smalledges)
    return newsg

# set brightest star to magnitude 0
def standardizemag(sg):
    mdf = pd.DataFrame.from_dict(dict(sg.nodes.data('m')), orient='index', columns=['m'])
    minmag = np.min(mdf.values)
    for node in sg.nodes():
        currmag = sg.nodes[node]['m']
        sg.add_node(node, mstd=currmag - minmag)
    return

# add invweight
def addinvweight(sg):
    for u, v, e in sg.edges(data=True):
        sg.add_edge(u, v, invweight = -e['weight'])
    return

# add edges with avmag
def addavmagedges(sg):
    for u, v in sg.edges():
        mstdu = sg.nodes[u]['m']
        mstdv = sg.nodes[v]['m']
        sg.add_edge(u, v, avmag=0.5 * (mstdu + mstdv))
    return

# add edges with maxmag
def addmaxmagedges(sg):
    for u, v in sg.edges():
        mstdu = sg.nodes[u]['m']
        mstdv = sg.nodes[v]['m']
        sg.add_edge(u, v, maxmag=np.max([mstdu,mstdv]))
    return

# add m_lr (likelihood ratio for magnitudes)
def add_m_lr(mf_lr, sg):
    for n in list(sg.nodes):
        sg.add_node(n, m_lr= mf_lr(sg.nodes[n]['m']))
    return

# add edges with min_m_lr
def add_min_m_lr_edges(sg):
    for u, v in sg.edges():
        mstdu = sg.nodes[u]['m_lr']
        mstdv = sg.nodes[v]['m_lr']
        sg.add_edge(u, v, min_m_lr=np.min([mstdu,mstdv]))
    return

def add_d_lr(df_lr, sg):
    for u, v in sg.edges():
        d = sg.edges[u, v]['d']
        sg.add_edge(u, v, d_lr=df_lr(d))
    return

# convert distance to similarity measure. l is weight for exponential transformation
def dist2sim(sg, att, l):
    for u, v, e in sg.edges(data=True):
        dist = e[att]
        sim = np.exp(-dist * l)
        newdic = {att + 'sim': sim}
        sg.add_edge(u, v, **newdic)
    return

def weightatt(sg, att, l):
    for u, v, e in sg.edges(data=True):
        worig = e[att] ** l
        newdic = {att + 'weight': worig}
        sg.add_edge(u, v, **newdic)
    return

# think that this scaling makes no difference
def scaleedgeweights(sg, att1, att2):
    return

# combine edge atts
def combineedgeatts(sg, atts, cname, ctype='sum'):
    if ctype == 'sum':
        cinit = 0
        cfunc = lambda x, y: x + y
    elif ctype == 'prod':
        cinit = 1
        cfunc = lambda x, y: x * y
    else:
        raise Exception('unexpected combination type: ' + ctype)

    for u, v, e in sg.edges(data=True):
        cval = cinit
        for att in atts:
            cval = cfunc(cval, e[att])
        if cval < 1e-20:
            cval = 1e-20
        newdic = {cname: cval}
        sg.add_edge(u, v, **newdic)

    return

# make weight attribute
def setweight(sg, att):
    for u, v, e in sg.edges(data=True):
        sg.add_edge(u, v, weight=e[att])
    return

# dilate node atts based on local neighborhood
def dilatenodes(sg, dnbgraph, att, dtype='orig'):
    # 'orig':  original value
    # 'scale_typ':  scale original value wrt to typical value in local neighborhood
    # 'lval' : set to value in local neighborhood
    # 'gval':  set to global value
    # 'scale_tot':  scale original value wrt to total value in local neighborhood

    alpha = 1.0
    newatt = att + '_scale'
    newattlocw = att + '_locw'
    # only include nodes above brightthresh -- these are the only ones with e['nn']
    globvals = [e[att] for _, e in sg.nodes(data=True) if e['nn']]
    if dtype == 'scale_tot':
        globval = np.sum(globvals)
    else:
        globval = np.median(globvals)

    for u, e in sg.nodes(data=True):
        currattval = e[att]
        if e['nn']:
            # make list of nodes that are within neighborhood of u
            uns = set([u] + [un for _, un in dnbgraph.edges(u)])
            nbvals = [sg.node[un][att] for un in uns]
            if dtype == 'scale_tot':
                nbval = np.sum(nbvals)
            else:
                nbval = np.median(nbvals)

            if dtype=='orig':
                newattval = currattval
            elif dtype =='scale_typ' or dtype == 'scale_tot':
                newattval = currattval * (globval/nbval)
            elif dtype == 'lval':
                newattval = nbval
            elif dtype =='gval':
                newattval = globval
            else:
                raise Exception('unknown dilation type')

            newdic = {newatt: newattval, newattlocw:(nbval/globval) ** alpha}
            sg.add_node(u, **newdic)
    return

# dilate edges based on local neighborhood
def dilateedges(sg, dnbgraph, att, newatt):
    globmedian = np.median([e[att] for u, v, e in sg.edges(data=True)])
    for u, v, e in sg.edges(data=True):
        currattval = e[att]
        # make list of nodes that are within neighborhood of u and v
        unb = set([u, v] + [u1 for _, u1 in dnbgraph.edges(u) if dnbgraph.has_edge(v, u1)])
        nbattvals = []
        for unbnode in unb:
            nbattvals  = nbattvals + [sg.get_edge_data(unbnode, nbn2)[att] for _, nbn2 in sg.edges(unbnode) if nbn2 in unb]
        nbmedian = np.median(nbattvals)
        newattval = currattval * (globmedian/nbmedian) ** 1
        newdic = {newatt: newattval}
        sg.add_edge(u, v, **newdic)
    return

# dilate edges based on extremely local neighborhood
def dilateedges_pp(sg, att, newatt):
    nodeatt = att + 'node'
    for n in sg.nodes():
        nedges = sg.edges(n)
        nvals = [sg[u][v][att] for u, v in nedges]
        if nvals:
            maxval = np.max(nvals)
            newdic = {nodeatt : maxval}
            sg.add_node(n, **newdic)

    for u, v, e in sg.edges(data=True):
        orig = e[att]
        new = orig / np.sqrt(sg.node[u][nodeatt] * sg.node[v][nodeatt])
        newdic = {newatt: new}
        sg.add_edge(u, v, **newdic)

    return

# alpha:  probability at each step that random walk continues
# thresh: cluster includes nodes with pagerank > thresh
def personalpr(sg, ref_node, alpha=0.9, thresh=0.05):
    nnode = sg.number_of_nodes()
    refnoded = dict(zip(range(nnode), itertools.repeat(0, nnode)))
    refnoded[ref_node] = 1

    pr = nx.pagerank(sg, alpha=alpha, personalization=refnoded, weight='weight')
    prframe = pd.DataFrame.from_dict(pr, orient='index')
    output_pr = prframe[0].values
    output_weighted = output_pr > thresh
    lc = np.where(output_weighted)
    return set(lc[0].flatten()), output_pr


# https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
# convert fragments (e.g. from personalized pagerank) to single graph with clusters

def frag_to_graph(l):
    G = nx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also implies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """ 
        treat `l` as a Graph and returns its edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

# score high if mcs contains something like each cluster in hcs
def scoreclusters(mcs, hcs, weights=[+1, 0, -1]):
    allscores = [0 for hc in hcs]
    for hi, hc in enumerate(hcs):
        allscores[hi] = scoresystem_wrtcluster(hc, mcs, weights)
                                        # bc first human cluster is empty
    #score = np.sum(allscores)
    score = np.mean(allscores)
    droponescore = np.sum(allscores) - np.min(allscores[1:])
    #allscores = [round(a,1) for a in allscores]

    return score, droponescore, allscores

# score high if hcs contains something like each cluster in mcs
def scoreclusters_mod(mcs, hcs, weights=[+1, 0, -1]):
    allscores = [0 for mc in mcs]
    mcs = llist2lset(mcs)
    for mi, mc in enumerate(mcs):
        allscores[mi] = scorecluster_wrtsystem(mc, hcs, weights)
                                        # bc first human cluster is empty
    #score = np.sum(allscores)
    score = np.mean(allscores)
    droponescore = 0
    #allscores = [round(a,1) for a in allscores]

    return score, droponescore, allscores

# return vector of scores for each culture
def scoreworld(mcs, hcslist,  weights=[+1, 0, -1], beta=10):
    allscores = [0 for hc in hcslist]
    for hi, hcs in enumerate(hcslist):
        # like recall
        a1, _, _ = scoreclusters_mod(mcs, hcs, weights)
        # like precision
        a2, _, _ = scoreclusters(mcs, hcs, weights)
        if (a1 + a2) == 0:
            allscores[hi] = 0
        # analogous to F-score
        else:
            allscores[hi] = (1 + beta**2) * a1 * a2 / (beta**2*a1 + a2)

    score = np.mean(allscores)
    droponescore = np.sum(allscores) - np.min(allscores[1:])
    #allscores = [round(a,1) for a in allscores]

    return score, droponescore, allscores

# return vector of scores for each culture
def scoreworld_arand(mcs, hcslist, starfile='', bval=[]):
    stargraph = makestargraph(starfile, bval)
    allcs = list(nx.connected_components(stargraph))
    cs = [c for c in allcs if len(c) > 1]
    bigc = cs[0]
    N = stargraph.number_of_nodes()
    mcs_vec = csys2vec(mcs, bigc, N)
    allscores = [0 for hc in hcslist]
    for hi, hcs in enumerate(hcslist):
        hcs_vec = csys2vec(hcs, bigc, N)
        allscores[hi] = sklearn.metrics.adjusted_rand_score(hcs_vec, mcs_vec)

    score = np.mean(allscores)
    droponescore = np.sum(allscores) - np.min(allscores[1:])
    return score, droponescore, allscores

# with default weights it's better for mc to miss stars in hc than to have extra stars
def scorecluster(mc, hc, weights):
    m = set(mc)
    h = set(hc)
    mnoth = m - h
    hnotm = h - m
    mandh = m & h
    score = weights[0] * len(mandh) + weights[1] * len(hnotm) + weights[2] * len(mnoth)
    maxscore = weights[0]*len(hc)
    # finalscore between -1 and 1.
    #finalscore = np.max([score/maxscore, -1])
    finalscore = np.max([score/maxscore, 0])
    return finalscore

def llist2lset(csystem):
    lset = [ set(c) for c in csystem]
    return lset

def csys2vec(csystem, bigc, N):
    all = np.zeros(N)
    for ci, c in enumerate(csystem):
        all[list(c)] = ci+1
    return(all[list(bigc)])

# score a system wrt a cluster
# c is gold standard -- want to give csystem a high score if it contains a cluster that captures c well

def scoresystem_wrtcluster(c, csystem, weights=[1,0,-1]):
    scores = []
    for sysc in csystem:
        if c.intersection(sysc):
            sysc_score = scorecluster(sysc, c, weights)
            scores.append(sysc_score)
    if scores:
        score = np.max(scores)
    else:
        score = 0
    return(score)


# score a cluster wrt a system
# csystem is gold standard -- want to give c a high score if it is like a cluster in csystem

def scorecluster_wrtsystem(c, csystem, weights=[1,0,-1]):
    scores = []
    for sysc in csystem:
        if c.intersection(sysc):
            sysc_score = scorecluster(c, sysc, weights)
            scores.append(sysc_score)
    if scores:
        score = np.max(scores)
    else:
        score = 0
    return(score)

def cluster_insystem(c, csystem, **kwargs):
    score = 0
    for sysc in csystem:
        if c == sysc:
            score = 1
            break
    return score

def cluster_subsetsystem(c, csystem, **kwargs):
    score = 0
    for sysc in csystem:
        if c.issubset(sysc):
            score = 1
            break
    return score

# scorefn could be 
#   scorecluster_wrtsystem
#   cluster_insystem
#   cluster_partofsystem

def consensus_modelscores(modelfile, csystems, scorefn, weights = [1,0,-1]):
    with open(modelfile, 'rb') as f_outfile:
        modelsys = pickle.load(f_outfile)
        modelsys = modelsys['mcss']

    scores = compare_humantomodelsys(modelsys, csystems, scorefn, weights)
    return(scores)


def compare_humantomodelsys(modelsys, csystems, scorefn, weights=[1,0,-1]):
    scores = [ [] for cs in csystems ]
    for sysi, csystemi in enumerate(csystems):
        cscores = []
        for c in csystemi:
            scores_wrtsys = []
            for sysj, csystemj in enumerate(modelsys):
                    scores_wrtsys.append(scorefn(c, csystemj, weights=weights))
            # find best match among model systems
            cscores.append(np.max(scores_wrtsys))
        scores[sysi] = cscores
    return(scores)

def consensusscores(csystems, scorefn, weights = [1,0,-1], summfun = np.mean):
    scores = [ [] for cs in csystems ]
    for sysi, csystemi in enumerate(csystems):
        cscores = [] 
        for c in csystemi:
            scores_wrtsys = []
            for sysj, csystemj in enumerate(csystems):
                # previously  didn't consider self-matches
                #if sysi != sysj:
                scores_wrtsys.append(scorefn(c, csystemj, weights=weights))
            cscores.append(summfun(scores_wrtsys))
        scores[sysi] = cscores
    return(scores)

def comparemodelhuman(mcs, humanfile, starfile, weights = [1,0,-1]):
    named = namedict(starfile)
    bw = 0.1
    histbins = np.arange(-1, 1 + bw, bw)
    allscores = []
    with open(humanfile, 'r') as f:
        hd = json.load(f)
        for i, c in enumerate(hd[0]):
            hsys = llist2lset(hd[1][i])
            sysscores = [scoresystem_wrtcluster(hc, mcs, weights=weights) for hc in hsys]
            print('*************')
            print(c)
            for j, sc in enumerate(sysscores):
              clstring = str(round(sc, 2)) + ' ' + '. '.join([named[s] for s in hd[1][i][j]] )
              print(clstring)

            plt.subplot(4, 6, i + 1)
            plt.hist(sysscores, bins=histbins)
            plt.title(c + ' ' + str(round(np.mean(sysscores), 2)))
            allscores.append(np.mean(sysscores))
        print("system mean = " + str(round(np.mean(allscores), 3)))
    plt.show()

# distance, brightness profiles based on edges in constellation minimum spanning trees
def profile_system(clusters, edges, sg):
    ds = []
    bs = []
    maxbs = []
    csizes = []
    allds = []
    allbs = []
    allmaxbs = []

    for c in clusters:
        if c:
            csizes.append(len(c))
            for s in c:
                if len(sg[s]): # XXX: temporarily consider nodes in nb graph only
                    bs.append(sg.nodes[s]['m'])

    for e in edges:
        if e:
            for u, v in e.edges:
                ds.append(sg.edges[u,v]['d'])
                maxbs.append(np.max([sg.nodes[u]['m'], sg.nodes[v]['m']]))

    for s in sg.nodes:
        if len(sg[s]):
            allbs.append(sg.nodes[s]['m'])

    for u, v in sg.edges:
        allds.append(sg.edges[u,v]['d'])
        allmaxbs.append(np.max([sg.nodes[u]['m'], sg.nodes[v]['m']]))

    return (csizes, bs, ds, maxbs, allbs, allds, allmaxbs)

def plot_profile(pf, fpref=''):
    cs = pf[0]
    bs = pf[1]
    ds = pf[2]
    maxbs = pf[3]
    allbs = pf[4]
    allds = pf[5]
    allmaxbs = pf[6]

    plt.subplot(241)
    plt.hist(cs, bins = np.arange(0.5, 25, 1))
    plt.title('sizes')

    plt.subplot(242)
    plt.hist(bs, bins = np.arange(0, 6.0, 0.25))
    #plt.hist(bs, bins = np.arange(0, 1.46, 0.05))
    plt.title('magnitudes')

    plt.subplot(243)
    plt.hist(ds, bins = np.arange(0, 24, 1.0))
    #plt.hist(ds, bins = np.arange(0, 1, 0.05))
    plt.title('edge distances')

    plt.subplot(244)
    plt.hist(maxbs, bins = np.arange(0, 6.0, 0.25))
    #plt.hist(maxbs, bins = np.arange(0, 1.46, 0.05))
    plt.title('edge magnitudes')

    plt.subplot(246)
    plt.hist(allbs, bins = np.arange(0, 6.0, 0.25))
    #plt.hist(allbs, bins = np.arange(0, 1.46, 0.05))
    plt.title('all magnitudes')

    plt.subplot(247)
    plt.hist(allds, bins = np.arange(0, 24, 1.0))
    #plt.hist(allds, bins = np.arange(0, 1, 0.05))
    plt.title('all edge distances')

    plt.subplot(248)
    plt.hist(allmaxbs, bins = np.arange(0, 6.0, 0.25))
    #plt.hist(allmaxbs, bins = np.arange(0, 1.46, 0.05))
    plt.title('all edge magnitudes')
    plt.tight_layout()

    if fpref:
       plt.savefig(fpref + '.pdf', bbox_inches="tight")

    plt.show()
    return

# empirical distributions based on attested clusters, edges and entire neighbourhood graph
def empirical_ds(clusters, edges, sg, maxmag = 4.5):

    pf = profile_system(clusters, edges, sg)

    bs_att =  [-b for b in pf[1]]
    bs_all = [-b for b in pf[4]]

    beta_att, loc_att, scale_att = halfgennorm.fit(bs_att, floc=-maxmag)
    beta_all, loc_all, scale_all = halfgennorm.fit(bs_all, floc=-maxmag)

    mfit_att = halfgennorm( beta=beta_att, loc=loc_att, scale=scale_att)
    mfit_all = halfgennorm( beta=beta_all, loc=loc_all, scale=scale_all)

    loc_att, scale_att = expon.fit(bs_att, floc=-maxmag)
    loc_all, scale_all = expon.fit(bs_all, floc=-maxmag)

    mfit_att = expon( loc=loc_att, scale=scale_att)
    mfit_all = expon( loc=loc_all, scale=scale_all)

    ds_att = pf[2]
    ds_all = pf[5]

    s_att, loc_att, scale_att = lognorm.fit(ds_att)
    s_all, loc_all, scale_all = lognorm.fit(ds_all)

    dfit_att = lognorm( s=s_att, loc=loc_att, scale=scale_att )
    dfit_all = lognorm( s=s_all, loc=loc_all, scale=scale_all)

    mf_att = lambda x : mfit_att.pdf( -x )
    mf_all = lambda x : mfit_all.pdf( -x )

    df_att = dfit_att.pdf
    df_all = dfit_all.pdf

    mf_lr = lambda x: mf_att(x) / mf_all(x)
    df_lr = lambda x: df_att(x) / df_all(x)

    return (mf_lr, df_lr, mf_att, mf_all, df_att, df_all)

def addtripquad(sg, wtripquad, tripflag):
    if tripflag:  # triples
        triples = all_triples(sg, sg)
        #tripscores = [ (np.abs(np.rad2deg(neighborhood.sphericalangle(t, sg)) - 180), t)  for t in triples]
        # angles in radians
        tripscores = [ (np.abs(np.pi - neighborhood.sphericalangle(t, sg)), t)  for t in triples]
        tripscores = [ (np.exp(-a * wtripquad), t)  for (a,t) in tripscores]
        tripscores = [ (np.min([sg.edges[i,j]['dmagprodweight'], sg.edges[j,k]['dmagprodweight']])*wa, (i,j,k)) for (wa, (i,j,k)) in tripscores ]
        for u, v in sg.edges:
            sg.add_edge(u, v, dmagprodtripquad=-np.inf)

        for (t, (i,j,k)) in tripscores:
            if sg.edges[i,j]['dmagprodtripquad'] < t:
                sg.add_edge(i, j, dmagprodtripquad=t)
            if sg.edges[j,k]['dmagprodtripquad'] < t:
                sg.add_edge(j, k, dmagprodtripquad=t)
    else: # quads: ie goodc model
        quads = all_quads(sg,sg)
        goodcscores = [ (neighborhood.goodc(q, sg), q) for q in quads]
        goodcscores= [ (np.exp(-a * wtripquad), q)  for (a,q) in goodcscores]
        goodcscores = [ (np.min([sg.edges[i,j]['dmagprodweight'], sg.edges[j,k]['dmagprodweight'], sg.edges[k,l]['dmagprodweight']])*wa, (i,j,k,l)) for (wa, (i,j,k,l)) in goodcscores]
        for u, v in sg.edges:
            sg.add_edge(u, v, dmagprodtripquad=-np.inf)
        for (t, (i,j,k,l)) in goodcscores:
            if sg.edges[i,j]['dmagprodtripquad'] < t:
                sg.add_edge(i, j, dmagprodtripquad=t)
            if sg.edges[j,k]['dmagprodtripquad'] < t:
                sg.add_edge(j, k, dmagprodtripquad=t)
            if sg.edges[j,k]['dmagprodtripquad'] < t:
                sg.add_edge(j, k, dmagprodtripquad=t)

    return

def triple_angles(g, sg):
    triples = all_triples(g,sg)
    angles = [ neighborhood.sphericalangle(t, sg)  for t in triples]
    return angles

def triple_distances(g, sg):
    triples = all_triples(g,sg)
    distances = [ neighborhood.tripledistance(t, sg)  for t in triples]
    return distances


def quad_goodc(g, sg):
    quads = all_quads(g,sg)
    goodcs = [ neighborhood.goodc(q, sg) for q in quads]
    return goodcs 

def all_triples(g, sg):
    triples = []
    for u, v in g.edges:
        for w in g[v]:
            if w != u:
                if u < w:  # so we don't end up with the same edge forwards and backwards
                    triples.append( (u,v,w) )
                else:
                    triples.append( (w,v,u) )
        for t in g[u]:
            if t != v:
                if t < v:
                    triples.append( (t,u,v) )
                else:
                    triples.append( (v,u,t) )

    # remove duplicates
    triples = np.unique(triples, axis=0)

    return(triples)

def all_quads(g, sg):
    triples = all_triples(g, sg)
    quads = []
    for u, v, w in triples:
        for x in g[w]:
            if x != u and x != v:
                if u < x:  
                    quads.append( (u,v,w,x) )
                else:
                    quads.append( (x,w,v,u) )

        for t in g[u]:
            if t != v and t != w:
                if t < w:
                    quads.append( (t,u,v, w) )
                else:
                    quads.append( (w,v,u, t) )

    # remove duplicates
    quads = np.unique(quads, axis=0)

    return(quads)


def edge_increment(sg, edges):
    nx.set_edge_attributes(sg, 0, 'ecount')
    for egraph in edges:
        if egraph:
            for u, v in egraph.edges:
                curr = sg.edges[u,v]['ecount']
                sg.add_edge(u, v, ecount=curr+1)
    return

def printedgecounts(sg, es, edgefile, att='ecount'):
    with open(edgefile, "w") as f:
        if not len(es):
            for u, v in sg.edges:
                count = sg.edges[u,v][att]
                if count > 0:
                    f.write(','.join([str(u), str(v), str(count) + '\n']))
        else:
            for eg in es:
                if eg:
                    for u, v in eg.edges:
                        count = sg.edges[u,v][att]
                        f.write(','.join([str(u), str(v), str(count) + '\n']))

def graph_from_tuplefile(file, attname='weight'):
    g = nx.Graph()
    data = pd.read_csv(file, header=None)
    tuples = list(data.itertuples(index=False, name=None))
    g.add_weighted_edges_from(tuples, weight=attname)

    return g

def compare_modelhumanedges(humanfile, stargraph, edgefile='', plotfile = '', modelname = 'model', printflag = 1):
    hg = graph_from_tuplefile(humanfile, attname = 'humanw')
    mg = nx.minimum_spanning_tree(stargraph)

    for (u, v, e) in hg.edges(data=True):
        mg.add_edge(u, v, humanw = e['humanw'])

    m = []
    h = []
    pairs = []
    # include all edges in human file, all edges in model file
    for (u, v, e) in mg.edges(data=True):
        thish = e['humanw'] if 'humanw' in e else 0
        thism = stargraph[u][v]['weight']

        m.append(thism)
        h.append(thish)
        pairs.append( (u,v))
    r = np.corrcoef(m, h)
    if printflag:
        print_corr(h, m, pairs, stargraph, modelname, plotfile)
    if edgefile  != '':
        with open(edgefile, "w") as f:
            for i, pair in enumerate(pairs):
                u, v = pair[0], pair[1]
                nu, nv = stargraph.nodes[u], stargraph.nodes[v]
                f.write(','.join([str(nu['h']), str(nv['h']), nu['name'], nv['name'], str(h[i]), str(m[i])])+'\n')

    return(r[0,1])


def compare_modelhumanedges_hfocused(humanfile, stargraph, printflag = 1):
    hg = graph_from_tuplefile(humanfile, attname = 'humanw')

    m = []
    h = []
    pairs = []
    for (u, v, e) in hg.edges(data=True):
        thish = e['humanw']
        thism = stargraph[u][v]['weight']
        m.append(thism)
        h.append(thish)
        pairs.append( (u,v))
    r = np.corrcoef(m, h)

    if printflag:
        print_corr(h, m, pairs, stargraph)

    return(r[0,1])


def compare_modelhumanedges_ml(humanfile, stargraph, printflag = 1):
    hg = graph_from_tuplefile(humanfile, attname = 'humanw')

    m = []
    h = []
    pairs = []
    ll = 0
    for (u, v, e) in hg.edges(data=True):
        thish = e['humanw']
        thism = stargraph[u][v]['weight']
        ll = ll + thish * np.log(thism)
        m.append(thism)
        h.append(thish)
        pairs.append( (u,v))

    return(ll)

def print_corr(h, m, pairs, stargraph, modelname, plotfile):
    rankh = rankdata(h)
    rankm = rankdata(m)
    r = np.corrcoef(m, h)
    rankr = np.corrcoef(rankm, rankh)

    est = sm.OLS(rankh, rankm).fit()
    residpairs = sorted(zip(est.resid, pairs))
    print('----- model likes but humans dont')
    for rp in residpairs[:10]:
        pair = rp[1]
        print(str(round(rp[0],2)) + ':' + stargraph.node[pair[0]]['name'] + ',' + stargraph.node[pair[1]]['name'])
    print('----- humans likes but model doesnt')
    for rp in residpairs[-20:]:
        pair = rp[1]
        print(str(round(rp[0],2)) + ':' + stargraph.node[pair[0]]['name'] + ',' + stargraph.node[pair[1]]['name'])

    ax = plt.subplot(231)
    sns.set_style("white")
    ax.scatter(m, h, marker='.', color='black' )
    ratio = 1.0
    ax.set_aspect(1.0 / ax.get_data_ratio() * ratio)
    ax.set_xlabel(modelname)
    ax.set_ylabel('human')
    ax.set_title('r = ' + str(round(r[0,1], 2)))

    if plotfile:
        plt.savefig(plotfile + '.eps', format="eps", bbox_inches="tight")

    if 0:
        ax = plt.subplot(332)
        ax.scatter(rankm, rankh)
        ratio = 1.0
        ax.set_aspect(1.0 / ax.get_data_ratio() * ratio)
        ax.set_xlabel('model')
        ax.set_ylabel('human')
        ax.set_title('r = ' + str(round(rankr[0,1], 2)))

    plt.show(block=False)

    return

# decrement with wrap around
def mydec(i, n):
    if i == 0:
        return n-1
    else:
        return i-1

# figure out which grid nodes each star should attach to
def startogrid(sps, decs, ras):
    stog = []
    nd = len(decs)
    nr = len(ras)
    dimatrix, rimatrix = np.meshgrid(np.arange(nd), np.arange(nr))
    # flipping the order is correct because of way nodes are ordered in grid graph
    dorder = np.ravel(dimatrix.T)
    rorder = np.ravel(rimatrix.T)

    for sp in sps:
        di = np.searchsorted(decs, sp[0])
        ri = np.searchsorted(ras, sp[1])
        if di == nd:
            di = 0
        if ri == nr:
            ri = 0
        dis = (mydec(di, nd), di)
        ris = (mydec(ri, nr), ri)
        fourpairs = np.array(np.meshgrid(dis, ris)).T.reshape(-1, 2)
        row = []
        for p in fourpairs:
            row.append( np.where( (dorder==p[0]) & (rorder==p[1]) )[0][0] )
        stog.append( row )

    return stog

def make_grid(sg, brightthresh, nd, nr):
    labels = np.reshape(np.arange(nd * nr), (nd, nr))
    ng = nx.grid_2d_graph(nd, nr, periodic=True)

    dstep = np.pi / nd
    rstep = 2 * np.pi / nr

    ras = np.arange(0, 2 * np.pi, rstep)
    decs = np.arange(-np.pi / 2, np.pi / 2, dstep)

    dimatrix, rimatrix = np.meshgrid(np.arange(nd), np.arange(nr))

    for (di, ri) in zip(np.ravel(dimatrix), np.ravel(rimatrix)):
        ng.add_node((di,ri), dec=decs[di], ra=ras[ri], m=6.1)
    # order inherited from sorted(ng.nodes())
    ng = nx.convert_node_labels_to_integers(ng, first_label=0, ordering='default', label_attribute=None)
    brightsg = sg.subgraph( [n for n,attrdict in sg.node.items() if attrdict['m'] < brightthresh ] )

    starposns = [(e['dec'], e['ra']) for _, e in brightsg.node.items()]

    stargridedges = startogrid(starposns, decs, ras)

    cg = nx.disjoint_union(ng,brightsg)
    ngridnodes = nd * nr
    for i, emap in enumerate(stargridedges):
        n = i + ngridnodes
        epairs = [(n, e) for e in emap]
        cg.add_edges_from(epairs)

    return cg, list(brightsg.nodes)

def kdeweights(sg, brightthresh, h, bfac, comb, peaks, kernel):

    if comb == 'max':
        fn = lambda x, y: np.maximum(x, y)
    else:
        fn = lambda x, y: x + y

    # see https://asterism.org/wp-content/uploads/2019/04/tut35-Magnitudes.pdf
    magfn = lambda x: (10 ** (-x/2.5)) ** bfac   # converts to luminosity then scales using bfac

    xtrain = [(e['dec'], e['ra']) for _, e in sg.node.items() if e['m'] < brightthresh]
    nns = [e['nn_scale'] for _, e in sg.node.items() if e['m'] < brightthresh]
    ms = [e['m_scale'] for _, e in sg.node.items() if e['m'] < brightthresh]
    lks = [e['nn_locw'] for _, e in sg.node.items() if e['m'] < brightthresh]
    if peaks == 'standard':
        weights = [ magfn(m) for m in ms]
        hs = [ h * nn for nn in nns]
    elif peaks == 'rescaled':
        weights = [ magfn(m) * nn ** 2 for (m, nn) in zip(ms,nns) ]
        hs = [ h * nn for nn in nns]
    elif peaks == 'nbstandard':
        weights = [ magfn(m) * lk ** 2 for (m, lk) in zip(ms,lks) ]
        hs = [ h * nn * lk for (nn, lk) in zip(nns, lks)]
    elif peaks == 'nbrescaled':
        weights = [ magfn(m) * (nn ** 2) * (lk**2) for (m, nn, lk) in zip(ms,nns, lks) ]
        hs = [ h * nn * lk for (nn, lk) in zip(nns, lks)]

    allnode = [(e['dec'], e['ra']) for _, e in sg.node.items() ]
    ksum = 0
    for dr, w, h in zip(xtrain, weights, hs):
        tree = sklearn.neighbors.BallTree(np.matrix(dr), metric='haversine')
        ksum = fn(ksum ,w * tree.kernel_density(allnode, h=h, kernel=kernel))

    for n in list(sg.nodes):
        sg.add_node(n, ksum=ksum[n])

    return

def make_rag(sg, thresh):
    rag = graph.rag.RAG()
    for (u, v, e) in sg.edges(data=True):
        w = 0 if sg.nodes[u]['ksum'] > thresh and sg.nodes[v]['ksum'] > thresh else 1
        rag.add_edge(u, v, weight=w)

    for n in rag.nodes():
            rag.node[n].update({'labels': [n]})

    return rag

def weight_boundary(graph, src, dst, n):
    """
    Handle merging of nodes of a region boundary region adjacency graph.

    This function computes the `"weight"` attribute of the edge between `n` and
    the node formed after merging `src` and `dst`.

    Parameters
    ----------
    graph : RAG
        The graph under consideration.
    src, dst : int
        The vertices in `graph` to be merged.
    n : int
        A neighbor of `src` or `dst` or both.

    Returns
    -------
    data : dict
        A dictionary with the "weight" attribute to be
        assigned for the merged node.

    """
    default = {'weight': 1}

    weight_src = graph[src].get(n, default)['weight']
    weight_dst = graph[dst].get(n, default)['weight']

    return {
        'weight': np.min([weight_src, weight_dst])
    }


def merge_boundary(graph, src, dst):
    """Call back called before merging 2 nodes.

    In this case we don't need to do any computation here.
    """
    pass

def labs2clusters(labs):
    _, b = np.unique(labs, return_inverse=True)
    cs = [{}]
    for i in np.arange(np.max(b)):
        cs = np.append(cs, set(np.where(b==i)[0]))
    return cs

def add_nn(stargraph):
    for n in stargraph.nodes():
        ds = [stargraph[u][v]['d'] for (u,v) in stargraph.edges(n)]
        nn = np.deg2rad(np.min(ds)) if ds else None
        stargraph.add_node(n, nn=nn)

def dropgrids(oldc, fixp, snnodes):
    newc = []
    ngrid = fixp['nrgrid'] * fixp['ndgrid']
    for c in oldc:
        c = [snnodes[ci-ngrid] for ci in list(c) if ci >  ngrid]
        if len(c)>1:
            newc.append(set(c))

    return newc

# add ecount
def addecount(sg, f=[], const=0, constval = 3):
    if not f:
        maxw = np.max([e['weight'] for u, v, e in sg.edges(data=True)])
        f = lambda x: 2 + np.floor(x * 28 / maxw)
    if const:
         f = lambda x: constval
    for u, v, e in sg.edges(data=True):
        sg.add_edge(u, v, ecount= f(e['weight']))
    return


def filterpairs(mc, me, sg):
    newmc = []
    newme = []
    for c, e in zip(mc, me):
        lc = list(c)
        if len(c) == 2 and (sg.nodes[lc[0]]['m'] > 3 or sg.nodes[lc[1]]['m'] > 3):
                pass
        else:
            newmc.append(c)
            newme.append(e)
    return newmc, newme

def printperculturescores(outfile, names, scores):
    with open(outfile, "w") as f:
        for i, n in enumerate(names):
            if n.endswith('stellarium'):
                n = n[:-11]
            for j, sc in enumerate(scores[i]):
                f.write(n +  ',' + str(j+1) + ',' + str(sc) + '\n')


# sort list of stars for display in tables

greekd = {
'Alp':'A',
'Bet':'B',
'Gam':'C',
'Del':'D',
'Eps':'E',
'Zet':'F',
'Eta':'G',
'The':'H',
'Iot':'I',
'Kap':'K',
'Lam':'L',
'Mu':'M',
'Nu':'N',
'Xi':'O',
'Omi':'P',
'Pi':'Q',
'Rho':'R',
'Sig':'S',
'Tau':'T',
'Ups':'U',
'Phi':'V',
'Chi':'W',
'Psi':'X',
'Ome':'Y'
}

pattern = re.compile('|'.join(greekd.keys()))
replacegreek = lambda s: pattern.sub(lambda x: greekd[x.group()], s)
uppers = lambda s: ''.join(reversed(re.sub('[^A-Z]', '', 'Z' + s)))
initialnums= lambda s: re.match('^\d+', '0'+s).group()

def sortnamelist(names):
    triples = [(uppers(replacegreek(s)), initialnums(s), s) for s in names]
    triples = sorted(triples)
    names_sorted = [c for (a,b,c) in triples]
    return(names_sorted)


def stargraph_listnodes(stargraph, outfile, hr2hip):
    with open(outfile, "w") as f:
        ramap = lambda x: np.pi - (x + 1.45 * np.pi) if x < 0.55 * np.pi else np.pi - (x - 0.55 * np.pi)
        for node, e in stargraph.nodes(data=True):
            ehstr = str(e['h'])
            if ehstr in hr2hip:
                hipstr = str(hr2hip[ehstr])
            else:
                hipstr = 'miss'
            f.write(','.join([str(e['i']), ehstr, hipstr, e['name'], str(round(ramap(e['ra']), 3)), str(round(e['dec'], 3)) + '\n']))

def colorgroups(stargraph):
    starcolormap = [('34DelOri', 'crimson'),
                    ('46EpsOri', 'navy'),
                    ('8Bet1Sco', 'darkgreen'),
                    ('35LamSco', 'darkgreen'),
                    ('20EpsSgr', 'darkorange'),
                    ('27GamCas', 'orange'),
                    ('50AlpUMa', 'darkslateblue'),
                    ('33LamUMa', 'darkslateblue'),
                    ('9IotUMa', 'darkslateblue'),
                    ('50AlpCyg', 'red'),
                    ('21EtaCyg', 'red'),
                    ('Alp1Cru', 'red'),
                    ('25EtaTau', 'orangered'),
                    ('87AlpTau', 'goldenrod'),
                    ('9AlpDel', 'olive'),
                    ('66AlpGem', 'mediumblue'),
                    ('27EpsGem', 'mediumblue'),
                    ('24GamGem', 'mediumblue'),
                    ('7DelCrv', 'dodgerblue'),
                    ('32AlpLeo', 'violet'),
                    ('17EpsLeo', 'violet'),
                    ('DelVel', 'mediumpurple'),
                    ('62EtaAqr', 'orchid'),
                    ('76DelAqr', 'orchid'),
                    ('5827', 'steelblue'),
                    ('25RhoBoo', 'coral'),
                    ('3AlpLyr', 'darkorange'),
                    ('10BetLyr', 'darkorange'),
                    ('53BetPeg', 'teal'),
                    ('42ZetPeg', 'teal'),
                    ('21AlpAnd', 'teal'),
                    ('AlpCrA', 'dodgerblue'),
                    ('5AlpCrB', 'steelblue'),
                    ('13AlpAri', 'orangered'),
                    ('Alp1Cen', 'deeppink'),
                    ('7BetUMi', 'lightgreen'),
                    ('11AlpLep', 'deeppink'),
                    ('13GamLep', 'deeppink'),
                    ('33AlpPer', 'royalblue'),
                    ('26BetPer', 'royalblue'),
                    ('45EpsPer', 'royalblue'),
                    ('38OmiPer', 'royalblue'),
                    ('13AlpAur', 'orange'),
                    ('37TheAur', 'orange'),
                    ('53AlpAql', 'orange'),
                    ('23BetDra', 'paleturquoise'),
                    ('14EtaDra', 'paleturquoise'),
                    ('47DelCnc', 'darkseagreen'),
                    ('42AlpCom', 'rosybrown'),
                    ('6Alp2Cap', 'steelblue'),
                    ('40GamCap', 'steelblue'),
                    ('68DelLeo', 'orange'),
                    ('25DelCMa', 'powderblue'),
                    ('1Pi3Ori', 'mediumblue'),
                    ('16AlpBoo', 'darkgoldenrod'),
                    ('26EpsSco', 'darkgreen'),
                    ('63EpsDra', 'violet'),
                    ('44ChiDra', 'violet'),
                    ('BetGru', 'steelblue'),
                    ('4BetTri', 'dodgerblue'),
                    ('86MuHer', 'navy'),
                    ('28BetSer', 'darkseagreen'),
                    ('27BetHer', 'mediumpurple'),
                    ('39OmiSgr', 'royalblue'),
                    ('AlpCol', 'orange'),
                    ('AlpMus', 'orange'),
                    ('ZetAra', 'royalblue'),
                    ('11EpsHya', 'navy'),
                    ('GamCen', 'royalblue'),
                    ('AlpLup', 'mediumslateblue')
                    ]
    clustercolors_g = collections.OrderedDict()
    graphnodesout = '../output/data/all_stargraph.txt'
    df = pd.read_csv(graphnodesout, names=['ouri', 'h', 'hip', 'name', 'ra', 'dec'])

    for (name, ccolor) in starcolormap:
         ourind =  df[df['name'] == name].ouri.values[0]
         clustercolors_g[ourind] =  colors.to_rgba(ccolor)

    return(clustercolors_g)

