import numpy as np
import networkx as nx
from spherecluster import SphericalKMeans
from networkx.algorithms.components.connected import connected_components
from skimage.future import graph
import functions.starhelpers as sh
import functions.parmaps as parmaps
from types import FunctionType
from copy import copy, deepcopy
import itertools

def copy_function(fn):
    return FunctionType(copy(fn.__code__), copy(fn.__globals__), name=fname,
                        argdefs=copy(fn.__defaults__),
                        closure=copy(fn.__closure__))

pm, invpm = parmaps.parmaps()

def rho2weights(rho):
    if np.isinf(rho):
        wb = 1
        wd = 0
    else:
        wb = rho / (rho + 1)
        wd = 1 / (rho + 1)
    return wb, wd

def topn_clusters(stargraph, nedges, alledge=0):
    # reverse sort
    alledgedata = -np.sort([-d['weight'] for (_, _, d) in stargraph.edges.data()])
    if nedges > len(alledgedata):
        nedges = len(alledgedata)

    if alledge:
        newsg = stargraph.copy()
    else:
        newsg = stargraph
    sh.dropsmalledges(newsg, 'weight', alledgedata[nedges-1], alledge=alledge)
    mcs = list(connected_components(newsg))
    mcs = [mc for mc in mcs if len(mc) > 1]
    return mcs, stargraph

# two parameter model: brightness, nedges
def model_bn(starfile, pars):
    brightthresh, nedges = pm['bmap'](pars[0]), pm['nmap'](pars[1])
    actualpars = (brightthresh, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    sh.setweight(stargraph, 'dsim')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars


# two parameter model: brightness, nedges with pp dilation
def model_bn_pp(starfile, pars):
    brightthresh, nedges = pm['bmap'](pars[0]), pm['nmap'](pars[1])
    actualpars = (brightthresh, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    sh.dilateedges_pp(stargraph, 'dsim', 'dsimpp')
    sh.setweight(stargraph, 'dsimpp')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars


# three parameter model: brightness, rho, nedges
def model_brn(starfile, pars): 
    brightthresh, rho, nedges, alledge = pm['bmap'](pars[0]), pm['rmap'](pars[1]), pm['nmap'](pars[2]), pars[3]
    actualpars = (brightthresh, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    #sh.standardizemag(stargraph)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)
    sh.combineedgeatts(stargraph, ['dsim', 'maxmagsim'], 'dmagprod', ctype='prod')
    sh.setweight(stargraph, 'dmagprod')
    topnc = topn_clusters(stargraph, nedges, alledge=alledge)
    return topnc[0], topnc[1], actualpars

# three parameter model: brightness, dnb, nedges
def model_bdn(starfile, pars): 
    brightthresh, dnb, nedges = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['nmap'](pars[2])
    actualpars = (brightthresh, dnb, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilateedges(stargraph, dnbgraph, 'dsim', 'dilatedsim')
    sh.setweight(stargraph, 'dilatedsim')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars

# four parameter model: brightness, dnb, rho, nedges
def model_bdrn_donly(starfile, pars):
    brightthresh, dnb, rho, nedges = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['rmap'](pars[2]), pm['nmap'](pars[3])
    actualpars = (brightthresh, dnb, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilateedges(stargraph, dnbgraph, 'dsim', 'dilatedsim')
    sh.dist2sim(stargraph, 'maxmag', wb)
    sh.combineedgeatts(stargraph, ['dilatedsim', 'maxmagsim'], 'dmagprod', ctype='prod')
    sh.setweight(stargraph, 'dmagprod')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars

# four parameter model: brightness, dnb, rho, nedges
def model_bdrn_bonly(starfile, pars):
    brightthresh, dnb, rho, nedges = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['rmap'](pars[2]), pm['nmap'](pars[3])
    actualpars = (brightthresh, dnb, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilateedges(stargraph, dnbgraph, 'maxmagsim', 'dilatemaxmagsim')
    sh.combineedgeatts(stargraph, ['dsim', 'dilatemaxmagsim'], 'dmagprod', ctype='prod')
    sh.setweight(stargraph, 'dmagprod')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars

# four parameter model: brightness, dnb, rho, nedges
def model_bdrn_both(starfile, pars, fp=[]):
    brightthresh, dnb, rho, nedges, alledge = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['rmap'](pars[2]), pm['nmap'](pars[3]), pars[4]
    actualpars = (brightthresh, dnb, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilateedges(stargraph, dnbgraph, 'maxmagsim', 'dilatemaxmagsim')
    sh.dilateedges(stargraph, dnbgraph, 'dsim', 'dilatedsim')
    sh.combineedgeatts(stargraph, ['dilatedsim', 'dilatemaxmagsim'], 'dmagprod', ctype='prod')
    sh.setweight(stargraph, 'dmagprod')
    topnc = topn_clusters(stargraph, nedges, alledge=alledge)
    return topnc[0], topnc[1], actualpars

def model_bdrn_both_tripquad(starfile, pars, fp=[]):
    brightthresh, dnb, rho, nedges, rho_tripquad, tripflag, alledge = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['rmap'](pars[2]), pm['nmap'](pars[3]), pm['rmap'](pars[4]), pars[5], pars[6]
    actualpars = (brightthresh, dnb, rho, nedges, rho_tripquad)

    # this is same as model_bdrn_both
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)

    if 1:  # dilation
        dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
        sh.dilateedges(stargraph, dnbgraph, 'maxmagsim', 'dilatemaxmagsim')
        sh.dilateedges(stargraph, dnbgraph, 'dsim', 'dilatedsim')
        sh.combineedgeatts(stargraph, ['dilatedsim', 'dilatemaxmagsim'], 'dmagprod', ctype='prod')
    else: # no dilation
        sh.combineedgeatts(stargraph, ['dsim', 'maxmagsim'], 'dmagprod', ctype='prod')

    # now add triples or quads
    wpairs, wtripquad = rho2weights(rho_tripquad)
    sh.weightatt(stargraph, 'dmagprod', wpairs)
    sh.addtripquad(stargraph, wtripquad, tripflag)

    sh.setweight(stargraph, 'dmagprodtripquad')
    topnc = topn_clusters(stargraph, nedges, alledge=alledge)
    return topnc[0], topnc[1], actualpars

# four parameter model: brightness, rho, dnb, nedges
def model_brdn(starfile, pars, fp=[]):
    brightthresh, dnb, rho, nedges = pm['bmap'](pars[0]), pm['dmap'](pars[1]), pm['rmap'](pars[2]), pm['nmap'](pars[3])
    actualpars = (brightthresh, dnb, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)
    sh.combineedgeatts(stargraph, ['dsim', 'maxmagsim'], 'dmagprod', ctype='prod')
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilateedges(stargraph, dnbgraph, 'dmagprod', 'dilatedmagprod')
    sh.setweight(stargraph, 'dilatedmagprod')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars


# three parameter model: brightness, rho, nedges with pp dilation
def model_brdn_pp(starfile, pars):
    brightthresh, rho, nedges = pm['bmap'](pars[0]), pm['rmap'](pars[1]), pm['nmap'](pars[2])
    actualpars = (brightthresh, rho, nedges)
    stargraph = sh.makestargraph(starfile, brightthresh)
    sh.addmaxmagedges(stargraph)
    sh.scaleedgeweights(stargraph, 'maxmag', 'd')
    wb, wd = rho2weights(rho)
    sh.dist2sim(stargraph, 'd', wd)
    sh.dist2sim(stargraph, 'maxmag', wb)
    sh.combineedgeatts(stargraph, ['dsim', 'maxmagsim'], 'dmagprod', ctype='prod')
    sh.dilateedges_pp(stargraph, 'dmagprod', 'dmagprodpp')
    sh.setweight(stargraph, 'dmagprodpp')
    topnc = topn_clusters(stargraph, nedges)
    return topnc[0], topnc[1], actualpars

# two parameter version of model_bdrn_both
def model_rn(starfile, pars, fp=[]):
    rho, nedges, alledge =  pars[0], pars[1], pars[2]
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_bdrn_both(starfile, (brightthresh, dnb, rho, nedges, alledge))
    return a, b, c

# two parameter version with no dnb
def model_rn_nod(starfile, pars, fp=[]):
    rho, nedges, alledge =  pars[0], pars[1], pars[2]
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_brn(starfile, (brightthresh, rho, nedges, alledge))
    return a, b, c

# one parameter version of model_bdrn_both
def model_r(starfile, pars, fp=[]):
    rho = pars[0]
    alledge = pars[1]
    nedges = fp['n']
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_bdrn_both(starfile, (brightthresh, dnb, rho, nedges, alledge))
    return a, b, c

def model_r_nodilation(starfile, pars, fp=[]):
    rho = pars[0]
    alledge = pars[1]
    nedges = fp['n']
    brightthresh = fp['b']
    a, b, c = model_brn(starfile, (brightthresh, rho, nedges, alledge))
    return a, b, c

# one parameter version of model_bdrn_both: this is the GC model
def model_n(starfile, pars, fp=[]):
    nedges = pars[0]
    alledge = pars[1]
    rho= fp['r']
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_bdrn_both(starfile, (brightthresh, dnb, rho, nedges, alledge))
    return a, b, c

def model_n_nobright(starfile, pars, fp=[]):
    nedges = pars[0]
    alledge = pars[1]
    rho= 0
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_bdrn_both(starfile, (brightthresh, dnb, rho, nedges, alledge))
    return a, b, c

def model_n_nodistance(starfile, pars, fp=[]):
    nedges = pars[0]
    alledge = pars[1]
    rho = np.inf
    brightthresh = fp['b']
    dnb = fp['d']
    a, b, c = model_bdrn_both(starfile, (brightthresh, dnb, rho, nedges, alledge))
    return a, b, c

def model_n_nodilation(starfile, pars, fp=[]):
    nedges = pars[0]
    alledge = pars[1]
    rho= fp['r']
    brightthresh = fp['b']
    a, b, c = model_brn(starfile, (brightthresh, rho, nedges, alledge))
    return a, b, c


# one parameter version of model_bdrn_both_tripquad
def model_r_triple(starfile, pars, fp=[]):
    rho_triple = pars[0]
    alledge = pars[1]
    nedges = fp['n']
    brightthresh = fp['b']
    dnb = fp['d']
    rho= fp['r']
    a, b, c = model_bdrn_both_tripquad(starfile, (brightthresh, dnb, rho, nedges, rho_triple, 1, alledge))
    return a, b, c

# one parameter version of model_bdrn_both_tripquad
def model_r_quad(starfile, pars, fp=[]):
    rho_goodc= pars[0]
    alledge = pars[1]
    nedges = fp['n']
    brightthresh = fp['b']
    dnb = fp['d']
    rho= fp['r']
    a, b, c = model_bdrn_both_tripquad(starfile, (brightthresh, dnb, rho, nedges, rho_goodc, 0, alledge))
    return a, b, c


# Elder-style models

def model_elder_md(starfile, pars, fp=[]):
    nedges = pm['nmap'](pars[0])
    alledge = pars[1]
    brightthresh = pm['bmap'](fp['b'])
    mf_lr = fp['mf_lr']
    df_lr = fp['df_lr']
    actualpars = (brightthresh, nedges)

    stargraph = sh.makestargraph(starfile, brightthresh)

    sh.add_m_lr(mf_lr, stargraph)
    sh.add_min_m_lr_edges(stargraph)
    sh.add_d_lr(df_lr, stargraph)
    sh.combineedgeatts(stargraph, ['min_m_lr', 'd_lr'], 'dm_lr', ctype='prod')
    sh.setweight(stargraph, 'dm_lr')
    topnc = topn_clusters(stargraph, nedges, alledge=alledge)

    return topnc[0], topnc[1], actualpars


# CODE models

def model_code(starfile, pars, fp=[]):
    thresh = pm['tmap'](pars[0])
    bfac = pm['bfmap'](pars[1])
    h = pm['hmap'](pars[2])

    brightthresh = pm['bmap'](fp['b'])
    dnb = pm['dmap'](fp['d'])

    actualpars = (thresh, bfac, h, brightthresh, dnb)
    ndgrid = fp['ndgrid']
    nrgrid = fp['nrgrid']

    stargraph = sh.makestargraph(starfile, brightthresh)
    #stargraph = allstarg.copy()
    #stargraph.remove_nodes_from([n for n, e in allstarg.nodes(data=True) if e['m'] >=brightthresh])
    sh.add_nn(stargraph)
    dnbgraph = sh.makestargraph_closenb(starfile, brightthresh, dnb)
    sh.dilatenodes(stargraph, dnbgraph, 'nn', dtype=fp['nnscaletyp'])
    sh.dilatenodes(stargraph, dnbgraph, 'm', dtype=fp['mscaletyp'])

    sgcopy = nx.create_empty_copy(stargraph)
    brightthresh = np.max(brightthresh)
    cg, starinds  = sh.make_grid(sgcopy, brightthresh, ndgrid, nrgrid)
    # kernel density estimate on grid
    sh.kdeweights(cg, brightthresh, h, bfac, fp['comb'], fp['peaks'], fp['kernel'])

    actualthresh = np.quantile([e['ksum'] for _, e in cg.nodes(data=True)], thresh)
    rag = sh.make_rag(cg, actualthresh)
    labelsin = np.arange(cg.number_of_nodes())
    # rag to find clusters
    labelsout = graph.merge_hierarchical(labelsin, rag, thresh=0.5, rag_copy=False,
                                         in_place_merge=True,
                                         merge_func=sh.merge_boundary,
                                         weight_func=sh.weight_boundary)
    starlabels = labelsout
    clusters = sh.labs2clusters(starlabels)


    # current plotting code needs edge weights
    sh.dist2sim(stargraph, 'd', 1)
    sh.setweight(stargraph, 'dsim')

    c = [ci for ci in clusters if len(ci) > 2 and sum(np.array(list(ci)) > fp['nrgrid'] * fp['ndgrid']) > 1]
    if 1:  # set to 0 for plotting over grid
        # for plotting wrt stargraph
        c = sh.dropgrids(c, fp, starinds)
        cg = stargraph

    return c, cg, actualpars

def model_code_1p(starfile, pars, fp=[]):
    pars = [pars[0], fp['bf'], fp['h']]

    c, cg, actualpars = model_code(starfile, pars, fp)

    return c, cg, actualpars


qpars = [['exponential', 'gaussian'], ['max','sum'], ['lval', 'gval'], ['standard', 'rescaled']]
allqvals = list(itertools.product(*qpars))

for i, qval in enumerate(allqvals):
    fname = 'model_code_' + str(i)
    def _fn(starfile, pars, fp):
        fp['kernel'] = qval[0]
        fp['comb'] = qval[1]
        fp['nnscaletyp'] = qval[2]
        fp['peaks'] = qval[3]
        print( (fp['kernel'], fp['comb'], fp['nnscaletyp'], fp['peaks'] ) )
        c, cg, actualpars = model_code(starfile, pars, fp)
        return c, cg, actualpars

    globals()[fname] = copy_function(_fn)


# kNN model
def cartesian_encoder(coord):
    """
    Input
    -----
        coord : numpy 2darray (size=(N, 2))

    Output
    ------
        out : numpy 2darray (size=(N, 3))
    """

    theta = coord[:, 0]  # lat [radians]
    phi = coord[:, 1]    # lon [radians]

    x = np.cos(phi) * np.cos(theta)
    y = np.sin(phi) * np.cos(theta)
    z = np.sin(theta)

    return np.concatenate([x.reshape(-1, 1), y.reshape(-1, 1), z.reshape(-1, 1)], axis=1)

def model_knn(starfile, pars, fp=[]):
    k = int(pm['kmap'](pars[0]))
    brightthresh = pm['bmap'](fp['b'])
    stargraph = sh.makestargraph(starfile, brightthresh)
    nodes = [n[1] for n in stargraph.nodes(data=True) if n[1]['m'] < np.max(brightthresh)]
    X = np.array([(n['dec'], n['ra']) for n in nodes])
    X_cart = cartesian_encoder(X)
    k = np.min([k, len(X_cart)])
    kmeans_labels = SphericalKMeans(k, n_init=100, max_iter=100).fit_predict(X_cart)
    nodeis = np.arange(len(kmeans_labels))
    mcs = [set([])]
    for j in np.arange(k):
        c = [nodes[i]['i'] for i in nodeis if kmeans_labels[i]== j]
        mcs.append(set(c))

    actualpars = (k, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    sh.setweight(stargraph, 'dsim')

    return mcs, stargraph, actualpars


def model_knn_thresh(starfile, pars, fp=[]):
    fp_tmp = deepcopy(fp)
    fp_tmp['b'] = [pars[1]]
    mcs, stargraph, actualpars = model_knn(starfile, pars, fp_tmp)

    return mcs, stargraph, actualpars


def model_singleton(starfile, pars, fp=[]):
    b = pm['b2map'](pars[0])
    brightthresh = pm['bmap'](fp['b'])
    stargraph = sh.makestargraph(starfile, brightthresh)
    mcs = [set([n[1]['i']]) for n in stargraph.nodes(data=True) if n[1]['m'] <= b]

    actualpars = (b, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    sh.setweight(stargraph, 'dsim')

    return mcs, stargraph, actualpars

def model_onebig(starfile, pars, fp=[]):
    b = pm['b2map'](pars[0])
    brightthresh = pm['bmap'](fp['b'])
    stargraph = sh.makestargraph(starfile, brightthresh)
    mcs = [n[1]['i'] for n in stargraph.nodes(data=True) if n[1]['m'] <= b]
    mcs = [set(mcs)]
    actualpars = (b, brightthresh)
    sh.dist2sim(stargraph, 'd', 1)
    sh.setweight(stargraph, 'dsim')

    return mcs, stargraph, actualpars
