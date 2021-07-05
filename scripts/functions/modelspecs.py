import numpy as np
import functions.models as models
import functions.parmaps as parmaps
import functions.starhelpers as sh
import functions.readgroups as readgroups


def modelspecs(starfile, attestedgroupsfile):
    pm, invpm = parmaps.parmaps()
    fixp = {'b':[invpm['bmap'](3.5), invpm['bmap'](4.0), invpm['bmap'](4.5)], 'd':invpm['dmap'](60.0),
            'n':invpm['nmap'](320),'r':invpm['rmap'](3.0), 'r_tq':invpm['rmap'](3.0), 'nrgrid':120, 'ndgrid':60, 'h':0.05,
            'bf':0.1, 'comb':'max', 'mscaletyp':'orig', 'nnscaletyp':'lval', 'kernel':'exponential',
            'peaks':'standard'}

    pgrids = {}
    pgrids['b'] = [invpm['bmap'](i) for i in np.arange(2.9, 4.2, 0.05)]
    # scalar version of 'b'
    pgrids['b2'] = [invpm['bmap'](i) for i in np.arange(1.5, 4.5, 0.1)]
    pgrids['n'] = [invpm['nmap'](i) for i in np.arange(1, 2001)] # for results
    pgrids['r'] = [invpm['rmap'](i) for i in np.append(np.array([np.inf, 0, 1, 2]), np.arange(3, 5, 0.1))]
    pgrids['r_tq'] = [invpm['rmap'](i) for i in np.array([np.inf, 0,1,2,3,4,5,7,9,11,13,15,17,19,21,23,25,27,30,35,40,80])]
    pgrids['d'] = [invpm['dmap'](i) for i in np.arange(30, 80, 5)]
    pgrids['t'] = np.arange(0.76, 0.86, 0.02)
    # old pgrids['bf'] = [0.001, 0.01, 0.1, 0.5]
    pgrids['bf'] = [0, 0.1, 0.5, 1, 2]
    pgrids['h'] = [0.01, 0.05, 0.1, 0.5]
    # could set max k based on size of largest system
    pgrids['k'] = [invpm['kmap'](i) for i in np.arange(1, 300)]

    # add mf_lr, df_lr here
    brightthresh = [3.5, 4.0, 4.5]
    stargraph = sh.makestargraph(starfile, brightthresh)
    cultures = readgroups.cultures(attestedgroupsfile)
    hc = [c for culture in cultures for c in culture]
    he = sh.mstedges(hc, stargraph)
    mf_lr, df_lr, *_  = sh.empirical_ds(hc, he, stargraph, 4.5)
    fixp['mf_lr']=mf_lr
    fixp['df_lr']=df_lr

    modelfs = {}
    parnamess = {}
    guesss = {}
    gridss = {}

    # ----------------------------------------------------------------------
    # Two parameter model: brightness, nedges
    modelfs['model_bn'] = models.model_bn
    parnamess['model_bn'] = ['b', 'n']
    gridss['model_bn'] = [np.arange(0.7, 0.72, 0.01), np.arange(0.5, 0.52, 0.01)]

    # ----------------------------------------------------------------------
    # Two parameter model: brightness, nedges with pp
    modelfs['model_bn_pp'] = models.model_bn_pp
    parnamess['model_bn_pp'] = ['b', 'n']
    guesss['model_bn_pp'] = [0.77, 0.62]

    # ----------------------------------------------------------------------
    # Three parameter model: brightness, rho, nedges
    modelfs['model_brn'] = models.model_brn
    parnamess['model_brn'] = ['b', 'r', 'n']
    guesss['model_brn'] = [1.08, 1.369, 0.905]

    # ----------------------------------------------------------------------
    # Three parameter model: brightness, dnb, nedges
    modelfs['model_bdn'] = models.model_bdn
    parnamess['model_bdn'] = ['b', 'd', 'n']
    guesss['model_bdn'] = [0.75, 1.33, 0.56]

    # ----------------------------------------------------------------------
    # Four parameter model: brightness, dnb, rho, nedges. Dilation on distance only before combination
    modelfs['model_bdrn_donly'] = models.model_bdrn_donly
    parnamess['model_bdrn_donly'] = ['b', 'd', 'r', 'n']
    guesss['model_bdrn_donly'] = [0.845, 1.6, 0.805, 0.85]

    # ----------------------------------------------------------------------
    # Four parameter model: brightness, dnb, rho, nedges. Dilation on brightness only before combination
    modelfs['model_bdrn_bonly'] = models.model_bdrn_bonly
    parnamess['model_bdrn_bonly'] = ['b', 'd', 'r', 'n']
    guesss['model_bdrn_bonly'] = [0.8844, 2.406, 0.476, 0.85]

    # ----------------------------------------------------------------------
    # Four parameter model: brightness, dnb, rho, nedges. Dilation on both before combination
    modelfs['model_bdrn_both'] = models.model_bdrn_both
    parnamess['model_bdrn_both'] = ['b', 'd', 'r', 'n']
    guesss['model_bdrn_both'] = [0.845, 2.0, 0.63, 0.84]
    # ----------------------------------------------------------------------
    # Four parameter model: brightness, dnb, rho, nedges. Dilation after combination
    modelfs['model_brdn'] = models.model_brdn
    parnamess['model_brdn'] = ['b', 'd', 'r', 'n']
    guesss['model_brdn'] = [0.845, 1.685265, 0.34058, 0.82487]
    # ----------------------------------------------------------------------
    # Three parameter model: brightness, rho, nedges
    modelfs['model_brdn_pp'] = models.model_brdn_pp
    parnamess['model_brdn_pp'] = ['b', 'r', 'n']
    guesss['model_brdn_pp'] = [1.08, 1.369, 0.905]
    # ----------------------------------------------------------------------
    # Two parameter model: rho, nedges (brightness, dnb treated as constants)
    modelfs['model_rn'] = models.model_rn
    parnamess['model_rn'] = ['r', 'n']
    guesss['model_rn'] = [ invpm['rmap'](1.61), invpm['nmap'](231) ] 

    # Two parameter model: rho, nedges (brightness, dnb treated as constants)
    modelfs['model_rn_nod'] = models.model_rn_nod
    parnamess['model_rn_nod'] = ['r', 'n']
    guesss['model_rn_nod'] = [0.776, 0.913] # for 3.5 r = 1.76

    # One parameter model: rho (brightness, dnb, nedges treated as constants)
    modelfs['model_r'] = models.model_r
    parnamess['model_r'] = ['r']
    guesss['model_r'] = [invpm['rmap'](3.0)] 

    # One parameter model: rho (brightness, dnb, nedges treated as constants)
    modelfs['model_r_nodilation'] = models.model_r_nodilation
    parnamess['model_r_nodilation'] = ['r']
    guesss['model_r_nodilation'] = [invpm['rmap'](3.0)] # produces better looking result

    # One parameter model: nedges (brightness, dnb, rhotreated as constants). This is the model called the GC model in the text
    modelfs['model_n'] = models.model_n
    parnamess['model_n'] = ['n']
    guesss['model_n'] = [invpm['nmap'](320)]

    # One parameter model: nedges (rho = 0)
    modelfs['model_n_nobright'] = models.model_n_nobright
    parnamess['model_n_nobright'] = ['n']
    guesss['model_n_nobright'] = [invpm['nmap'](320)]

    # One parameter model: nedges (rho = 0)
    modelfs['model_n_nodistance'] = models.model_n_nodistance
    parnamess['model_n_nodistance'] = ['n']
    guesss['model_n_nodistance'] = [invpm['nmap'](320)]

    # One parameter model: nedges (no dilaton)
    modelfs['model_n_nodilation'] = models.model_n_nodilation
    parnamess['model_n_nodilation'] = ['n']
    guesss['model_n_nodilation'] = [invpm['nmap'](320)]

    # One parameter model: rho
    modelfs['model_r_triple'] = models.model_r_triple
    parnamess['model_r_triple'] = ['r_tq']
    guesss['model_r_triple'] = [invpm['rmap'](1.0)]

    # One parameter model: rho
    modelfs['model_r_quad'] = models.model_r_quad
    parnamess['model_r_quad'] = ['r_tq']
    guesss['model_r_quad'] = [invpm['rmap'](1.0)]

    # One parameter model: nedges 
    modelfs['model_elder_md'] = models.model_elder_md
    parnamess['model_elder_md'] = ['n']
    guesss['model_elder_md'] = [invpm['nmap'](320)]


    # CODE model: nedges (brightness, dnb, rhotreated as constants)
    modelfs['model_code'] = models.model_code
    parnamess['model_code'] = ['t', 'bf', 'h']
    guesss['model_code'] =  [0.70, 0.01, 0.05]

    # CODE model: nedges (brightness, dnb, rhotreated as constants)
    modelfs['model_code_1p'] = models.model_code_1p
    parnamess['model_code_1p'] = ['t']
    guesss['model_code_1p'] =  [0.79] # 3.5

    for suff in np.arange(16):
        mname = 'model_code_' + str(suff)
        modelfs[mname]=getattr(models, mname)
        parnamess[mname] = ['t', 'bf', 'h']
        guesss[mname] =  [0.82, 2, 0.1]

    # kNN model
    modelfs['model_knn'] = models.model_knn
    parnamess['model_knn'] = ['k']
    guesss['model_knn'] =  [300] 

    # kNN threshold model
    modelfs['model_knn_thresh'] = models.model_knn_thresh
    parnamess['model_knn_thresh'] = ['k', 'b2']
    guesss['model_knn_thresh'] =  [115, invpm['b2map'](3.1)]

    # singleton model
    modelfs['model_singleton'] = models.model_singleton
    parnamess['model_singleton'] = ['b2']
    guesss['model_singleton'] =  [invpm['b2map'](4.5)]

    # onebig model
    modelfs['model_onebig'] = models.model_onebig
    parnamess['model_onebig'] = ['b2']
    guesss['model_onebig'] =  [invpm['b2map'](4.5)]

    return fixp, pgrids, modelfs, parnamess, guesss, gridss


def printnames(mname):
    pn = {'model_n':'GC',
          'model_n_nobright':'GC (no brightness)',
          'model_n_nodistance':'GC (no distance)',
          'model_n_nodilation':'GC (no scaling)',
          }
    if mname in pn:
        pname = pn[mname]
    else:
        pname = mname

    return pname

