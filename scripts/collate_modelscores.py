import numpy as np
import pandas as pd
import pickle

modelnames = [ 'model_knn', 'model_knn_thresh', 'model_n', 'model_onebig', 'model_singleton',
              'model_n_nobright', 'model_n_nodistance', 'model_n_nodilation']
scorefns= ['F', 'F10', 'adj_rand']

codemodelnames = ['model_code_' + str(i) for i in np.arange(15)]
codemodelnames = ['model_code_13']  # best performing CODE model

outdir = '../output/results/'
collatedfile= '../output/results/allmodelscores.csv'

allscores = []
allbestpars = []

allscoresopt = []
allbestparsopt = []
for scoref in scorefns:
    for modelname in modelnames:
        outfile = outdir + modelname + '_' + scoref + '.pickle'
        with open(outfile, 'rb') as f_outfile:
            outd = pickle.load(f_outfile)
            bestscore = np.max(outd['scores'])
            bestpars = outd['pvals'][np.argmax(outd['scores'])]
            allscores.append( (modelname, scoref, bestscore))
            allbestpars.append( bestpars )

    for modelname in codemodelnames:
        outfile = outdir + modelname + '_' + scoref + '_opt.pickle'
        with open(outfile, 'rb') as f_outfile:
            outd = pickle.load(f_outfile)
            bestscore = -outd['fun']
            bestpars = outd['bestpars']
            allscoresopt.append( (modelname, scoref, bestscore))
            allbestparsopt.append( bestpars )


allscores = allscores + allscoresopt

df = pd.DataFrame(allscores, columns=['model', 'scorefn', 'score'])
df.to_csv(collatedfile, index=False)

allbestpars= pd.DataFrame(allbestpars)

