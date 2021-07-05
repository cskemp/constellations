import networkx as nx
import json
import functions.starhelpers as sh

def minnonzero ( a ):
  return min(i for i in a if i > 0)

# returns:
#   clusters[i] = list of vertices
#   edges[i] = sparse graph specifying edges for cluster i

# infile has stars with hrids -- so invhd needed

def group( infile, invhd, hrnames=[]):
    with open(infile) as f:
            lines = f.readlines()

    # strip newline
    ccs = [x.strip() for x in lines] 
    # strip [ and ] at start and end
    ccs = [x[1:-1] for x in ccs] 
    ccs = [map(int, x.split(',')) for x in ccs]

    ncluster = len(ccs)
    
    clusters = [set([]) for i in range(ncluster+1)]

    for i, cc in enumerate(ccs):
        cc = list(cc)
        for j in cc:
            # invhd['2061'] is index of node with hrid 2061
            if j not in invhd:
                if hrnames:
                    print("missing star: " + hrnames[str(j)])
                else:
                    print("missing star: " + str(j))
            else:
                set.add(clusters[i+1], invhd[j])

    return clusters


# infile has stars with our ids -- so invhd not needed

def cultures( infile, hrnames=[]):
    with open(infile, 'r') as f:
        clists = json.load(f)

    clists = [sh.llist2lset(hs) for hs in clists]

    return clists
