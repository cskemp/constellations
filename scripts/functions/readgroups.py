import networkx as nx
import json
import numpy as np
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


def readbayer( infile, stargraph, invbayerd, ename ='d'):
    allcs = [[]]
    alles = [[]]
    with open(infile) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '#':
                continue
            else:
              line = line.strip('\n')
              fields = line.split(',')
              cs = set()
              for i in range(3, len(fields), 2):
                print(fields[i] + '--' + fields[i+1])
                s1 = invbayerd[fields[i]]
                s2 = invbayerd[fields[i+1]]
                cs.add(s1)
                cs.add(s2)
                stargraph.add_edge(s1, s2, weight=ename)
              cs = list(cs)
              #es = nx.minimum_spanning_tree(stargraph.subgraph(cs), weight=ename)
              es = stargraph.subgraph(cs)
              allcs.append(cs)
              alles.append(es)
    return(allcs, alles)
