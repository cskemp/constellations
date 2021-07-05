import numpy as np
import sklearn.neighbors
import scipy.spatial
from scipy.sparse import lil_matrix
from angles import sep, bear

def angsep(s1, s2):
    """ Determine separation in degrees between two celestial objects 
        arguments are RA and Dec in radians. 
    """
    # calculate scalar product for determination
    # of angular separation

    ra1rad = s1[0]
    dec1rad = s1[1]

    ra2rad = s2[0]
    dec2rad = s2[1]

    return np.rad2deg(sep(ra1rad, dec1rad, ra2rad, dec2rad))


def polar2cartesian(rd):
  r = rd[0]
  d = rd[1]
  x = np.cos(d) * np.cos(r)  
  y = np.cos(d) * np.sin(r)
  z = np.sin(d) 
  return (x,y,z)


def delaunaytriangulate(stars):
  pcs = [(r,d) for (h,r,d,m,n) in stars]
  rcs = list(map( polar2cartesian, pcs ))
  rcs = np.array(rcs)
  tri = scipy.spatial.ConvexHull(rcs)
  return tri

# make delaunay neighborhood
def makeneighborhood(stars, thresh):
    if len(thresh) > 1:
        first = makeneighborhood(stars, [thresh[0]])
        rest = makeneighborhood(stars, thresh[1:])
        return scipy.sparse.lil_matrix.maximum(first, rest)
    else:
        thresh = thresh[0]

    bstars= [ (h,r,d,m,i)  for (i, (h,r,d,m,n)) in enumerate(stars) if m <= thresh] 

    n = len(stars)
    rcs = delaunaytriangulate(bstars)
    nb = lil_matrix( (n, n) )
    for s in rcs.simplices:
        i1 = bstars[s[0]][4]
        i2 = bstars[s[1]][4]
        i3 = bstars[s[2]][4]
        nb[i1,i2] = 1
        nb[i2,i3] = 1
        nb[i1,i3] = 1

    return nb

# minimal spanning tree neighborhood
def makeneighborhood_mst(stars, thresh):
    if len(thresh) > 1:
        first = makeneighborhood_mst(stars, [thresh[0]])
        rest = makeneighborhood_mst(stars, thresh[1:])
        return scipy.sparse.lil_matrix.maximum(first, rest)
    else:
        thresh = thresh[0]

    delaunay = makeneighborhood(stars, [thresh])

    [rs, cs] = delaunay.nonzero()
    for (r,c) in zip(rs, cs):
        s1 = (stars[r][1], stars[r][2]) 
        s2 = (stars[c][1], stars[c][2]) 
        delaunay[r,c] = angsep(s1, s2)

    mst = scipy.sparse.csgraph.minimum_spanning_tree(delaunay)

    return mst

def scoreneighborhood(nb, clusters, edges, stars):
    ncluster = len(clusters)
    edgeinnb = 0
    edgeoutnb = 0
    for i in range(1, ncluster):
        clusterstars = clusters[i] 
        ris, cis = edges[i].nonzero()
        for s1, s2 in zip(ris, cis):
            if nb[clusterstars[s1], clusterstars[s2]] or nb[clusterstars[s2], clusterstars[s1]]:
                edgeinnb += 1
            else:
                edgeoutnb += 1

    return edgeinnb, edgeoutnb

def makecloseedges_nb(stars, thresh, maxd):
    thresh = np.max(thresh)
    bstars= [ (h,r,d,m,c, i)  for (i, (h,r,d,m,c)) in enumerate(stars) if m <= thresh] 
    n = len(stars)
    nb = lil_matrix( (n, n) )

    brightn = len(bstars)
    rd = np.matrix([ (d,r) for (h,r,d,m,c,i) in bstars ])
    tree = sklearn.neighbors.BallTree(rd, metric='haversine')

    i = 0
    nns = tree.query_radius(rd, r=np.deg2rad(maxd))
    for i in range(brightn):
        iorig = bstars[i][5]
        for nni in nns[i]:
            nniorig = bstars[nni][5]
            if i != nni:
                nb[iorig,nniorig] = 1

    return nb

def sphericalangle(triple, sg):

    rds = [ (sg.nodes[i]['ra'], sg.nodes[i]['dec']) for i in triple] 

    a = np.deg2rad(angsep(rds[0], rds[1]))
    b = np.deg2rad(angsep(rds[1], rds[2]))
    c = np.deg2rad(angsep(rds[0], rds[2]))

    # based on https://en.wikipedia.org/wiki/Spherical_law_of_cosines

    C1 = np.arccos( (np.cos(c) - np.cos(a)*np.cos(b)) / (np.sin(a)*np.sin(b) ) )

    if np.isnan(C1):
        print("isnan in spherical angle: ",  a, b, c)

#    turn1 = turndir(rds[0], rds[1], rds[2])
#    if turn1 > 0:
#        C1 = 2 * np.pi - C1

    return C1


def tripledistance(triple, sg):

    rds = [ (sg.nodes[i]['ra'], sg.nodes[i]['dec']) for i in triple]

    a = np.deg2rad(angsep(rds[0], rds[1]))
    b = np.deg2rad(angsep(rds[1], rds[2]))

    if a < b:
        ratio = b/a
    else:
        ratio = a/b

    return ratio


# if walking from s1 to s2 and then on to s3, does one turn right or left at s2?

def turndir(s1, s2, s3):
    ra1rad = s1[0]
    dec1rad = s1[1]

    ra2rad = s2[0]
    dec2rad = s2[1]

    ra3rad = s3[0]
    dec3rad = s3[1]

    b2 = bear(ra1rad, dec1rad, ra2rad, dec2rad)
    b3 = bear(ra1rad, dec1rad, ra3rad, dec3rad)

    tsign = 0
    if (b2 <= np.pi):   # no wrap around issgues
        if (b2 <= b3 and b3 <= b2 + np.pi):
            tsign = 1
    else:
        tsign = 1
        if (b3 > b2 - np.pi and b3 < b2):
            tsign = 0

    return tsign


def goodc(quad, sg):
    rds = [ (sg.nodes[i]['ra'], sg.nodes[i]['dec']) for i in quad]

    b_01 = bear(rds[0][0], rds[0][1], rds[1][0], rds[1][1])
    b_12 = bear(rds[1][0], rds[1][1], rds[2][0], rds[2][1])
    b_23 = bear(rds[2][0], rds[2][1], rds[3][0], rds[3][1])

    d1 = b_12 - b_01
    d2 = b_23 - b_12

    alld2 = [d2 - 2*np.pi, d2, d2 + 2*np.pi]
    diffs = [d1 - d2prime for  d2prime in alld2]

    return np.min(np.abs(diffs))

def goodcorig(quad, sg):

    ang1 = sphericalangle(quad[0:3], sg)
    ang2 = sphericalangle(quad[1:4], sg)

    return np.abs(ang1 - ang2)

