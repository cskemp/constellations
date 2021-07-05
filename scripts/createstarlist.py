import numpy as np
import ephem
import collections
import json
import functions.neighborhood as neighborhood
import functions.starhelpers as sh

# follow chain of substitutions
def finalmap(i, subdict):
    if not(i in subdict):
        return i
    else:
        return finalmap(subdict[i], subdict)

with open("../output/data/hrnames.txt", "rb") as f:
    hrnames = json.load(f)

hnamed = dict()

for (hn, name) in hrnames:
    hnamed[hn] = name

# Read data from the Yale Bright Star Catalogue
# Downloaded from http://tdc-www.harvard.edu/catalogs/bsc5.html bsc5.dat.gz [ASCII catalog]

ybspath = "../data/bsc5.dat"
ybsdict = collections.defaultdict(dict)

with open(ybspath) as ybs:
    for line in ybs:
        linehrid = line[0:4]
        linemag= line[102:107]
        # J2000 
        lineRAh = line[75:77]
        lineRAm = line[77:79]
        lineRAs = line[79:83]
        lineDEPM = line[83:84]
        lineDEd = line[84:86]
        lineDEm = line[86:88]
        lineDEs = line[88:90]
        try: 
            hrid = int(linehrid)
            mag = float(linemag)
            RAh = float(lineRAh)
            RAm = float(lineRAm)
            RAs = float(lineRAs)
            RA = (RAh + RAm / 60.0 + RAs / 3600) * 15
            DEd = float(lineDEd)
            DEm = float(lineDEm)
            DEs = float(lineDEs)
            Dec = DEd + DEm / 60.0 + DEs / 3600
            if lineDEPM == "-":
                Dec *= -1
            star = ephem.FixedBody()
            # set these attributes in radians
            star._ra = np.deg2rad(RA)
            star._dec = np.deg2rad(Dec)
            star._epoch= ephem.Date(ephem.J2000)
            star.compute()
            ybsdict[hrid] = {'ra':float(star.a_ra), 'dec':float(star.a_dec), 'mag':mag, 'cids':[]}
        except ValueError:
            print("HRID " + str(hrid) + " missing information")

# manual alterations

ybsdict[5958]['mag'] = 10  # recurring nova: listed in YBS with mag = 2.0
# Mira: a variable star
# Wikipedia: "In the particular case of Mira, its increases in brightness
# take it up to about magnitude 3.5 on average" 
ybsdict[681]['mag'] = 3.5  # variable star

# make input for grouping models
Star = collections.namedtuple('Star', 'hrid ra dec mag')
stars = [Star(*[key, v['ra'], v['dec'], v['mag']]) for key, v in ybsdict.items()]
# include only stars visible with naked eye
minmag = 6.5
stars = [(h,r,d,m) for (h,r,d,m) in stars if  m <= minmag]
# add starnames
stars = [(h,r,d,m, hnamed[str(h)]) if str(h) in hnamed else (h,r,d,m,'XXX') for (h,r,d,m) in stars]

removesubs = []

# combine stars with same position (e.g. two stars in Gam Vir)
posdict = collections.defaultdict(dict)
for i, (h,r,d,m, name) in enumerate(stars):
    starpos = (round(r,5), round(d,5))
    if starpos in posdict:
        m1 = posdict[starpos][1]
        m2 = m
        combinedmag = sh.magsum(m1, m2)
        test1 = stars.pop(i)  # pop current item 
        prevind = posdict[starpos][0]
        (oldh, oldr, oldd, oldm, oldn) = stars.pop(prevind)  # pop prev item
        stars.insert(prevind, (oldh, oldr, oldd, combinedmag, oldn) )
        posdict[starpos] = (prevind, combinedmag)
        print("combining " + str(h) + " " + str(oldh)) 
        removesubs.append( (h, oldh) )
    else:
        posdict[starpos] = (i, m)
        
# drop items that are very close

nb = neighborhood.makeneighborhood(stars, [minmag])
G, _ = sh.makedistgraph(nb, stars)

remove = []
faintmag = 800
newmagdict = collections.defaultdict(dict)
for v1,v2 in G.edges():
    if G[v1][v2]['d'] < 0.1:
        m1= G.node[v1]['m']
        m2= G.node[v2]['m']
        # combined magnitudes for close stars
        combinedmag = sh.magsum(m1, m2)
        if m1 > m2:
            remove.append(v1)
            G.node[v2]['m'] = combinedmag
            G.node[v1]['m'] = faintmag
            print("removing: " + str(G.node[v1]['h']) + " for " + str(G.node[v2]['h']) )  
            removesubs.append( (G.node[v1]['h'], G.node[v2]['h']) )
            newmagdict[G.node[v2]['h']] = combinedmag
        else:
            remove.append(v2)
            G.node[v1]['m'] = combinedmag
            G.node[v2]['m'] = faintmag
            print("removing: " + str(G.node[v2]['h']) + " for " + str(G.node[v1]['h']) ) 
            removesubs.append( (G.node[v2]['h'], G.node[v1]['h']) )
            newmagdict[G.node[v1]['h']] = combinedmag

include = np.setdiff1d(range(len(stars)), remove)
oldstars = stars
stars = [oldstars[i] for i in include]

for i, (h,r,d,m, name) in enumerate(stars):
    if h in newmagdict:
        (oldh, oldr, oldd, oldm, oldn) = stars.pop(i)
        stars.insert(i, (oldh, oldr, oldd, newmagdict[h], oldn) )

with open("../output/data/stars.txt", "w") as f:
    json.dump(stars, f)

with open("../output/data/stars.csv", "w") as f:
    for s in stars:
        f.write(','.join(map(str, s)) + '\n')


typeconvert= lambda x: x.item() if isinstance(x, np.int32) else x

# convert to int so we can export to json format
removesubs = [(typeconvert(a), typeconvert(b)) for (a,b) in removesubs]

removed = dict()
for (a,b) in removesubs:
    removed[a] = b

finalmaps = [(a, finalmap(a, removed)) for (a,b) in removesubs]

with open("../output/data/removals.txt", "w") as f:
    json.dump(finalmaps, f)
