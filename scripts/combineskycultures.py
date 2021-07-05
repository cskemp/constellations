import json
import os
import functions.starhelpers as sh

# combine attested groups from a number of cultures

infiles = ['../output/data/boorong_processed.txt', 
'../output/data/india_processed.txt',
'../output/data/indomalay_processed.txt',
'../output/data/marshall_processed.txt',
'../output/data/pacariqtambo_processed.txt',
'../output/data/anutan_stellarium_processed.txt',
'../output/data/arabic_moon_stations_stellarium_processed.txt',
'../output/data/belarusian_stellarium_processed.txt',
'../output/data/chinese_stellarium_processed.txt',
'../output/data/dakota_stellarium_processed.txt',
'../output/data/egyptian_stellarium_processed.txt',
'../output/data/inuit_stellarium_processed.txt',
'../output/data/lokono_stellarium_processed.txt',
'../output/data/macedonian_stellarium_processed.txt',
'../output/data/maori_stellarium_processed.txt',
'../output/data/mulapin_stellarium_processed.txt',
'../output/data/navajo_stellarium_processed.txt',
'../output/data/norse_stellarium_processed.txt',
'../output/data/ojibwe_stellarium_processed.txt',
'../output/data/romanian_stellarium_processed.txt',
'../output/data/sami_stellarium_processed.txt',
'../output/data/siberian_stellarium_processed.txt',
'../output/data/tongan_stellarium_processed.txt',
'../output/data/tukano_stellarium_processed.txt',
'../output/data/tupi_stellarium_processed.txt',
'../output/data/vanuatu_netwar_stellarium_processed.txt',
'../output/data/western_stellarium_processed.txt']

starfile = "../output/data/stars.txt"

with open(starfile, "rb") as f:
   stars = json.load(f)

clists = [[] for f in infiles]
cnames= ['' for f in infiles]

for filei, infile in enumerate(infiles):
    clist = clists[filei]
    cname = os.path.basename(infile)
                         # strip suffix
    cnames[filei]= cname[0:-14]
    with open(infile, "r") as f:
        for line in f:
            ids = line.strip()
            # strip [ and ] at start and end
            ids = ids[1:-1]
            ids = list(map(int, ids.split(',')))
            if len(ids) > 0:
                clist.append(ids)

# all human data
cout = '../output/data/allhuman.txt'

hd = (cnames, clists)

with open(cout, 'w') as f:
        json.dump(hd, f)

invhd = sh.invhdict(starfile)

# all human data with our ids
cout = '../output/data/allhuman_ourids.txt'
allmags = []
for ii, cll in enumerate(clists):
    for j, cl in enumerate(cll):
        cll[j] = [invhd[s] for s in cl if s in invhd]
        allmags = allmags + [(cnames[ii], stars[invhd[s]][3]) for s in cl if s in invhd]

hd = (cnames, clists)
with open(cout, 'w') as f:
        json.dump(hd, f)

cout = '../output/data/allhuman_mags.txt'
with open(cout, 'w') as f:
    for pr in allmags:
        f.write(','.join(map(str, pr))+ '\n')

# just attested asterisms
cout = '../output/data/allhuman_opt_ourids.txt'
with open(cout, 'w') as f:
    json.dump(clists, f)

stargraph = sh.makestargraph(starfile, [3.5, 4.0, 4.5])

clists_bright= [ list(map(list, sh.brightclusters(clist, stargraph, 4.5))) for clist in clists]

# all human data with our ids
cout_bright = '../output/data/allhuman_bright_ourids.txt'
hd_bright = (cnames, clists_bright)
with open(cout_bright, 'w') as f:
        json.dump(hd_bright, f)

# just attested asterisms
cout_bright = '../output/data/allhuman_opt_bright_ourids.txt'
with open(cout_bright, 'w') as f:
    json.dump(clists_bright, f)

# all nodes in graph with our ids, hr ids, hip ids, names
graphnodesout = '../output/data/all_stargraph.txt'

with open('../output/data/hr2hip.txt', 'r') as f:
    hr2hip  = dict(json.load(f))

sh.stargraph_listnodes(stargraph, graphnodesout, hr2hip)

