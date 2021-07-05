import json
import os, glob

def read():
    hip2hrfile = os.path.abspath(os.path.join('..', 'output/data/hip2hr.txt'))
    with open(hip2hrfile) as hrf:
        hip2hr = dict(json.load(hrf))

    removalsfile = os.path.abspath(os.path.join('..', 'output/data/removals.txt'))
    with open(removalsfile) as rf:
        removals = dict(json.load(rf))
    removaldict = dict(removals)

    files = glob.iglob(os.path.abspath(os.path.join('..', 'data/stellarium/*_stellarium.txt')))

    for f in files:
        cname = os.path.basename(f)
        cname = cname[:-15]
        cout = os.path.abspath(os.path.join('..', 'output/data/' + cname + "_stellarium_processed.txt"))

        cs = []
        with open(f) as stellc:
            with open(cout, "w") as stellout:
                for line in stellc:
                    if line[0] == '#':
                        continue
                    else:
                        c = line.split()
                        c = c[2:]
                        c = [hip2hr[int(e)] for e in c if int(e) in hip2hr and hip2hr[int(e)] != '']
                        c = [str(removaldict[int(e)]) if int(e) in removaldict else e for e in c]
                        c = list(set(c))
                        starstr = ','.join(c)
                        if len(c) > 1:
                            stellout.write('[' + starstr + ']\n')
