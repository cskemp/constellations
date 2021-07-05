import os
import json

def merge(lsts):
  sets = [set(lst) for lst in lsts if lst]
  merged = 1
  while merged:
    merged = 0
    results = []
    while sets:
      common, rest = sets[0], sets[1:]
      sets = []
      for x in rest:
        if x.isdisjoint(common):
          sets.append(x)
        else:
          merged = 1
          common |= x
      results.append(common)
    sets = results
  return sets


def readmerge(cin):
    # use this for combined culture plot
    prefix = os.path.splitext(os.path.basename(cin))[0]
    cout = os.path.abspath(os.path.join('..', 'output/data/' + prefix + "_processed.txt"))

    cs=[]
    with open(cin) as ccon:
        for line in ccon:
            if line[0] == '#':
              continue
            else:
              ids = line.strip()
              # strip [ and ] at start and end
              ids = ids[1:-1]
              ids = list(map(int, ids.split(',')))
              if len(ids) > 1:
                  cs.append(ids)
    
    cs = merge(cs)

    with open(cout, 'w') as f:
       for c in cs:
          f.write(str(list(c))+'\n')



def read(cflag):

    removalsfile = os.path.abspath(os.path.join('..', 'output/data/removals.txt'))
    with open(removalsfile) as rf:
        removals= dict(json.load(rf))
    removaldict = dict(removals)

    cin = os.path.abspath(os.path.join('..', 'data/' + cflag + "_const.txt"))
    cout = os.path.abspath(os.path.join('..', 'output/data/' + cflag + "_processed.txt"))

    cs=[]
    with open(cin) as ccon:
        for line in ccon:
            if line[0] == '#':
              continue
            else:
              ids = line.strip()
              # strip [ and ] at start and end
              ids = ids[1:-1]
              ids = list(map(int, ids.split(',')))
              ids = [removaldict[i] if i in removaldict else i for i in ids]

              if len(ids) > 1:
                  cs.append(ids)
    
    # no longer merge rows that have shared ids
    # cs = merge(cs)

    with open(cout, 'w') as f:
       for c in cs:
          f.write(str(list(c))+'\n')
            

