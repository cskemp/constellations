import csv
import json

hrnames = []

with open('../data/hygfull.csv', 'r') as csvfile:
    hygreader = csv.reader(csvfile, delimiter=',')
    next(hygreader) # skip header
    for row in hygreader:
        hnum = row[3]
        hnum = hnum.replace(' ', '')
        name = row[5]
        if not(hnum ==  ''):
            name = name.replace(' ', '')
            if not(name):
                name = hnum
            # print(hnum  + ' ' + name)
            hrnames.append( (hnum, name) )              

# manual additions
hrnames.append( ('5506', 'EpsBoo') )   # 5505 and 5506 form a binary star
hrnames.append( ('5834', 'Zet2CrB') )

with open("../output/data/hrnames.txt", "w") as f:
    json.dump(hrnames, f)


