import csv
import json

hip2hr = []
hr2hip = []

with open('../data/hygfull.csv', 'r') as csvfile:
    hygreader = csv.reader(csvfile, delimiter=',')
    next(hygreader) # skip header
    for row in hygreader:
        hipnum = int(row[1])
        hnum = row[3]
        hnum = hnum.replace(' ', '')
        hip2hr.append( (hipnum, hnum) )
        hr2hip.append( (hnum, hipnum) )

with open("../output/data/hip2hr.txt", "w") as f:
    json.dump(hip2hr, f)

with open("../output/data/hr2hip.txt", "w") as f:
    json.dump(hr2hip, f)


