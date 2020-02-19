import sys
import os.path
import argparse
import re
import csv

dd = {}

with open('/scicomp/home/ydn3/SNP_Wrap_Python/testing/fauxEpidemiology.tsv', 'rU') as treeIDs:
    for line in treeIDs:
        if line.strip():
            lineList = re.split(r'\t+', line)
            ## return all but first element of array
            tempLine = lineList[1:]
            ## remove \n from all elements using map(lambda), then convert back to list
            dd[lineList[0]] = list(map(lambda s: s.strip(), tempLine))

## print key and comma-delimited list elements
#print('strain1690 -> ', end="")
## use splat operator *
#print(*dd['strain1690'], sep=',')

traction = {}

print('isolate,City,isTractionCity,body site,date of isolation')

## change dict, sorted by key, to make new field, isTractionCity
for key in sorted(dd.keys()):
    #print(dd[key][0])
    city = []
    if(re.search(pattern='Traction:', string=dd[key][0], flags=re.IGNORECASE)):
        city = dd[key][0].split(':')
        print(key + "," + city[1].strip() + ",Y," + dd[key][1] + "," + dd[key][2])
    elif(re.search(pattern='Static:', string=dd[key][0], flags=re.IGNORECASE)):
        city = dd[key][0].split(':')
        print(key + "," + city[1].strip() + ",N," + dd[key][1] + "," + dd[key][2])

