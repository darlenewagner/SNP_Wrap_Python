#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess
from collections import defaultdict

## Reads an arbitrary number of snp.tsv files to compute union of SNPs positions
## The seven columns should have the following headers:
## CHROM   POS     REF     ALT     QUAL    INFO-DP CONSENS

## Function: Union 

def Union(list1, list2):
    return list(set().union(list1, list2))

## Simple filename extractor function

def getIsolateStr(filePathString):
    splitStr = re.split(pattern='/', string=filePathString)
    fileNameIdx = len(splitStr) - 1
    isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
    return isolateString[0]

SampleSheet = dict()
samples = []

### Create a triple-nested dictionary

for argument in sys.argv[1:]:
    columns = defaultdict(list)
    header = []
    samples.append(argument)
    with open(argument, 'r') as table:
        csvTable = csv.DictReader(table, delimiter="\t")
        table.seek(0)
        isolateStr = getIsolateStr(argument)
        index = 0
        Lines = dict()
        #columns["Position"]=defaultdict(list)
        for row in csvTable:
            tempStr = {'Ref_ID' : row['Ref_ID'], 'Position' : row['Position'], 'Ref' : row['Ref'], 'Alt' : row['Alt'], 'Qual' : row['Qual'], 'Coverage' : row['Coverage'], 'Fraction' : row['Fraction']}
            # print(tempStr)
            #isolateStr = isolateStr + "_" + str(index)
            ## print(isolateStr)
            Lines[str(index)] = tempStr
            index = index + 1

        SampleSheet[argument] = Lines


myList = []
unionList = []

#for keys, values in Lines.items():
#    values['Position']

i = 0
while i < len(samples):
    for keys1,values1 in SampleSheet[samples[i]].items():
        myList.append(values1['Position'])
    unionList = Union(unionList, myList)
    #print(len(unionList))
    myList = []
    i = i + 1

print("Union of SNPs positions size is:  ", len(unionList))
        
    #for keys, values in SampleSheet[x].items():
    #    print(values)

