#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

## Reads two seven-column .txt, .csv, or .tsv files to compute intersection of SNPs positions
## The seven columns should have the following headers:
## CHROM   POS     REF     ALT     QUAL    INFO-DP CONSENS

## Function: A closure for .tsv or .csv extension checking

def tsv_check(expected_ext1, expected_ext2, expected_ext3, openner):
    def extension(filename):
        if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
            raise ValueError()
        return openner(filename)
    return extension

## Function: Intersection (default)

def Intersection(list1, list2):
    return sorted(set(list1).intersection(list2))

## Function: Union 

def Union(list1, list2):
    return list(set().union(list1, list2))

## Function: Complement

def Compl(list1, intersec):
    snpSet1 = set(list1)
    snpSet2 = set(intersec)
    return list(snpSet1 - snpSet2)

logger = logging.getLogger("intersectAndComplement_colTabSNPs.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions from two plain text, single-column files', usage="intersectAndComplement_colTabSNPs.py snpFile1.tsv snpFile2.tsv")

parser.add_argument("listFile1", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile2", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--sampleType1', '-s1', default='', help="Experimental group for listFile1, if applicable")

parser.add_argument('--sampleType2', '-s2', default='', help="Experimental group for listFile2, if applicable")

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions?")

parser.add_argument('--outputType', '-o', default='S', choices=['P', 'C', 'D', 'I', 'S'], help="'P' = output a lollipop plot of SNPs positions, 'S' = output summary of intersection/union, 'C' = output single column, 'I' = output intersection in original tabular format, or 'D' output non-intersection snps in tabular format.")

args = parser.parse_args()

columns1 = defaultdict(list)
columns2 = defaultdict(list)
header1 = []
header2 = []

fileName1 = os.path.basename(args.listFile1.name)
fileName2 = os.path.basename(args.listFile2.name)

## truncate sampleType1 and sampleType2 by splitting at whitespace 
## and showing only first word in input string

sample1 = args.sampleType1.split(r'\s+')
sample2 = args.sampleType2.split(r'\s+')

fileName1 = fileName1 + " " + sample1[0]
fileName2 = fileName2 + " " + sample2[0]

with open(args.listFile1.name, 'r') as table1:
    csvTable1 = csv.DictReader(table1, delimiter="\t")
    header1 = list(list(csvTable1)[0].keys())
    table1.seek(0)
    for row in csvTable1:
        for (key, val) in row.items():
            columns1[key].append(val)

with open(args.listFile2.name, 'r') as table2:
    csvTable2 = csv.DictReader(table2, delimiter="\t")
    header2 = list(list(csvTable2)[0].keys())
    table2.seek(0)
    for row in csvTable2:
        for (key2, val2) in row.items():
            columns2[key2].append(val2)

#print(header1)

snpsPos1 = [x for x in header1 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
refID = [x for x in header1 if re.search(r'(\bCHROM\b|\bRef(\s|_)ID\b)', x)]
refBase = [x for x in header1 if re.search(r'(\bRef\b|\bREF\b)', x)]
altBase = [x for x in header1 if re.search(r'(\bAlt\b|\bALT\b)', x)]
stats = [x for x in header1 if re.search(r'(\bQual\b|\bQUAL\b)', x)]
depth = [x for x in header1 if re.search(r'(\bINFO\-DP\b|\bCoverage\b)', x)]
consensus = [x for x in header1 if re.search(r'(\bCONSENS\b|\bFraction\b)', x)]

try:
    temp = snpsPos1[0]
except (IndexError):
    print("SNPs position header, 'POS', not found in {}".format(args.listFile1.name))
    sys.exit()

snpsPos2 = [x for x in header2 if re.search(r'(\bPOS\b|\bPosition\b)', x)]


try:
    temp = snpsPos2[0]
except (IndexError):
    print("SNPs position header, 'POS', not found in {}".format(args.listFile2.name))
    sys.exit()


#print(columns1[snpsPos1[0]])
snpsFile1 = columns1[snpsPos1[0]]
snpsFile1.pop(0)
snpsFile1 = [int(y) for y in snpsFile1]

#print(columns2[snpsPos2[0]])
snpsFile2 = columns2[snpsPos2[0]]
snpsFile2.pop(0)
snpsFile2 = [int(z) for z in snpsFile2]


if(args.union == 'Y'):
    allPositions = Union(snpsPos1[0], snpsPos2[0])
else:
    allPositions = Intersection(snpsFile1, snpsFile2)

uniqueFile1 = Compl(snpsFile1, allPositions)
uniqueFile2 = Compl(snpsFile2, allPositions)


if(args.outputType == 'C'): ### output single column of positions, intersection or union
    for item in allPositions:
        print(item)
elif(args.outputType == 'D'): ### output five-column .tsv of non-intersection
    try:
        columns1[refID[0]].pop(0)
        columns1[refBase[0]].pop(0)
        columns1[altBase[0]].pop(0)
        columns1[stats[0]].pop(0)
        columns1[depth[0]].pop(0)
        columns1[consensus[0]].pop(0)
        columns2[refID[0]].pop(0)
        columns2[refBase[0]].pop(0)
        columns2[altBase[0]].pop(0)
        columns2[stats[0]].pop(0)
        columns2[depth[0]].pop(0)
        columns2[consensus[0]].pop(0)
    except (IndexError):
        print("Improperly-formatted headers in {}".format(args.listFile1.name))
        sys.exit()
    print("Unique to {}".format(fileName1))
    print(refID[0] + "\t" + snpsPos1[0] + "\t" + refBase[0] + "\t" + altBase[0] + "\t" + stats[0] + "\t" + depth[0] + "\t" + consensus[0])
    for item in uniqueFile1:
        idx = 0
        while(idx < len(columns1[snpsPos1[0]])):
            if( item == int(columns1[snpsPos1[0]][idx]) and ( float(columns1[consensus[0]][idx]) > 0.25)):
                print(columns1[refID[0]][idx] + "\t" + columns1[snpsPos1[0]][idx] + "\t" + columns1[refBase[0]][idx] + "\t" + columns1[altBase[0]][idx] + "\t" + columns1[stats[0]][idx] + "\t" + columns1[depth[0]][idx]+ "\t" + columns1[consensus[0]][idx])
            idx = idx + 1
    print("Unique to {}".format(fileName2))
    print(refID[0] + "\t" + snpsPos1[0] + "\t" + refBase[0] + "\t" + altBase[0] + "\t" + stats[0] + "\t" + depth[0] + "\t" + consensus[0])
    for item in uniqueFile2:
        idx = 0
        while(idx < len(columns2[snpsPos2[0]]) ):
            if( item == int(columns2[snpsPos2[0]][idx] ) and ( float(columns2[consensus[0]][idx]) > 0.25)):
                print(columns2[refID[0]][idx] + "\t" + columns2[snpsPos2[0]][idx] + "\t" + columns2[refBase[0]][idx] + "\t" + columns2[altBase[0]][idx] + "\t" + columns2[stats[0]][idx] + "\t" + columns2[depth[0]][idx] + "\t" + columns1[consensus[0]][idx])
            idx = idx + 1

elif(args.outputType == 'I'): ### output five-column .tsv format
    try:
        columns1[refID[0]].pop(0)
        columns1[refBase[0]].pop(0)
        columns1[altBase[0]].pop(0)
        columns1[stats[0]].pop(0)
        columns1[depth[0]].pop(0)
        columns1[consensus[0]].pop(0)
        columns2[refID[0]].pop(0)
        columns2[refBase[0]].pop(0)
        columns2[altBase[0]].pop(0)
        columns2[stats[0]].pop(0)
        columns2[depth[0]].pop(0)
        columns2[consensus[0]].pop(0)
    except (IndexError):
        print("Improperly-formatted headers in {}".format(args.listFile1.name))
        sys.exit()
    print("Shared SNPs in {}".format(fileName2))
    print(refID[0] + "\t" + snpsPos1[0] + "\t" + refBase[0] + "\t" + altBase[0] + "\t" + stats[0] + "\t" + depth[0] + "\t" + consensus[0])
    for item in allPositions:
        idx = 0
        while(idx < len(columns1[snpsPos1[0]])):
            if( (item == int(columns1[snpsPos1[0]][idx]) and (float(columns1[consensus[0]][idx]) > 0.25) )):
                print(columns1[refID[0]][idx] + "\t" + columns1[snpsPos1[0]][idx] + "\t" + columns1[refBase[0]][idx] + "\t" + columns1[altBase[0]][idx] + "\t" + columns1[stats[0]][idx] + "\t" + columns1[depth[0]][idx] + "\t" + columns1[consensus[0]][idx])
            idx = idx + 1
        
    print("Shared SNPs in {}".format(fileName2))
    print(refID[0] + "\t" + snpsPos1[0] + "\t" + refBase[0] + "\t" + altBase[0] + "\t" + stats[0] + "\t" + depth[0] + "\t" + consensus[0])
    for item in allPositions:
        idx = 0
        while(idx < len(columns2[snpsPos1[0]])):
            if( (item == int(columns2[snpsPos1[0]][idx]) and (float(columns2[consensus[0]][idx]) > 0.25) )):
                print(columns2[refID[0]][idx] + "\t" + columns2[snpsPos1[0]][idx] + "\t" + columns2[refBase[0]][idx] + "\t" + columns2[altBase[0]][idx] + "\t" + columns2[stats[0]][idx] + "\t" + columns2[depth[0]][idx] + "\t" + columns2[consensus[0]][idx])
            idx = idx + 1

elif((args.outputType == 'S') and (args.union == 'N')):  ### output verbose description of intersection
    print("{} has {} unique snps and {} has {} unique snps.".format(fileName1, str(len(list(uniqueFile1))), fileName2, str(len(list(uniqueFile2)))) )
    if(args.union == 'Y'):
        print("Their union contains {} snps.".format(str(len(allPositions))))
    else:
        print("Their intersection contains {} snps.".format(str(len(allPositions))))
elif((args.outputType == 'P') and (args.union == 'N')):  ### output scatterplots
    try:
        columns1[refID[0]].pop(0)
        columns1[refBase[0]].pop(0)
        columns1[altBase[0]].pop(0)
        columns1[stats[0]].pop(0)
        columns1[depth[0]].pop(0)
        columns1[consensus[0]].pop(0)
        columns2[refID[0]].pop(0)
        columns2[refBase[0]].pop(0)
        columns2[altBase[0]].pop(0)
        columns2[stats[0]].pop(0)
        columns2[depth[0]].pop(0)
        columns2[consensus[0]].pop(0)
    except (IndexError):
        print("Improperly-formatted headers in {}".format(args.listFile1.name))
        sys.exit()

    intStats1 = []
    intStats2 = []
    coverage1 = []
    coverage2 = []

    for item in allPositions:
        idx = 0
        while(idx < len(columns1[snpsPos1[0]])):
            if( (item == int(columns1[snpsPos1[0]][idx]) and (float(columns1[consensus[0]][idx]) > 0.25) )):
               intStats1.append( int(columns1[stats[0]][idx]) )
               coverage1.append( int(columns1[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions:
        idx = 0
        while(idx < len(columns2[snpsPos1[0]])):
            if( (item == int(columns2[snpsPos1[0]][idx]) and (float(columns2[consensus[0]][idx]) > 0.25) )):
                intStats2.append( int(columns2[stats[0]][idx]) )
                coverage2.append( int(columns2[depth[0]][idx]) )
            idx = idx + 1
    
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(8.5,11))
    axes[0].scatter(allPositions, intStats1, s=15, c='b', marker="s", label=fileName1)
    axes[0].scatter(allPositions, intStats2, s=15, c='r', marker="o", label=fileName2)
    axes[0].set_title("Quality of Shared SNPs")
    axes[1].set_xlabel('Reference Position')
    axes[0].set_ylabel('SNP Quality')
    axes[1].scatter(allPositions, coverage1, s=15, c='b', marker="s")
    axes[1].scatter(allPositions, coverage2, s=15, c='r', marker="o")
    axes[1].set_title("Depth of Shared SNPs")
    axes[1].set_xlabel('Reference Position')
    axes[1].set_ylabel('SNP Depth')
    fig.legend(loc='lower right')
    fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_pandas/SNP_SARS-CoV-2_scatter.png')
