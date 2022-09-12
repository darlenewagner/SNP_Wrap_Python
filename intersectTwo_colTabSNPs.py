#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess
from collections import defaultdict

## Reads two five-column or six-column .txt, .csv, or .tsv files to compute intersection of SNPs positions

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

logger = logging.getLogger("intersectTwo_colTabSNPs.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions from two plain text, single-column files', usage="intersectTwo_colTabSNPs.py snpFile1.tsv snpFile2.tsv")

parser.add_argument("listFile1", type=tsv_check('.tab', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile2", type=tsv_check('.tab', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions?")

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'C', 'T'], help="'S' = output summary of intersection/union, 'C' = output single column, or 'T' = output in original tabular format.")

args = parser.parse_args()

columns1 = defaultdict(list)
columns2 = defaultdict(list)
header1 = []
header2 = []

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
depth = [x for x in header1 if re.search(r'\bQUAL\b', x)]
quality = [x for x in header1 if re.search(r'\bINFO-DP\b', x)]

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

uniqueFile1 = len(snpsFile1) - len(allPositions)
uniqueFile2 = len(snpsFile2) - len(allPositions)


if(args.outputType == 'C'): ### output single column of positions, intersection or union
    for item in allPositions:
        print(item)
elif(args.outputType == 'T'): ### output six-column .tsv format
    try:
        columns1[refID[0]].pop(0)
        columns1[refBase[0]].pop(0)
        columns1[altBase[0]].pop(0)
        columns1[depth[0]].pop(0)
        columns1[quality[0]].pop(0)
    except (IndexError):
        print("Improperly-formatted headers in {}".format(args.listFile1.name))
        sys.exit()
    print(refID[0] + "\t" + snpsPos1[0] + "\t" + refBase[0] + "\t" + altBase[0] + "\t" + quality[0] + "\t" + depth[0])
    for item in allPositions:
        idx = 0
        while(idx < len(columns1[snpsPos1[0]])):
            if( item == int(columns1[snpsPos1[0]][idx] ) ):
                print(columns1[refID[0]][idx] + "\t" + columns1[snpsPos1[0]][idx] + "\t" + columns1[refBase[0]][idx] + "\t" + columns1[altBase[0]][idx] + "\t" + columns1[depth[0]][idx] + "\t" + columns1[quality[0]][idx])
            idx = idx + 1
elif((args.outputType == 'S') and (args.union == 'N')):  ### output verbose description of intersection
    print("{} has {} unique snps and {} has {} unique snps.".format(args.listFile1.name, str(uniqueFile1), args.listFile2.name, str(uniqueFile2)) )
    if(args.union == 'Y'):
        print("Their union contains {} snps.".format(str(len(allPositions))))
    else:
        print("Their intersection contains {} snps.".format(str(len(allPositions))))
