#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess
import pandas as pd
import numpy as np
import seaborn as sns
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


logger = logging.getLogger("intersectSNPs_BoxplotNoroSeriesIndiv.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions from two plain text, single-column files', usage="intersectSNPs_BoplotNoroSeriesIndiv.py miSeq1.snp.tsv iSeq1.snp.tsv miSeq2.snp.tsv iSeq2.snp.tsv miSeq3.snp.tsv iSeq3.snp.tsv miSeq4.snp.tsv iSeq4.snp.tsv miSeq5.snp.tsv iSeq5.snp.tsv miSeq6.snp.tsv iSeq6.snp.tsv miSeq7.snp.tsv iSeq7.snp.tsv miSeq8.snp.tsv iSeq8.snp.tsv")

parser.add_argument("listFile1", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile2", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile3", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile4", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile5", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile6", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile7", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile8", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile9", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile10", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile11", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile12", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile13", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile14", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile15", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile16", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--sampleType1', '-s1', default='', help="Experimental group for listFile1, if applicable")

parser.add_argument('--sampleType2', '-s2', default='', help="Experimental group for listFile2, if applicable")

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions?")

parser.add_argument('--outputType', '-o', default='S', choices=['P', 'C', 'I', 'S'], help="'P' = output a boxplot of SNPs positions, 'S' = output summary of intersection/union, 'C' = output single column, 'I' = output intersection in original tabular format.")

args = parser.parse_args()

columns1 = defaultdict(list)
columns2 = defaultdict(list)
columns3 = defaultdict(list)
columns4 = defaultdict(list)
columns5 = defaultdict(list)
columns6 = defaultdict(list)
columns7 = defaultdict(list)
columns8 = defaultdict(list)
columns9 = defaultdict(list)
columns10 = defaultdict(list)
columns11 = defaultdict(list)
columns12 = defaultdict(list)
columns13 = defaultdict(list)
columns14 = defaultdict(list)
columns15 = defaultdict(list)
columns16 = defaultdict(list)

header1 = []
header2 = []
header3 = []
header4 = []
header5 = []
header6 = []
header7 = []
header8 = []
header9 = []
header10 = []
header11 = []
header12 = []
header13 = []
header14 = []
header15 = []
header16 = []

fileName1 = os.path.basename(args.listFile1.name)
fileName2 = os.path.basename(args.listFile2.name)
fileName3 = os.path.basename(args.listFile3.name)
fileName4 = os.path.basename(args.listFile4.name)
fileName5 = os.path.basename(args.listFile5.name)
fileName6 = os.path.basename(args.listFile6.name)
fileName7 = os.path.basename(args.listFile7.name)
fileName8 = os.path.basename(args.listFile8.name)
fileName9 = os.path.basename(args.listFile9.name)
fileName10 = os.path.basename(args.listFile10.name)
fileName11 = os.path.basename(args.listFile11.name)
fileName12 = os.path.basename(args.listFile12.name)
fileName13 = os.path.basename(args.listFile13.name)
fileName14 = os.path.basename(args.listFile14.name)
fileName15 = os.path.basename(args.listFile15.name)
fileName16 = os.path.basename(args.listFile16.name)

## truncate sampleType1 and sampleType2 by splitting at whitespace 
## and showing only first word in input string

sample1 = args.sampleType1.split(r'\s+')
sample2 = args.sampleType2.split(r'\s+')

#fileName1 = fileName1 
#fileName2 = fileName2 
#fileName3 = fileName3
#fileName4 = fileName4

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

with open(args.listFile3.name, 'r') as table3:
    csvTable3 = csv.DictReader(table3, delimiter="\t")
    header3 = list(list(csvTable3)[0].keys())
    table3.seek(0)
    for row in csvTable3:
        for (key3, val3) in row.items():
            columns3[key3].append(val3)

with open(args.listFile4.name, 'r') as table4:
    csvTable4 = csv.DictReader(table4, delimiter="\t")
    header4 = list(list(csvTable4)[0].keys())
    table4.seek(0)
    for row in csvTable4:
        for (key4, val4) in row.items():
            columns4[key4].append(val4)

with open(args.listFile5.name, 'r') as table5:
    csvTable5 = csv.DictReader(table5, delimiter="\t")
    header5 = list(list(csvTable5)[0].keys())
    table5.seek(0)
    for row in csvTable5:
        for (key5, val5) in row.items():
            columns5[key5].append(val5)

with open(args.listFile6.name, 'r') as table6:
    csvTable6 = csv.DictReader(table6, delimiter="\t")
    header6 = list(list(csvTable6)[0].keys())
    table6.seek(0)
    for row in csvTable6:
        for (key6, val6) in row.items():
            columns6[key6].append(val6)

with open(args.listFile7.name, 'r') as table7:
    csvTable7 = csv.DictReader(table7, delimiter="\t")
    header7 = list(list(csvTable7)[0].keys())
    table7.seek(0)
    for row in csvTable7:
        for (key7, val7) in row.items():
            columns7[key7].append(val7)

with open(args.listFile8.name, 'r') as table8:
    csvTable8 = csv.DictReader(table8, delimiter="\t")
    header8 = list(list(csvTable8)[0].keys())
    table8.seek(0)
    for row in csvTable8:
        for (key8, val8) in row.items():
            columns8[key8].append(val8)

with open(args.listFile9.name, 'r') as table9:
    csvTable9 = csv.DictReader(table9, delimiter="\t")
    header9 = list(list(csvTable9)[0].keys())
    table9.seek(0)
    for row in csvTable9:
        for (key9, val9) in row.items():
            columns9[key9].append(val9)

with open(args.listFile10.name, 'r') as table10:
    csvTable10 = csv.DictReader(table10, delimiter="\t")
    header10 = list(list(csvTable10)[0].keys())
    table10.seek(0)
    for row in csvTable10:
        for (key10, val10) in row.items():
            columns10[key10].append(val10)

with open(args.listFile11.name, 'r') as table11:
    csvTable11 = csv.DictReader(table11, delimiter="\t")
    header11 = list(list(csvTable11)[0].keys())
    table11.seek(0)
    for row in csvTable11:
        for (key11, val11) in row.items():
            columns11[key11].append(val11)

with open(args.listFile12.name, 'r') as table12:
    csvTable12 = csv.DictReader(table12, delimiter="\t")
    header12 = list(list(csvTable12)[0].keys())
    table12.seek(0)
    for row in csvTable12:
        for (key12, val12) in row.items():
            columns12[key12].append(val12)

with open(args.listFile13.name, 'r') as table13:
    csvTable13 = csv.DictReader(table13, delimiter="\t")
    header13 = list(list(csvTable13)[0].keys())
    table13.seek(0)
    for row in csvTable13:
        for (key13, val13) in row.items():
            columns13[key13].append(val13)

with open(args.listFile14.name, 'r') as table14:
    csvTable14 = csv.DictReader(table14, delimiter="\t")
    header14 = list(list(csvTable14)[0].keys())
    table14.seek(0)
    for row in csvTable14:
        for (key14, val14) in row.items():
            columns14[key14].append(val14)

with open(args.listFile15.name, 'r') as table15:
    csvTable15 = csv.DictReader(table15, delimiter="\t")
    header15 = list(list(csvTable15)[0].keys())
    table15.seek(0)
    for row in csvTable15:
        for (key15, val15) in row.items():
            columns15[key15].append(val15)

with open(args.listFile16.name, 'r') as table16:
    csvTable16 = csv.DictReader(table16, delimiter="\t")
    header16 = list(list(csvTable16)[0].keys())
    table16.seek(0)
    for row in csvTable16:
        for (key16, val16) in row.items():
            columns16[key16].append(val16)


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
snpsPos3 = [x for x in header3 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos4 = [x for x in header4 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos5 = [x for x in header5 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos6 = [x for x in header6 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos7 = [x for x in header7 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos8 = [x for x in header8 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos9 = [x for x in header9 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos10 = [x for x in header10 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos11 = [x for x in header11 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos12 = [x for x in header12 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos13 = [x for x in header13 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos14 = [x for x in header14 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos15 = [x for x in header15 if re.search(r'(\bPOS\b|\bPosition\b)', x)]
snpsPos16 = [x for x in header16 if re.search(r'(\bPOS\b|\bPosition\b)', x)]

try:
    temp = snpsPos2[0]
except (IndexError):
    print("SNPs position header, 'POS', not found in {}".format(args.listFile2.name))
    sys.exit()


#print(columns1[snpsPos1[0]])
snpsFile1 = columns1[snpsPos1[0]]
snpsFile1.pop(0)
snpsFile1 = [int(y) for y in snpsFile1]

snpsFile2 = columns2[snpsPos2[0]]
snpsFile2.pop(0)
snpsFile2 = [int(z) for z in snpsFile2]

snpsFile3 = columns3[snpsPos3[0]]
snpsFile3.pop(0)
snpsFile3 = [int(z) for z in snpsFile3]

snpsFile4 = columns4[snpsPos4[0]]
snpsFile4.pop(0)
snpsFile4 = [int(w) for w in snpsFile4]

snpsFile5 = columns5[snpsPos5[0]]
snpsFile5.pop(0)
snpsFile5 = [int(z) for z in snpsFile5]

snpsFile6 = columns6[snpsPos6[0]]
snpsFile6.pop(0)
snpsFile6 = [int(z) for z in snpsFile6]

snpsFile7 = columns7[snpsPos7[0]]
snpsFile7.pop(0)
snpsFile7 = [int(z) for z in snpsFile7]

snpsFile8 = columns8[snpsPos8[0]]
snpsFile8.pop(0)
snpsFile8 = [int(z) for z in snpsFile8]

snpsFile9 = columns9[snpsPos9[0]]
snpsFile9.pop(0)
snpsFile9 = [int(y) for y in snpsFile9]

snpsFile10 = columns10[snpsPos10[0]]
snpsFile10.pop(0)
snpsFile10 = [int(z) for z in snpsFile10]

snpsFile11 = columns11[snpsPos11[0]]
snpsFile11.pop(0)
snpsFile11 = [int(z) for z in snpsFile11]

snpsFile12 = columns12[snpsPos12[0]]
snpsFile12.pop(0)
snpsFile12 = [int(w) for w in snpsFile12]

snpsFile13 = columns13[snpsPos13[0]]
snpsFile13.pop(0)
snpsFile13 = [int(z) for z in snpsFile13]

snpsFile14 = columns14[snpsPos14[0]]
snpsFile14.pop(0)
snpsFile14 = [int(z) for z in snpsFile14]

snpsFile15 = columns15[snpsPos15[0]]
snpsFile15.pop(0)
snpsFile15 = [int(z) for z in snpsFile15]

snpsFile16 = columns16[snpsPos16[0]]
snpsFile16.pop(0)
snpsFile16 = [int(z) for z in snpsFile16]


if(args.union == 'Y'):
    allPositions = Union(snpsPos1[0], snpsPos2[0])
    allPositions2 = Union(snpsPos3[0], snpsPos4[0])
    allPositions3 = Union(snpsPos5[0], snpsPos6[0])
    allPositions4 = Union(snpsPos7[0], snpsPos8[0])
else:
    allPositions = Intersection(snpsFile1, snpsFile2)
    allPositions2 = Intersection(snpsFile3, snpsFile4)
    allPositions3 = Intersection(snpsFile5, snpsFile6)
    allPositions4 = Intersection(snpsFile7, snpsFile8)
    allPositions5 = Intersection(snpsFile9, snpsFile10)
    allPositions6 = Intersection(snpsFile11, snpsFile12)
    allPositions7 = Intersection(snpsFile13, snpsFile14)
    allPositions8 = Intersection(snpsFile15, snpsFile16)

uniqueFile1 = Compl(snpsFile1, allPositions)
uniqueFile2 = Compl(snpsFile2, allPositions)


if(args.outputType == 'C'): ### output single column of positions, intersection or union
    for item in allPositions:
        print(item)
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
        columns3[refID[0]].pop(0)
        columns3[refBase[0]].pop(0)
        columns3[altBase[0]].pop(0)
        columns3[stats[0]].pop(0)
        columns3[depth[0]].pop(0)
        columns3[consensus[0]].pop(0)
        columns4[refID[0]].pop(0)
        columns4[refBase[0]].pop(0)
        columns4[altBase[0]].pop(0)
        columns4[stats[0]].pop(0)
        columns4[depth[0]].pop(0)
        columns4[consensus[0]].pop(0)
        columns5[refID[0]].pop(0)
        columns5[refBase[0]].pop(0)
        columns5[altBase[0]].pop(0)
        columns5[stats[0]].pop(0)
        columns5[depth[0]].pop(0)
        columns5[consensus[0]].pop(0)
        columns6[refID[0]].pop(0)
        columns6[refBase[0]].pop(0)
        columns6[altBase[0]].pop(0)
        columns6[stats[0]].pop(0)
        columns6[depth[0]].pop(0)
        columns6[consensus[0]].pop(0)
        columns7[refID[0]].pop(0)
        columns7[refBase[0]].pop(0)
        columns7[altBase[0]].pop(0)
        columns7[stats[0]].pop(0)
        columns7[depth[0]].pop(0)
        columns7[consensus[0]].pop(0)
        columns8[refID[0]].pop(0)
        columns8[refBase[0]].pop(0)
        columns8[altBase[0]].pop(0)
        columns8[stats[0]].pop(0)
        columns8[depth[0]].pop(0)
        columns8[consensus[0]].pop(0)
        columns9[refID[0]].pop(0)
        columns9[refBase[0]].pop(0)
        columns9[altBase[0]].pop(0)
        columns9[stats[0]].pop(0)
        columns9[depth[0]].pop(0)
        columns9[consensus[0]].pop(0)
        columns10[refID[0]].pop(0)
        columns10[refBase[0]].pop(0)
        columns10[altBase[0]].pop(0)
        columns10[stats[0]].pop(0)
        columns10[depth[0]].pop(0)
        columns10[consensus[0]].pop(0)
        columns11[refID[0]].pop(0)
        columns11[refBase[0]].pop(0)
        columns11[altBase[0]].pop(0)
        columns11[stats[0]].pop(0)
        columns11[depth[0]].pop(0)
        columns11[consensus[0]].pop(0)
        columns12[refID[0]].pop(0)
        columns12[refBase[0]].pop(0)
        columns12[altBase[0]].pop(0)
        columns12[stats[0]].pop(0)
        columns12[depth[0]].pop(0)
        columns12[consensus[0]].pop(0)
        columns13[refID[0]].pop(0)
        columns13[refBase[0]].pop(0)
        columns13[altBase[0]].pop(0)
        columns13[stats[0]].pop(0)
        columns13[depth[0]].pop(0)
        columns13[consensus[0]].pop(0)
        columns14[refID[0]].pop(0)
        columns14[refBase[0]].pop(0)
        columns14[altBase[0]].pop(0)
        columns14[stats[0]].pop(0)
        columns14[depth[0]].pop(0)
        columns14[consensus[0]].pop(0)
        columns15[refID[0]].pop(0)
        columns15[refBase[0]].pop(0)
        columns15[altBase[0]].pop(0)
        columns15[stats[0]].pop(0)
        columns15[depth[0]].pop(0)
        columns15[consensus[0]].pop(0)
        columns16[refID[0]].pop(0)
        columns16[refBase[0]].pop(0)
        columns16[altBase[0]].pop(0)
        columns16[stats[0]].pop(0)
        columns16[depth[0]].pop(0)
        columns16[consensus[0]].pop(0)

    except (IndexError):
        print("Improperly-formatted headers in {}".format(args.listFile1.name))
        sys.exit()

    intStats1 = []
    intStats2 = []
    intStats3 = []
    intStats4 = []
    intStats5 = []
    intStats6 = []
    intStats7 = []
    intStats8 = []
    intStats9 = []
    intStats10 = []
    intStats11 = []
    intStats12 = []
    intStats13 = []
    intStats14 = []
    intStats15 = []
    intStats16 = []


    coverage1 = []
    coverage2 = []
    coverage3 = []
    coverage4 = []
    coverage5 = []
    coverage6 = []
    coverage7 = []
    coverage8 = []
    coverage9 = []
    coverage10 = []
    coverage11 = []
    coverage12 = []
    coverage13 = []
    coverage14 = []
    coverage15 = []
    coverage16 = []


    for item in allPositions:
        idx = 0
        while(idx < len(columns1[snpsPos1[0]])):
            if( (item == int(columns1[snpsPos1[0]][idx]) and (float(columns1[consensus[0]][idx]) > 0.15) )):
               intStats1.append( int(columns1[stats[0]][idx]) )
               coverage1.append( int(columns1[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions:
        idx = 0
        while(idx < len(columns2[snpsPos1[0]])):
            if( (item == int(columns2[snpsPos1[0]][idx]) and (float(columns2[consensus[0]][idx]) > 0.15) )):
                intStats2.append( int(columns2[stats[0]][idx]) )
                coverage2.append( int(columns2[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions2:
        idx = 0
        while(idx < len(columns3[snpsPos1[0]])):
            if( (item == int(columns3[snpsPos1[0]][idx]) and (float(columns3[consensus[0]][idx]) > 0.15) )):
                intStats3.append( int(columns3[stats[0]][idx]) )
                coverage3.append( int(columns3[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions2:
        idx = 0
        while(idx < len(columns4[snpsPos1[0]])):
            if( (item == int(columns4[snpsPos1[0]][idx]) and (float(columns4[consensus[0]][idx]) > 0.15) )):
                intStats4.append( int(columns4[stats[0]][idx]) )
                coverage4.append( int(columns4[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions3:
        idx = 0
        while(idx < len(columns5[snpsPos1[0]])):
            if( (item == int(columns5[snpsPos1[0]][idx]) and (float(columns5[consensus[0]][idx]) > 0.15) )):
                intStats5.append( int(columns5[stats[0]][idx]) )
                coverage5.append( int(columns5[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions3:
        idx = 0
        while(idx < len(columns6[snpsPos1[0]])):
            if( (item == int(columns6[snpsPos1[0]][idx]) and (float(columns6[consensus[0]][idx]) > 0.15) )):
                intStats6.append( int(columns6[stats[0]][idx]) )
                coverage6.append( int(columns6[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions4:
        idx = 0
        while(idx < len(columns7[snpsPos1[0]])):
            if( (item == int(columns7[snpsPos1[0]][idx]) and (float(columns7[consensus[0]][idx]) > 0.15) )):
                intStats7.append( int(columns7[stats[0]][idx]) )
                coverage7.append( int(columns7[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions4:
        idx = 0
        while(idx < len(columns8[snpsPos1[0]])):
            if( (item == int(columns8[snpsPos1[0]][idx]) and (float(columns8[consensus[0]][idx]) > 0.15) )):
                intStats8.append( int(columns8[stats[0]][idx]) )
                coverage8.append( int(columns8[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions5:
        idx = 0
        while(idx < len(columns9[snpsPos9[0]])):
            if( (item == int(columns9[snpsPos5[0]][idx]) and (float(columns9[consensus[0]][idx]) > 0.15) )):
               intStats9.append( int(columns9[stats[0]][idx]) )
               coverage9.append( int(columns9[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions5:
        idx = 0
        while(idx < len(columns10[snpsPos1[0]])):
            if( (item == int(columns10[snpsPos1[0]][idx]) and (float(columns10[consensus[0]][idx]) > 0.15) )):
                intStats10.append( int(columns10[stats[0]][idx]) )
                coverage10.append( int(columns10[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions6:
        idx = 0
        while(idx < len(columns11[snpsPos1[0]])):
            if( (item == int(columns11[snpsPos1[0]][idx]) and (float(columns11[consensus[0]][idx]) > 0.15) )):
                intStats11.append( int(columns11[stats[0]][idx]) )
                coverage11.append( int(columns11[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions6:
        idx = 0
        while(idx < len(columns12[snpsPos1[0]])):
            if( (item == int(columns12[snpsPos1[0]][idx]) and (float(columns12[consensus[0]][idx]) > 0.15) )):
                intStats12.append( int(columns12[stats[0]][idx]) )
                coverage12.append( int(columns12[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions7:
        idx = 0
        while(idx < len(columns13[snpsPos1[0]])):
            if( (item == int(columns13[snpsPos1[0]][idx]) and (float(columns13[consensus[0]][idx]) > 0.15) )):
                intStats13.append( int(columns13[stats[0]][idx]) )
                coverage13.append( int(columns13[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions7:
        idx = 0
        while(idx < len(columns14[snpsPos1[0]])):
            if( (item == int(columns14[snpsPos1[0]][idx]) and (float(columns14[consensus[0]][idx]) > 0.15) )):
                intStats14.append( int(columns14[stats[0]][idx]) )
                coverage14.append( int(columns14[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions8:
        idx = 0
        while(idx < len(columns15[snpsPos1[0]])):
            if( (item == int(columns15[snpsPos1[0]][idx]) and (float(columns15[consensus[0]][idx]) > 0.15) )):
                intStats15.append( int(columns15[stats[0]][idx]) )
                coverage15.append( int(columns15[depth[0]][idx]) )
            idx = idx + 1

    for item in allPositions8:
        idx = 0
        while(idx < len(columns16[snpsPos1[0]])):
            if( (item == int(columns16[snpsPos1[0]][idx]) and (float(columns16[consensus[0]][idx]) > 0.15) )):
                intStats16.append( int(columns16[stats[0]][idx]) )
                coverage16.append( int(columns16[depth[0]][idx]) )
            idx = idx + 1



    file1 = re.sub(r'\.snp\.tsv', r'', fileName1)
    file2 = re.sub(r'\.snp\.tsv', r'', fileName2)
    file3 = re.sub(r'\.snp\.tsv', r'', fileName3)
    file4 = re.sub(r'\.snp\.tsv', r'', fileName4)
    file5 = re.sub(r'\.snp\.tsv', r'', fileName5)
    file6 = re.sub(r'\.snp\.tsv', r'', fileName6)
    file7 = re.sub(r'\.snp\.tsv', r'', fileName7)
    file8 = re.sub(r'\.snp\.tsv', r'', fileName8)
    file9 = re.sub(r'\.snp\.tsv', r'', fileName9)
    file10 = re.sub(r'\.snp\.tsv', r'', fileName10)
    file11 = re.sub(r'\.snp\.tsv', r'', fileName11)
    file12 = re.sub(r'\.snp\.tsv', r'', fileName12)
    file13 = re.sub(r'\.snp\.tsv', r'', fileName13)
    file14 = re.sub(r'\.snp\.tsv', r'', fileName14)
    file15 = re.sub(r'\.snp\.tsv', r'', fileName15)
    file16 = re.sub(r'\.snp\.tsv', r'', fileName16)


    set1 = pd.DataFrame({'SNP_Depth': coverage1, 'Sample' : np.repeat(['sample-1'], len(coverage1)), 'Instrument' : np.repeat(['MiSeq'], len(coverage1))})
    set2 = pd.DataFrame({'SNP_Depth': coverage2, 'Sample' : np.repeat(['sample-1'], len(coverage2)), 'Instrument' : np.repeat(['iSeq'], len(coverage2))})
    set3 = pd.DataFrame({'SNP_Depth': coverage3, 'Sample' : np.repeat(['sample-2'], len(coverage3)), 'Instrument' : np.repeat(['MiSeq'], len(coverage3))})
    set4 = pd.DataFrame({'SNP_Depth': coverage4, 'Sample' : np.repeat(['sample-2'], len(coverage4)), 'Instrument' : np.repeat(['iSeq'], len(coverage4))})
    set5 = pd.DataFrame({'SNP_Depth': coverage5, 'Sample' : np.repeat(['sample-3'], len(coverage5)), 'Instrument' : np.repeat(['MiSeq'], len(coverage5))})
    set6 = pd.DataFrame({'SNP_Depth': coverage6, 'Sample' : np.repeat(['sample-3'], len(coverage6)), 'Instrument' : np.repeat(['iSeq'], len(coverage6))})
    set7 = pd.DataFrame({'SNP_Depth': coverage7, 'Sample' : np.repeat(['sample-4'], len(coverage7)), 'Instrument' : np.repeat(['MiSeq'], len(coverage7))})
    set8 = pd.DataFrame({'SNP_Depth': coverage8, 'Sample' : np.repeat(['sample-4'], len(coverage8)), 'Instrument' : np.repeat(['iSeq'], len(coverage8))})
    set9 = pd.DataFrame({'SNP_Depth': coverage9, 'Sample' : np.repeat(['sample-5'], len(coverage9)), 'Instrument' : np.repeat(['MiSeq'], len(coverage9))})
    set10 = pd.DataFrame({'SNP_Depth': coverage10, 'Sample' : np.repeat(['sample-5'], len(coverage10)), 'Instrument' : np.repeat(['iSeq'], len(coverage10))})
    set11 = pd.DataFrame({'SNP_Depth': coverage11, 'Sample' : np.repeat(['sample-6'], len(coverage11)), 'Instrument' : np.repeat(['MiSeq'], len(coverage11))})
    set12 = pd.DataFrame({'SNP_Depth': coverage12, 'Sample' : np.repeat(['sample-6'], len(coverage12)), 'Instrument' : np.repeat(['iSeq'], len(coverage12))})
    set13 = pd.DataFrame({'SNP_Depth': coverage13, 'Sample' : np.repeat(['sample-7'], len(coverage13)), 'Instrument' : np.repeat(['MiSeq'], len(coverage13))})
    set14 = pd.DataFrame({'SNP_Depth': coverage14, 'Sample' : np.repeat(['sample-7'], len(coverage14)), 'Instrument' : np.repeat(['iSeq'], len(coverage14))})
    set15 = pd.DataFrame({'SNP_Depth': coverage15, 'Sample' : np.repeat(['sample-8'], len(coverage15)), 'Instrument' : np.repeat(['MiSeq'], len(coverage15))})
    set16 = pd.DataFrame({'SNP_Depth': coverage16, 'Sample' : np.repeat(['sample-8'], len(coverage16)), 'Instrument' : np.repeat(['iSeq'], len(coverage16))})


#    data1 = pd.DataFrame({'MiSeq' : intStats1, 'iSeq' : intStats2, 'MiSeq' : intStats3, 'iSeq' : intStats4})
    data2 = pd.concat([set1, set2, set3, set4, set5, set6, set7, set8, set9, set10, set11, set12, set13, set14, set15, set16])
#    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, figsize=(8.5,11))
    fig = plt.figure(figsize =(8.5, 11))
#    plt.boxplot( data1, whis=50, labels=ticks )
#    plt.title("Quality of Shared SNPs")
#    plt.ylabel('SNP Quality')
    sns.set(font_scale=1.2)
    sns.boxplot( x=data2['Instrument'], y=data2['SNP_Depth'], dodge=True )


#    sns.boxplot( x=data2['Sample'], y=data2['SNP_Depth'], hue=data2['Instrument'], dodge=True )
    #plt.title("Coverage of Shared SNPs")
    #plt.ylabel('SNP Depth')
    plt.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_pandas/SNP_SARS-CoV-2_seabornTrimmoPooled.png')
