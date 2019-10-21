#!/usr/bin/python

import sys
import os.path
import argparse
import re
import logging
import warnings
import csv
import subprocess

## Reads two LyveSET out.filteredMatrix.tsv files to output SNP positions held in common (intersection)
## Also can output union of SNPs with argument --union Y

## Function: A closure for .tsv or .csv extension checking

def tsv_check(expected_ext1, expected_ext2, openner):
    def extension(filename):
        if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
            raise ValueError()
        return openner(filename)
    return extension

## Function: Intersection (default)

def Intersection(list1, list2):
    return sorted(set(list1).intersection(list2))

## Note that python sorted() coverts set back to list


## Function: Union or Full Outer Join

def Union(list1, list2):
    final_list = list(set().union(list1, list2))
    return final_list

logger = logging.getLogger("filteredMatrixIntersect.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions', usage="filteredMatrixIntersect.py LyveSETproj1/msa/out.filteredMatrix.tsv LyveSETproj2/msa/out.filteredMatrix.tsv --noNs [Y/N] --union [Y/N]")

parser.add_argument("matrix1", type=tsv_check('.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("matrix2", type=tsv_check('.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions")

parser.add_argument('--noNs', '-n', default='Y', choices=['Y','N'], help="Calculate union for fully informative positions only")


args = parser.parse_args()

table1 = open(args.matrix1.name, 'r')

table2 = open(args.matrix2.name, 'r')

snpMatrix1 = csv.reader(table1, delimiter='\t')

snpMatrix2 = csv.reader(table2, delimiter='\t')

textMatrix1 = [row for row in snpMatrix1]

textMatrix2 = [row for row in snpMatrix2]

## First, parse textMatrix1, retrieve Chromosome and position, then concatenate them

idx = 0

snpPosition1 = []

while idx < len(textMatrix1):
    hasAmbig = 0
    if(not(re.search(pattern='^#', string=textMatrix1[idx][0]))):
        for item in textMatrix1[idx]:
            if((re.search(pattern='N', string=item)) and (args.noNs == 'Y')):
               hasAmbig = 1
        if(hasAmbig == 0):
            snpPosition1.append(textMatrix1[idx][0] + "_" + textMatrix1[idx][1])
    idx = idx + 1

## Next, parse textMatrix2, concatenate Chromosome and position, compare each to snpPosition elements

idx = 0

snpPosition2 = []

while idx < len(textMatrix2):
    hasAmbig = 0
    if(not(re.search(pattern='^#', string=textMatrix2[idx][0]))):
        for item in textMatrix2[idx]:
            if((re.search(pattern='N', string=item)) and (args.noNs == 'Y')):
                hasAmbig = 1
        if(hasAmbig == 0):
            snpPosition2.append(textMatrix2[idx][0] + "_" + textMatrix2[idx][1]) 
    idx = idx + 1



if(args.union == 'Y'):
    allPositions = Union(snpPosition1, snpPosition2)
    for item in allPositions:
        print(item)
else:
    allPositions = Intersection(snpPosition1, snpPosition2)
    for item in allPositions:
        print(item)





        


