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


logger = logging.getLogger("columnOfSingleMatrix.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions', usage="columnOfSingleMatrix.py LyveSETproj1/msa/out.filteredMatrix.tsv --noNs [Y/N]")

parser.add_argument("matrix1", type=tsv_check('.tsv', '.csv', argparse.FileType('r')))

#parser.add_argument("matrix2", type=tsv_check('.tsv', '.csv', argparse.FileType('r')))

#parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions")

parser.add_argument('--noNs', '-n', default='Y', choices=['Y','N'], help="Calculate union for fully informative positions only")


args = parser.parse_args()

table1 = open(args.matrix1.name, 'r')

#table2 = open(args.matrix2.name, 'r')

snpMatrix1 = csv.reader(table1, delimiter='\t')

#snpMatrix2 = csv.reader(table2, delimiter='\t')

textMatrix1 = [row for row in snpMatrix1]

#textMatrix2 = [row for row in snpMatrix2]

## First, parse textMatrix1, retrieve Chromosome and position, then concatenate them

idx = 0

snpPosition = []

while idx < len(textMatrix1):
    hasAmbig = 0
    if(not(re.search(pattern='^#', string=textMatrix1[idx][0]))):
        for item in textMatrix1[idx]:
            if((re.search(pattern='N', string=item)) and (args.noNs == 'Y')):
               hasAmbig = 1
        if(hasAmbig == 0):
            if(re.search(pattern='::', string=textMatrix1[idx][0])):
                newString = re.sub('::', '_', textMatrix1[idx][0])
               snpPosition.append(newString)
            else:
                snpPosition.append(textMatrix1[idx][0] + "_" + textMatrix1[idx][1])
    idx = idx + 1


for item in snpPosition:
    print(item)

        


