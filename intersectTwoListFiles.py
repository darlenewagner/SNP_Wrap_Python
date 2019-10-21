#!/usr/bin/python

import sys
import os.path
import argparse
import re
import logging
import warnings
import csv
import subprocess

## Reads two single-column .txt or .tsv files to compute intersection

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

logger = logging.getLogger("intersectTwoListFiles.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions from two plain text, single-column files', usage="intersectTwoListFiles.py listFile1.tsv listFile2.tsv")

parser.add_argument("listFile1", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument("listFile2", type=tsv_check('.txt', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions")

args = parser.parse_args()

table1 = open(args.listFile1.name, 'r')

table2 = open(args.listFile2.name, 'r')

## list comprehension followed by map(lambda), map returns an iterator
textList1 = [row for row in table1]
textList1 = map(lambda s: s.strip(), textList1)

textList2 = [row for row in table2]
textList2 = map(lambda s: s.strip(), textList2)

#for item in textList1:
#    print(item)
if(args.union == 'Y'):
    allPositions = Union(textList1, textList2)
else:
    allPositions = Intersection(textList1, textList2)

for item in allPositions:
    print(item)



