#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess

## Reads an arbitrary number of snp.tsv files to compute union of SNPs positions
## The seven columns should have the following headers:
## CHROM   POS     REF     ALT     QUAL    INFO-DP CONSENS

## Function: Union 

def Union(list1, list2):
    return list(set().union(list1, list2))

for argument in sys.argv[1:]:
    columns = defaultdict(list)
    header = []

    with open(argument.name, 'r') as table:
        csvTable = csv.DictReader(table1, delimiter="\t")
        header = list(list(csvTable1)[0].keys())
        table.seek(0)
        for row in csvTable:
            for (key, val) in row.items():
                columns[key].append(val)

