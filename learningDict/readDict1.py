import sys
import os.path
import argparse
import re
import csv

dd = {}

with open('/scicomp/home/ydn3/SNP_Wrap_Python/testing/threeColumnDict.tsv', 'rU') as treeIDs:
    for line in treeIDs:
        if line.strip():
            key, value1, value2 = re.split(r'\t+', line)
            dd[key] = value1 + ", " + value2.rstrip()

print(dd)
