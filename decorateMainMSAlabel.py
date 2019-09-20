#!/usr/bin/python

import sys
import os.path
import argparse
import re
import csv
from io import StringIO
import pandas as pd
import numpy
#import statistics

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
    def extension(filename):
        if not filename.lower().endswith(expected_ext):
            raise ValueError()
        return openner(filename)
    return extension

## Function: Filename extractor from filepath
def getIsolateID(filePathString):
    splitStr = re.split(pattern='/', string=filePathString)
    fileNameIdx = len(splitStr) - 1
    isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx]) 
    if(len(isolateString[0]) < 10):
        isolateString = re.split(pattern='\.', string=splitStr[0])
    return isolateString[0]


parser = argparse.ArgumentParser(description='Create a subset from SNP matrix based upon input list', usage="decorateMainMSAlabel.py filepath1/MSA.fasta filepath2/isolateMetadata.tsv > output.tsv")

## read MSA file in fasta format
parser.add_argument("file1", type=ext_check('.fasta', argparse.FileType('r')))

## read CSID and metadata table file
parser.add_argument("file2", type=ext_check('.tsv', argparse.FileType('r')))

args = parser.parse_args()

## create a filehandle for fasta MSA
msaHandle = open(args.file1.name)

## create a table from CSID metadata using pandas
df = pd.read_csv(open(args.file2.name), sep='\t')

idx = 0

tempIDs = df['csid'].tolist()

## loop for appending CSID to tempid
for line in msaHandle:
    if(re.search(r'^>', string=line)):
        #print all but first character of fasta header line and chomp newline
        searchID = line[1:len(line)].rstrip()
        #findID = searchID.split('_')
        #print(findID[0])
        j = 0
        for item in tempIDs:
            if(item == searchID):
                print(">" + item + "_" + str(df.at[j, 'shortID']))
            j = j + 1
    elif(re.search(r'^[A|T|G|C|N]', string=line, flags=re.IGNORECASE)):
        print(line)
    idx = idx + 1


### Exploring pandas dataframe access operations
#print(list(df.columns.values))
#print(df[:4])
#print(df.at[1, 'csid'])
#print(len(df['tempid'].tolist()))

