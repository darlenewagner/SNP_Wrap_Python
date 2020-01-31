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


parser = argparse.ArgumentParser(description='Create a subset from SNP matrix based upon input list', usage="decorateCaurisMSAlabel.py filepath1/MSA.fasta filepath2/isolateMetadata.csv > output.tsv")

## read MSA file in fasta format
parser.add_argument("file1", type=ext_check('.fasta', argparse.FileType('r')))

## read BLIND_ID and metadata table file
parser.add_argument("file2", type=ext_check('.csv', argparse.FileType('r')))

args = parser.parse_args()

## create a filehandle for fasta MSA
msaHandle = open(args.file1.name)

## create a table from BLIND_ID metadata using pandas
df = pd.read_csv(open(args.file2.name), sep=',')

idx = 0

tempIDs = df['BLIND_ID'].tolist()

for line in msaHandle:
    if(re.search(r'^>', string=line)):
        searchID = line.split("_")
        j = 0
        for item in tempIDs:
            if(item == searchID[1].rstrip()):
                print( ">" + str(df.at[j, 'Original_ID'] + " " + str(df.at[j, 'Country'])))
            j = j + 1
    elif(re.search(r'^[A|T|G|C|N]', string=line, flags=re.IGNORECASE)):
        print(line, end="")

