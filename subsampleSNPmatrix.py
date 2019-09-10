#!/usr/bin/python

import sys
import os.path
import argparse
import re
import csv
from io import StringIO
import pandas as pd
import numpy
import statistics

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


parser = argparse.ArgumentParser(description='Create a subset from SNP matrix based upon input list', usage="subsampleSNPmatrix.py  filepath1/SNPmatrix.tsv filepath2/isolateList.txt > output.tsv")

## read SNP matrix file name
parser.add_argument("file1", type=ext_check('.tsv', argparse.FileType('r')))

## read isolate list file name
parser.add_argument("file2", type=ext_check('.txt', argparse.FileType('r')))

args = parser.parse_args()

## open SNP matrix file
table = numpy.loadtxt(args.file1.name, skiprows=1, usecols=range(1,112)) 

with open(args.file1.name) as f:
    first_line = f.readline()

#df = pd.read_csv(args.file1.name, sep='\t')
#arr = numpy.array(df)

## open isolate list file
filehandle2 = open(args.file2.name)

clade = [line.strip() for line in filehandle2]

#print(clade)

clade.sort()

IDlist = first_line.split("\t")
IDlist.pop(0)

#print(IDlist[0])
#print(table[2,1])

# Xs and Ys store indices from original SNPs matrix
# where SNP IDlist[i] == isolateList[j]
Xs = []
Ys = []

# Populate Xs and Ys with SNP IDlist[i] == isolateList[j]
for i in range(0, len(IDlist)):
        for j in range(0, len(clade)):
                if(IDlist[i] == clade[j]):
                        #print(IDlist[i])
                        Xs.append(i)
                        Ys.append(i)


print(".", end="")

for h in range(len(Xs)):
        print("\t", IDlist[Xs[h]], end="")
print()

snpsList = []

## print out a new SNPs matrix based upon isolateList
for ii in range(0, len(Xs)):
        rowStart = 0;
        zeroOut = 0
        for jj in range(0, len(Ys)):
                if(rowStart == 0):
                        print(IDlist[Xs[ii]], end="")
                        print("\t", table[Xs[ii], Ys[jj]], end="")
                        rowStart = 1
                        if(table[Xs[ii], Ys[jj]] == 0):
                                zeroOut = 1
                        else:
                                snpsList.append(table[Xs[ii], Ys[jj]])
                else:
                        if(zeroOut == 0):
                                print("\t", table[Xs[ii], Ys[jj]], end="")
                        else:
                                print("\t0.0", end="")
                        if(table[Xs[ii], Ys[jj]] == 0):
                                zeroOut = 1
                        elif(zeroOut == 0):
                                snpsList.append(table[Xs[ii], Ys[jj]])
        print()


print("Cluster Median SNPs:\t",statistics.median(snpsList))

