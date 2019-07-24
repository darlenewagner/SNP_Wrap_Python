#!/usr/bin/python

import sys
import os.path
import argparse
import re
import logging
import warnings
import csv
import glob

## Function: Test for readable directory
def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

## Function: Filename extractor from filepath
def getIsolateID(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
	if(len(isolateString[0]) < 10):
		isolateString = re.split(pattern='\.', string=splitStr[0])
	return isolateString[0]


logger = logging.getLogger("mummerPairSNPs.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='reads .fa files in input folder', usage="catConsenseGeneflow.py --inDir input_folder")

parser.add_argument('--inDir', type=readable_dir, action='store')

args = parser.parse_args()

from os import listdir
from os.path import isfile, join

#print(args.inDir)

# first option, list comprehension
#onlyFiles = [f for f in listdir(args.inDir) if isfile(join(args.inDir, f))]\

# second option, glob

onlyFiles = glob.glob(args.inDir + "*snp-consensus.fa")

#print(onlyFiles[0])

fastaLines = []

## loop with a Python iterator for each input file
## for an iterator, each time same operation is performed, "next" result is given

for fas in onlyFiles:
	with open (fas, 'rt') as my_file_handle:
		for fasLine in my_file_handle:
			fastaLines.append(fasLine)
	## special print, cat multifasta without internal '>'
	## print(fastaLines[0] + fastaLines[1])
	print(">" + getIsolateID(fas))
	genomeStr = ''
	for l in fastaLines:
		if(re.match(r'^(A|T|C|G|N)+', l, re.IGNORECASE)):
			genomeStr = genomeStr + l.rstrip()
	print(genomeStr)
	fastaLines = []
	

