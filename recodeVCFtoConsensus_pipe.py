#!/usr/bin/python

import sys
import os
import argparse
import logging
import re
import csv
import numbers
import time

def readable_dir(prospective_dir):
    if not os.path.isdir(prospective_dir):
        raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
    if os.access(prospective_dir, os.R_OK):
        if( not prospective_dir.endswith("/") ):
            prospective_dir = prospective_dir + "/"
        return prospective_dir
    else:
        raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

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

logger = logging.getLogger("recodeVCFtoConsensus.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='A wrapper for vcftools, bgzip, tabix, and bcftools consensus', usage="recodeVCFtoConsensus.py filepath/input.fasta filepath/reference.fasta --outDir outfilefolder/ --strict [Y/N]")

## input file, must be unzipped vcf
parser.add_argument("filename",type=ext_check('.vcf', argparse.FileType('r')))

## reference file, requires associated .fai
parser.add_argument("reference",type=ext_check('.fasta', argparse.FileType('r')))

## output folder
parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

## for basic filter, on/off
parser.add_argument('--filter', '-fil', default='N', choices=['Y','N'], help="Remove SNP positions which do not have 'PASS'")

## for most stringent filter
parser.add_argument('--strict', '-st', default='N', choices=['Y','N'], help="Set thin to 5 bp, --thin 5, than Geneflow default")

args = parser.parse_args()

## output folder
outFolder = args.outDir

## input vcf file
vcfRaw = args.filename.name

## reference fasta file
mapReference = args.reference.name

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

logger.info("Parameters loaded.")

origWD = os.getcwd()

## Simple filename extractor function

def getIsolateStr(filePathString, format):
    splitStr = re.split(pattern='/', string=filePathString)
    fileNameIdx = len(splitStr) - 1
    fileString = splitStr[fileNameIdx]
    if(format == 1):
        isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
        fileString = isolateString[0]
    return fileString

print(getIsolateStr(vcfRaw, 1))

outputPath = outFolder + getIsolateStr(vcfRaw, 1)

logger.info("begin vcf filtering.")
if(args.filter == 'Y'):
    os.system("vcftools --vcf {} --recode --recode-INFO-all --out {} --remove-filtered-all --remove-indels --minGQ 50".format(vcfRaw, outputPath))
elif(args.strict == 'Y'):
    os.system("vcftools --vcf {} --recode --recode-INFO-all --out {} --remove-filtered-all --remove-indels --minGQ 50 --minQ 40 --min-meanDP 10 --thin 5".format(vcfRaw, outputPath))
else:
    os.system("vcftools --vcf {} --recode --recode-INFO-all --out {} --minGQ 50".format(vcfRaw, outputPath))

logger.info("filtering complete, begin indexing.")

vcfRecode = outputPath + ".recode.vcf"

gzRecode = vcfRecode + ".gz"

os.system("bgzip -c {} > {}".format(vcfRecode, gzRecode))

time.sleep(0.75)

os.system("tabix -fp vcf {}".format(gzRecode))

logger.info("indexing complete, begin building consensus fasta.")

fastaRecode = outputPath + ".recode.fasta"

os.system("bcftools consensus -f {} {} > {}".format(mapReference, gzRecode, fastaRecode))

logger.info("consensus fasta built.")

