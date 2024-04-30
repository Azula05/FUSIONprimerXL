#!/usr/bin/env python3

"""
This script imports sequences based on the specified chromosomal positions with FASTAhack using the created index file and fasta in assets.
USE 0 BASED COUNTING!!! --> ENSEMBL 
"""

# # import all libraries
import argparse
from Bio import SeqIO
import math
import pybedtools
import re

# # get info on BSJ seq
parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-n', nargs=1, required=True, help='the nr of nucleotides surrounding the BP at each side', metavar='length')
parser.add_argument('-i', nargs=1, required=True, help='input fusionRNA bed file, 0-based')
parser.add_argument('-s', choices=['yes','no'], nargs=1, required=True, help='spliced? yes or no')
parser.add_argument('-e', nargs=1, required=True, help='bed file with canonical exons')
parser.add_argument('-t', nargs=1, required=True, help='list with canonical transcripts')


args = parser.parse_args()
temp_length = int(args.n[0])
input_bed = open(args.i[0])
splice = args.s[0]

# # read first (and only) line of input bed file
fusionRNA = input_bed.read()

# retrieve chrom start end info
fusion_chrom = fusionRNA.split()[0]
fusion_start = int(fusionRNA.split()[1]) #+ 1 # change to 1-based system
fusion_end = int(fusionRNA.split()[2])
fusion_ID = fusionRNA.split()[3]


# add parameters for circ annotation
annotation = open('annotation_' + circ_ID + '.txt', 'w')
annotation.write(circ_ID + "**")

