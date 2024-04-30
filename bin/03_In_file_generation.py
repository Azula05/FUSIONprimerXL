#!/usr/bin/env python3

############################################################################################################
############################## 1. Importing necessary libraries ############################################
############################################################################################################

import argparse

############################################################################################################
################################### 2. Parse the arguments for primer design ###############################
############################################################################################################
# note: default values are specified in the nextflow pipeline (params) and can be overwritten in the 
# nextflow.config file file

## Arguments:
parser = argparse.ArgumentParser(description='Giving all arguments to the FUSIONprimerXL pipeline scripts')
# input sequence file
parser.add_argument('-i','--input_seq',nargs=1, required=True, help='path to the file containing the fusion sequences')
# nr of primers
parser.add_argument('-p','--primer3_nr',nargs=1, type=int, required=True, help='number of primers to be designed by primer3')
# the minimum number of base pairs between the 3' ends of any two left primers
parser.add_argument('-n','--primer3_diff',nargs=1, type=int, required=True, help='nr of difference between primers (primer3: primer3_diff)')
# minimum melt temperature of the primers
parser.add_argument('-a','--min_tm',nargs=1, type= int, required=True, help='min TM (primer3: PRIMER_MIN_TM, default: 58)')
# optimal melt temperature of the primers
parser.add_argument('-c','--opt_tm',nargs=1, type= int, required=True,help='opt TM (primer3: PRIMER_OPT_TM, default: 59)')
# maximum melt temperature of the primers
parser.add_argument('-b','--max_tm',nargs=1, type= int, required=True, help='max TM (primer3: PRIMER_MAX_TM, default: 60)')
# maximum difference in melt temperature between the primers
parser.add_argument('-d', '--diff_tm', nargs=1, type= int, required=True, help='TM diff (primer3: PRIMER_MAX_DIFF_TM, default: 2)')
# minimum GC content of the primers
parser.add_argument('-e','--min_gc',nargs=1, type= int, required=True, help='min GC (primer3: PRIMER_MIN_GC, default: 30)')
# optimal GC content of the primers
parser.add_argument('-g','--opt_gc',nargs=1, type= int, required=True, help='opt GC (primer3: PRIMER_OPT_GC, default: 50)')
# maximum GC content of the primers
parser.add_argument('-f','--max_gc',nargs=1, type= int, required=True, help='max GC (primer3: PRIMER_MAX_GC, default: 80)')
# minimum amplicon length
parser.add_argument('-j','--amp_min',nargs=1, type= int, required=True, help='min amplicon length (primer3: PRIMER_MIN_SIZE, default: 50)')
# maximum amplicon length
parser.add_argument('-k','--amp_max',nargs=1, type= int, required=True, help='max amplicon length (primer3: PRIMER_MAX_SIZE, default: 200)')

## Parse the arguments:
args = parser.parse_args()

## Check some impossible values:
if args.amp_max[0] == 0:
    exit("Please provide a value for the maximum amplicon length (-k or --amp_max) different than 0")
if args.amp_min[0] == 0:
    exit("Please provide a value for the minimum amplicon length (-j or --amp_min) different than 0")
if args.primer3_nr[0] == 0:
    exit("Please provide a value for the number of primers to be designed by primer3 (-p or --primer3_nr) different than 0")
if args.amp_min[0] > args.amp_max[0]:
    exit("Please provide a value for the minimum amplicon length (-j or --amp_min) smaller than the maximum amplicon length (-k or --amp_max)")
if args.min_tm[0] > args.opt_tm[0] or args.min_tm[0] > args.max_tm[0] or args.opt_tm[0] > args.max_tm[0]:
    exit("Please provide values for the melt temperatures (-a or --min_tm, -c or --opt_tm, -b or --max_tm) in increasing order")
if args.min_gc[0] > args.opt_gc[0] or args.min_gc[0] > args.max_gc[0] or args.opt_gc[0] > args.max_gc[0]:
    exit("Please provide values for the GC content (-e or --min_gc, -g or --opt_gc, -f or --max_gc) in increasing order")

# get path to the input file
input_file = open(args.input_seq[0])
# get the number of primers to be designed by primer3
nr = args.primer3_nr[0]
# primer3_diff
diff = args.primer3_diff[0]

############################################################################################################
####################################### 3. processing ######################################################
############################################################################################################

## make a txt file for primer design
output = open("input_primer3.txt", "w")


# read the input file
fusion = input_file.read()
""" 
Output from the example input file:
chr1    16606   17055
chr1    17232   17742
chr1    17525   18061
"""

