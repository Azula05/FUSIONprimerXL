#!/usr/bin/env python3

import argparse
import os
import datetime
# This script is used to split up the input file into different fusion sequences
# The sequence can eighter be provided in a text file or can be retrieved

####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
parser = argparse.ArgumentParser(description='Give input to the FUSIONprimerXL pipeline')
# input fusion seq file
parser.add_argument('-t', '--input_type', nargs=1, choices=['file','retrieve'], required=True, help='input type: file or retrieve')
parser.add_argument('-i', '--input_seq', nargs=1, required=True, help='input fusion seq file')
parser.add_argument('-o', '--output_dir', nargs=1, required=True, help='output directory')
args = parser.parse_args()

####################################################################################################
####################################  2. Retrieve sequence  ########################################
####################################################################################################
if args.input_type[0] == 'retrieve':





####################################################################################################
####################################  3. Processing  ###############################################
####################################################################################################
# save current working directory
Base = os.getcwd()
# store the path to the output folder with current date and time
currentDT = datetime.datetime.now()
output_folder = os.path.join(args.output_dir[0], 'fusion_files_' + currentDT.strftime("%Y-%m-%d_%H-%M") + '/')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
print(output_folder)

## Input file open
input_file = open(args.input_seq[0])

# Count the number of lines in the input file (equal to fusion transcripts)
count_lines = 0
for line in input_file:
	count_lines +=1

## close input file
input_file.close()

## Create the format for the fusion file
fusion_nr = len(str(count_lines))
fusion_nr = "fusion{:0" + str(fusion_nr) +"d}"

##  Open the input file again
input_file = open(args.input_seq[0])

# move to the output folder
os.chdir(output_folder)

## Create all fusion file and an epmty list to check for duplicates
all_fusion_file = open('all_fusion.txt', 'w')
all_fusion = []

# Loop over the input file and create individual files for each fusion
ID = 0
for fusion in input_file:
    # Create the ID for the fusion file
    ID_str = fusion_nr.format(ID)
    # apppend line to list
    all_fusion.append(fusion)
    # Create individual fusion file
    ind_fusion_file = open(ID_str + ".txt", "w")
    ind_fusion_file.write(fusion + '\t' + ID_str + '\n')
    ind_fusion_file.close()
    # Append fusion to all fusion file
    all_fusion_file.write(ID_str + '\t' + fusion)
    
    ID += 1

# check for duplicates
def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

if checkIfDuplicates(all_fusion):
	raise SystemExit('One or more fusion sequences are present more than once in your input file. Please make sure each fusion sequence is unique.')

## Close all files
all_fusion_file.close()
input_file.close()

# move back to the base directory
os.chdir(Base)
