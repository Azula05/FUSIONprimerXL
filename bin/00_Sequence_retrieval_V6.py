#!/usr/bin/env python3

# This script imports sequences of 2 genes from a text file
# uses the start and end position to get part of the sequence
# and writes the sequence to a new file to be used downstream

"""
Concept of the file to provide
use 1 from beginning of transcipt
use 0 to end of transcript
transcript 1  transcript 2	Start T 1 End T 1	Start T 2 End T 2
Example:
NM_000546     NM_000492  1 112 210 0    (TP53,CFTR, random breakpoint)
NM_001371558  NM_000518	 1 51	 1001 1201 ( CXCL13,HBB random not a real fusion gene)
The positions are as you would normally count not first one is 0
example first 50 and last 77
Line added to the file
Part Gene1 * Part Gene2
With * representing the breakpoint
"""
# Requires biopython and os
from Bio import Entrez
from Bio import SeqIO
import os

# Hardcoded values for testing (needs to be replaced with a file input)
Accessions = ["NM_000546", "NM_000492","NM_001371558","NM_000518"]
ID = 0
counter = 0
num = 0
starts = [1, 210, 1, 100]
start = []
ends = [112, 0, 51, 0]
end = []
# Always tell NCBI who you are (your email address)
# A.N.Other@example.com
Entrez.email = "arne.blom@hotmail.com"

# Change the positions from 1 based to 0 based counting (1 in start is from beginning, 0 in end is to end)
for i in starts:
    if i == 0:
        i = 1
    else:
        i -= 1
    start.append(i)
for i in ends:
    if i == 0:
        i = None
    end.append(i)

# A control to check if the beginning and end is correct
def check_start_end(start_pos, end_pos, ID):
    if end_pos == None:
        return
    else:
        if start_pos > end_pos:
            exit("Error: Start position is greater than end position for fusion pair: " + str(ID+1))
        elif start_pos == end_pos:
            exit("Error: Start position is equal to end position for fusion pair: " + str(ID+1))

if len(starts) != len(ends):
    exit("Error: The number of start and end positions are not equal " + str(len(starts)) + " start positions and " + str(len(ends)) + " end positions were provided")


# Open the sequence output file
sequence_file = open("input_seq.txt", "w")

# Loop through the accessions
for Accession in Accessions:
    Accession = str(Accession)
    try:
        # Fetch the genbank record
        handle = Entrez.efetch(db="nucleotide", id=Accession, rettype="gb", retmode="text")
        sequence_data = handle.read()
        handle.close()
        # write the record to a temporary file
        output_file = open("sequence.gb", "w")
        output_file.write(sequence_data)
        output_file.close()
    except:
        print("Error: Accession could not be found not found")
        print("Tried to access: " + Accession + "\nAt NCBI database nucleotide\nPlease change the accession number and try again")
        exit()
    # parse the genbank record file (nucleotide database)
    input_file = open("sequence.gb", "r")
    sequence_record = SeqIO.read(input_file, "genbank")
        # get the sequence
    sequence = sequence_record.seq
    # First sequence
    if counter % 2 == 0:
        # start and end position of the sequence
        if num != 0:
            num+=1
        start_pos = start[num]
        end_pos = end[num]
        check_start_end(start_pos, end_pos, ID)
        # Slice the sequence
        sequence1 = sequence[start_pos:end_pos]
        # extract the header information
        first_Accession = Accession 
        defenition1 = sequence_record.description
        source1 = sequence_record.annotations.get("source", None)
        # write the first header line
        sequence_file.write("> DEFENITION ACCESSION 1: " + defenition1 + "; SOURCE: " + source1 + "\n")
    # Second sequence
    else:
        # start and end position of the sequence
        num+=1
        start_pos = start[num]
        end_pos = end[num]
        check_start_end(start_pos, end_pos, ID)
        # Slice the sequence
        sequence2 = sequence[start_pos:end_pos]
        # extract the header information
        defenition2 = sequence_record.description
        source2 = sequence_record.annotations.get("source", None)
        # write the second header line
        sequence_file.write("> DEFENITION ACCESSION 2: " + defenition2 + "; SOURCE: " + source2 + "\n")
        # write the third header line
        sequence_file.write("> Fusion" + str() + "_" + str(first_Accession) + "_" + str(Accession) +"\n")
        # write the sequences to the file
        sequence_file.write(str(sequence1))
        sequence_file.write("*")
        sequence_file.write(str(sequence2))
        sequence_file.write("\n")
        # Terminal progress update
        print("Retrieved sequence for: " + first_Accession +" and " + Accession)
        # Increase fusion gene ID
        ID+=1
    # close the input file (genbank record)
    input_file.close()
    # remove the temporary genbank record file
    os.remove("sequence.gb")
    counter+=1
# close the finished output file
sequence_file.close()






