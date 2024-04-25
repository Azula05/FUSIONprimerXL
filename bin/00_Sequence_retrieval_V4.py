#!/usr/bin/env python3

# This script imports sequences of 2 genes from a text file
# uses the start and end position to get part of the sequence
# and writes the sequence to a new file to be used downstream

"""
Concept of the file to provide
Gene1	Gene2	End gene 1	Start gene 2
Example:
7157     1080  112 210     (TP53,CFTR, random breakpoint)
10563	  3043	51	1001	(random not a real fusion gene)
The positions are as you would normally count not first one is 0
example first 50 and last 77
Line added to the file
Part Gene1 * Part Gene2
With * representing the breakpoint
"""

from Bio import Entrez
from Bio import SeqIO
import os

gene_ids = ["7157", "1080","10563","3043"]
ID = 0
counter = 0
position = [51,1001,112, 210]
positions = []
for num in position:
    num -= 1
    positions.append(num)


# Always tell NCBI who you are (your email address)
Entrez.email = "arne.blom@hotmail.com"

# Open the sequence output file
sequence_file = open("input_seq.txt", "w")

for gene_id in gene_ids:
    gene_id = gene_id
    # Fetch the genbank record
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
    sequence_data = handle.read()
    handle.close()
    # write the sequence to a file
    output_file = open("sequence.gb", "w")
    output_file.write(sequence_data)
    output_file.close()
    # Read the sequence from the file
    input_file = open("sequence.gb", "r")
    sequence_record = SeqIO.read(input_file, "genbank")
    sequence = sequence_record.seq
    cut = positions[counter]
    if counter % 2 == 0:
        sequence1 = sequence[:cut]
        first_ID = gene_id 
    else:
        sequence2 = sequence[cut:]
        sequence_file.write("> Fusion" + str() + "_" + str(first_ID) + "_" + str(gene_id) +"\n")
        sequence_file.write(str(sequence1))
        sequence_file.write("*")
        sequence_file.write(str(sequence2))
        sequence_file.write("\n")
        print("Retrieved sequence for gene ID: " + first_ID +" and " + gene_id)
        ID+=1
        
    input_file.close()
    # remove the file
    os.remove("sequence.gb")
    counter+=1

sequence_file.close()






