#!/usr/bin/env python3

# This script imports sequences of 2 genes from a text file
# uses the start and end position to get part of the sequence
# and writes the sequence to a new file to be used downstream

"""
Concept of the file to provide
Gene1	Gene2	Start1	End1	Start2	End2
Example:
10563	10564	100	200	300	400

Line added to the file
Part Gene1 * Part Gene2
With * representing the breakpoint
"""

from Bio import Entrez
from Bio import SeqIO
import os

gene_ids = ["10563", "3043"]
ID=0

# Always tell NCBI who you are (your email address)
Entrez.email = "arne.blom@hotmail.com"

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
    print(sequence,"\n\n")
    input_file.close()
    # remove the file
    os.remove("sequence.gb")


####################################################################################################################
# THE CODE
""" 
# Get the sequence of a gene from NCBI
# Replace with your desired gene ID
gene_id = "10563"



# Use Entrez.efetch to retrieve the sequence in GenBank format
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
print(sequence)
input_file.close()

# remove the file
os.remove("sequence.gb") """





