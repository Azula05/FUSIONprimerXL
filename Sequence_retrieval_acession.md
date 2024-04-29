# Sequence retrieval by accession number
## Description
This script imports sequences of specified accession numbers from the NCBI database, then uses the end position of the first sequence and the start position of the second sequence specified by the user. The script uses 1 based counting as found on NCBI, meaning that the first nucleotide will be 1 instead of 0. The program will use the defined  breakpoint and take 150 nucleotides from each side of the breakpoint. If no 150 nucleotides are available the program will take as many nucleotides as possible. It then joins the parts of these sequences separated by * which indicates the breakpoint.  These sequences are then written to a file in a defined output directory to a text file called input_seq.txt, which can than be fed to other processes in the pipeline for processing. The output file contains a header to check some general information about the sequences that where combined and the fusion-sequence itself.

## Format of the input
The script requires 2 accessions per fusion, the breakpoint and the users email which is used by NCBI to contact you if excessive use would harm there servers performance.
The provided input file should be a csv file with the following format:
email address 
Accession 1, Accession2, End of sequence 1, Start of sequence 2

**Note**: The start and end positions are 1 based counting, meaning that the first position is 1 and the last position is the length of the sequence

**Example**:
```
A.N.Other@example.com
NM_000546, NM_000492, 153, 10
NM_001371558, NM_000518, 60, 601
```
**Output:**
```
> DEFENITION ACCESSION 1: Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA; SOURCE: Homo sapiens (human)
> DEFENITION ACCESSION 2: Homo sapiens CF transmembrane conductance regulator (CFTR), mRNA; SOURCE: Homo sapiens (human)
> Fusion_NM_000546_NM_000492
AAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCC*CTTTGGCATTAGGAGCTTGAGCCCAGACGGCCCTAGCAGGGACCCCAGCGCCCGAGAGACCATGCAGAGGTCGCCTCTGGAAAAGGCCAGCGTTGTCTCCAAACTTTTTTTCAGCTGGACCAGACCAATTTTGAGGAAAGGATACAGACA
> DEFENITION ACCESSION 1: Homo sapiens C-X-C motif chemokine ligand 13 (CXCL13), transcript variant 2, mRNA; SOURCE: Homo sapiens (human)
> DEFENITION ACCESSION 2: Homo sapiens hemoglobin subunit beta (HBB), mRNA; SOURCE: Homo sapiens (human)
> Fusion_NM_001371558_NM_000518
ACAGCCTGGACTCAGAGCTCAAGTCTGAACTCTACCTCCAGACAGAATGAAGTTCATCTC*TAATAAAAAACATTTATTTTCATTGCAA
```
## Executing
```bash
python ./00_Sequence_retrieval_accession.py -c <path to csv input> -o <path to output location>
```
Which should produce the following terminal output:
```
Retrieved sequence for: NM_000546 and NM_000492
Retrieved sequence for: NM_001371558 and NM_000518
```

## Possible errors
The path to the input file is incorrect or the file does not exist.
```
Error: Input file does not exist
```
The provided input file is not a csv file:
```
Error: Input file is not a csv file
```
The output directory does not exist
```
Error: Output directory does not exist
```
Start position 0:
```
Error: Start position cannot be 0, please use 1 based counting
```
Number of end and start positions do not match
```
Error: The number of start and end positions are not equal 1 start positions and 2 end positions were provided
```
Incorrect accession number
```
Error: Accession could not be found not found
Tried to access: NL_03
At NCBI database nucleotide
Please change the accession number and try again
```
Start position of the second sequence is larger than the sequence length
```
Error: The start position is larger than the sequence length for NM_000518
```

