#!/usr/bin/python3
"""
This script filters primers based on specificity, SNPs, folding template and folding amplicon.
It takes the following arguments:
-A: input file with fusionRNA + ID (input_filter_fusionID.bed)
-P: intput file with all primers
-l: the nr of nucleotides surrounding the BP at each side
-s: the output snp file
-t: the output folding file from the template
-a: the output folding file from the amplicon
-b: the output specificity file
-p: specificity filter (strict, loose)
-f: filtering SNPs, should be strict or loose (off, strict, loose)
And outputs a file with filtered primers and a log file with the results of the filtering proces
The log file contains the following columns:
fusion_ID, chrom1, end, chrom2, start, design, primer_found, total_primers, passed, failed_spec, failed_SNP, failed_str_temp, failed_str_amp
"""

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-A", nargs = 1, required=True, help='input file with fusionRNA + ID (input_filter_fusionID.bed)')
parser.add_argument("-P", nargs = 1, required=True, help='intput file with all primers')
parser.add_argument('-l', nargs=1, required=True, help='the nr of nucleotides surrounding the BP at each side')
parser.add_argument('-s', nargs=1, required=True, help='the output snp file')
parser.add_argument('-t', nargs=1, required=True, help='the output folding file from the template')
parser.add_argument('-a', nargs=1, required=True, help='the output folding file from the amplicon')
parser.add_argument('-b', nargs=1, required=True, help='the output specificity file')
parser.add_argument('-p', nargs=1,choices=['strict','loose'], required=True, help='specificity filter (strict, loose)')
parser.add_argument('-f', nargs=1,choices=['off','strict','loose'], required=True, help='filtering SNPs, should be strict or loose (off, strict, loose)')

args = parser.parse_args()

all_fusion = open(args.A[0])
length = int(args.l[0])
SNP_filter = args.f[0]
spec_filter = args.p[0]


# get general info fusionRNA

fusionRNA = all_fusion.read()

chrom1 = fusionRNA.split()[0]
end = int(fusionRNA.split()[1])
chrom2 = fusionRNA.split()[2]
start = int(fusionRNA.split()[3])
fusion_ID = fusionRNA.split()[4]



# get specificity info from file
spec = open(args.b[0])
avoid_spec = []
all_lines = []

for line in spec:
	all_lines.append(line)
spec.close()


# rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6
for i in range(0, len(all_lines) - 1, 2):
	fwd_spec = all_lines[i].split()
	rev_spec = all_lines[i+1].split()
	if fusion_ID == all_lines[i].split("_")[0]:
		# get the number of mismatches
		fwd_MM, rev_MM = 0, 0
		if len(fwd_spec) > 7:
			fwd_MM = fwd_spec[7].count('>')
			print(fwd_spec)
		if len(rev_spec) > 7:
			rev_MM = rev_spec[7].count('>')
			print(rev_spec)
		# check if the number of mismatches is below the threshold
		if fwd_MM > 0 or rev_MM > 0:
			if spec_filter == 'loose':
				if fwd_MM + rev_MM > 5:
					avoid_spec.append(all_lines[i].split("_")[2])
			if spec_filter == 'strict':
				if (fwd_MM + rev_MM > 3) or (fwd_MM == 1 & rev_MM == 2) or (fwd_MM == 2 & rev_MM == 1):
					avoid_spec.append(all_lines[i].split("_")[2])
		"""else:
			avoid_spec.append(all_lines[i].split("_")[2])"""


print(avoid_spec)



# get folding template from file
fold_temp = open(args.t[0])

skip = fold_temp.readline()
skip = fold_temp.readline()
skip = fold_temp.readline()
fold_temp_avoid = fold_temp.readline()

# make a list from the positions
if fold_temp_avoid != '[]':
	fold_temp_avoid = fold_temp_avoid.replace('[', "").replace(']', '').split(', ')
	fold_temp_avoid = [ int(x) for x in fold_temp_avoid ]


# snps
SNPs = open(args.s[0]).readline()

SNP_avoid = 'none'

# make a list from the positions
if SNPs != '[]':
	SNPs = SNPs.replace('[', "").replace(']', '').split(', ')
	SNP_avoid = [ int(x) for x in SNPs ]


# get amplicon folding info into dicts
amp_fold = open(args.a[0])
amp_fold_avoid = {}
amp_fold_avoid_pos = {}

for lin in amp_fold:
	primer_ID = str(lin.split()[0].split("_")[1])
	deltag = float(lin.rstrip().split()[1])
	str_pos = lin.rstrip().split('\t')[2]
	amp_fold_avoid[primer_ID] = deltag
	amp_fold_avoid_pos[primer_ID] = str_pos
	
amp_fold.close()


# filter all primers and write to file

all_primers = open(args.P[0])
primer_file = open("all_primers/filtered_primers_" + fusion_ID + '_' + chrom1 + '_' + str(end) + '_' + chrom2 + '_' + str(start) + "_.txt", "a")

# to make log file
total_primers = 0
passed = 0
failed_spec = 0
failed_SNP = 0
failed_str_temp = 0
failed_str_amp = 0
design = 1
primer_found = 0

# check if file is empty
if os.path.getsize(args.P[0]) == 0:
	design = 0

# loop over all primers
for primer in all_primers:
	# add 1 to primer count
	total_primers += 1
	# get info from primer
	primer_ID = str(primer.split()[1])
	FWD_pos = int(primer.split()[4])
	FWD_len = int(primer.split()[5])
	REV_pos = int(primer.split()[6])
	REV_len = int(primer.split()[7])

	# create ranges from the forward and reverse positions
	FWD_poss = range(FWD_pos, FWD_pos + FWD_len)
	REV_poss = range(REV_pos - REV_len, REV_pos)
	filter_str = ""

	# check specificity
	if primer_ID in avoid_spec:
		filter_str = filter_str + "FAIL_specificity_"

	# check snps
	
	if (SNP_filter != 'off') & (SNP_avoid != 'none'):
		if SNP_filter == 'strict':
			if any(x in FWD_poss for x in SNP_avoid):
				filter_str = filter_str + "FAIL_snp_F_"
			if any(x in REV_poss for x in SNP_avoid):
				filter_str = filter_str + 'FAIL_snp_R_'
			

		if SNP_filter == 'loose':
			
			FWD_poss_5end = FWD_poss[0:int(len(FWD_poss)/2)]
			REV_poss_5end = REV_poss[int(len(REV_poss)/2):]

			if any(x in FWD_poss_5end for x in SNP_avoid):
				filter_str = filter_str + "FAIL_snp_F_"

			if any(x in REV_poss_5end for x in SNP_avoid):
				filter_str = filter_str + 'FAIL_snp_R_'

	# check folding template
	if any(x in FWD_poss for x in fold_temp_avoid):
		filter_str = filter_str + 'FAIL_fold_template_F_'
	if any(x in REV_poss for x in fold_temp_avoid):
		filter_str = filter_str + 'FAIL_fold_template_R_'

	# check folding amplicon

	fold_amp_avoid = amp_fold_avoid_pos["primer" + str(primer_ID)]

	if fold_amp_avoid != '[]':
		# generate list
		fold_amp_avoid = fold_amp_avoid.replace('[', "").replace(']', '').split(', ')
		fold_amp_avoid = [ int(x) for x in fold_amp_avoid ]

	# deltaG check (-9 is a good value)
	if amp_fold_avoid["primer" + str(primer_ID)] < -15:
		filter_str = filter_str + 'FAIL_fold_amplicon_'
	elif (amp_fold_avoid["primer" + str(primer_ID)] < -5) & any(x in range(FWD_len + 1) for x in fold_amp_avoid):
		filter_str = filter_str + 'FAIL_fold_amplicon_'
	elif (amp_fold_avoid["primer" + str(primer_ID)] < -5) & any(x in range(REV_pos - REV_len - FWD_pos, REV_pos - FWD_pos + 1) for x in fold_amp_avoid):
		filter_str = filter_str + 'FAIL_fold_amplicon_'



	# if all tests succeded => primer pair passed filters
	if filter_str == "":
		filter_str = "PASS_"
		primer_found = 1

	# gather info to add to log file
	if filter_str == "PASS_":
		passed += 1
	if "FAIL_specificity" in filter_str:
		failed_spec += 1
	if "FAIL_snp" in filter_str:
		failed_SNP += 1
	if "FAIL_fold_template" in filter_str:
		failed_str_temp += 1
	if "FAIL_fold_amplicon" in filter_str:
		failed_str_amp += 1


	primer_file.write(primer.rstrip() + "\t" + filter_str[0:len(filter_str)-1] + "\n")
# end of primer loop

primer_file.close()
all_primers.close()

# print log file
log_file = open("log_file_" + fusion_ID + ".txt", "w")
log_file.write("fusion_ID\tchrom1\tend\tchrom2\tstart\tdesign\tprimer_found\ttotal_primers\tpassed\tfailed_spec\tfailed_SNP\tfailed_str_temp\tfailed_str_amp\n")
log_file.close()
log_file = open("log_file_" + fusion_ID + ".txt", "a")
log_file.write(fusion_ID + '\t' + chrom1 + '\t' + str(end) + '\t' + chrom2 + '\t' + str(start) + '\t' + 
	str(design) + "\t" + str(primer_found) + "\t" + str(total_primers) + "\t" + 
	str(passed) + '\t' + str(failed_spec) + "\t" + str(failed_SNP) + "\t" + 
	str(failed_str_temp) + "\t" + str(failed_str_amp) + '\n' )
log_file.close()
