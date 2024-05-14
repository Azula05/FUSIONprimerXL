#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

/*
====================================================================================================
Pipeline: FUSIONprimerXL
Description: A Nextflow pipeline for the design of fusion transcript-specific primers. 
License: MIT
Copyright (c) 2021 Ghent University
====================================================================================================
*/
/*
====================================================================================================
DEFAULT PARAMETERS (can be overwrittien in config file)
====================================================================================================
*/
// Input file
params.input_bed = "$baseDir/input/input_fusionRNAs.bed"
// Bowtie index
params.index_bowtie = "$baseDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "hg38_cdna"
// Fastahack index
params.index_fasta = "$baseDir/assets/GRCh38/index_fastahack"
params.index_fasta_name = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// Primer3 settings
params.primer_settings = "$baseDir/assets/primer3_settings.txt"
// chromosome sizes
params.chrom_file = "$baseDir/assets/GRCh38/chrom_sizes_GRCh38.txt"
// known exons
params.known_exons = "$baseDir/assets/GRCh38/Known_exons_GRCh38.111.bed"
// ENST_list
params.list_ENST = "$baseDir/assets/GRCh38/ENST_list_GRCh38.txt"
// output directory
params.output_dir = "$baseDir/output/"
// splice
params.splice = 'yes'
// Primer3 parameters:
params.primer3_diff = 1
params.primer3_nr = 20
// Primer Tm parameters:
params.min_tm = 58
params.opt_tm = 59
params.max_tm = 60
params.diff_tm = 2
// Primer GC parameters:
params.min_gc = 30
params.opt_gc = 50
params.max_gc = 80
// Amplicon length parameters:
params.amp_min = 50
params.amp_max = 0 // this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max
params.temp_l = 150

// filters
params.snp_filter = 'strict'
params.spec_filter = 'strict'
params.upfront_filter = "yes"
// SNP database
params.snp_url = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb'

// display help message
params.help = false

// required parameters
input_bed = file(params.input_bed)
chrom_file = file(params.chrom_file)
index_bowtie = file(params.index_bowtie)
index_fasta = file(params.index_fasta)
known_exons = file(params.known_exons)
list_ENST = file(params.list_ENST)


/*
====================================================================================================
HELP MESSAGE
====================================================================================================
*/

def helpMessage() {
	log.info"""
	
	Usage:
	
	The typical command for running the pipeline is as follows:
	nextflow run FUSIONprimerXL.nf -profile example 
	
	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	
	Mandatory pipeline arguments:
	--input_bed			path to input file with fusionRNAs in bed format (0-based annotation)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
	--index_fasta		path to directory that contains the fastahack genome and index file
	--index_fasta_name	the name of the fastahack genome file

	Optional pipeline arguments:
	--splice			when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
	--primer_settings	path to file with primer3 settings (see primer3 manual)
	--chrom_file		file containing all chromosome sizes (to validate bed file)
	--known_exons		bed file containing exon annotation
	--list_ENST			file containing ENST numbers of canonical transcripts or transcripts of interest (this file can also be left empty)
	--primer3_diff		the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr		the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm			minimum melt temperature of the primers (default: 58)
	--max_tm			maximum melt temperature of the primers(default: 60)
	--opt_tm			optimal melt temperature of the primers(default: 59)
	--diff_tm			maximum difference in melt temperature between the primers(default: 2)
	--min_gc			minimum GC contect of the primers (default: 30)
	--max_gc			maximum GC contect of the primers(default: 80)
	--opt_gc			optimal GC contect of the primers(default: 50)
	--amp_min			minimum amplicon length (default: 60)
	--amp_max			maximum amplicon length (default: 0)
	--temp_l			the number of nucleotides on each side of the breakpoint that will be used for the template (example 150 => template of 300 nts in total)
	--spec_filter		when set to 'strict', a maximum of 2 MM is allowed; when set to 'loose', a maximum of 5 combined MM is allowed
	--snp_filter		when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer; when set to 'off', no filter is applied
	--snp_url			when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	--upfront_filter	when set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed
	--output_dir		path to directory where the output files will be saved


	"""
}
// help message
if (params.help) {
	helpMessage()
	exit 0
}

/*
====================================================================================================
CHECK IF ALL PARAMETERS ARE PROVIDED CORRECTLY
====================================================================================================
*/
// --input_seq
if (!file(params.input_bed).exists()) {exit 1, "Input bed file not found: ${params.input_bed}"}
// --out_dir
if (!file(params.output_dir).exists()) {exit 1, "Output directory not found: ${params.output_dir}"}
// --index_bowtie
if (!file(index_bowtie).exists()) {exit 1, "Index file not found: ${index_bowtie}"}
// --primer_settings
if (!file(params.primer_settings).exists()) {exit 1, "Primer3 settings file not found: ${params.primer_settings}"}
// --chrom_file
if (!file(chrom_file).exists()) {exit 1, "Chromosome file not found: ${chrom_file}"}
// -- splice
if (params.splice != "yes" && params.splice != 'no'){
	exit 1, "Invalid splicing option: ${params.splice}. Valid options: 'yes','no'."}
if (params.splice == "yes"){if (!file(params.known_exons).exists()) {exit 1, "Known exons file not found: ${params.known_exons}"}}
if (params.splice == "yes" && !(params.list_ENST == "none")) {if (!file(params.list_ENST).exists()) {exit 1, "ENST list file not found: ${params.list_ENST}"}}
// --primer3_diff
if (!params.primer3_diff.toString().isNumber()){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
if (params.primer3_diff.toInteger() < 0){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
// --primer3_nr
if (!params.primer3_nr.toString().isNumber()){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
// --temp_l
if (!params.temp_l.toString().isNumber()){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}
if (params.temp_l.toInteger() < 50 || params.temp_l.toInteger() > 500){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}
// --snp_filter
if (params.snp_filter != "strict" && params.snp_filter != 'loose' && params.snp_filter != 'off'){
	exit 1, "Invalid SNP filter: ${params.snp_filter}. Valid options: 'strict','loose'."}
// --spec_filter
if (params.spec_filter != "strict" && params.spec_filter != 'loose'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}
// --upfront_filter
if (params.upfront_filter != "yes" && params.upfront_filter != 'str' && params.upfront_filter != 'snp' && params.upfront_filter != 'no'){
	exit 1, "Invalid SNP filter: ${params.upfront_filter}. Valid options: 'yes','str','snp','no'."}
//	--min_tm
if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
//	--max_tm
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
//	--opt_tm
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
//	--diff_tm
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
//	--min_gc
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
//	--max_gc
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
//	--opt_gc
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
//	--amp_min
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
//	--amp_max
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer >= 0, input 0 to adjust depending on temp_l."}

//	checking logic
if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}

/*
====================================================================================================
RUN INFO
====================================================================================================
*/
log.info """\
==============================================
FUSIONprimerXL pipeline
==============================================
OncoRNALab - Marieke Vromman / Arne Blom
Github -
Docker - 
==============================================
your input file: ${params.input_bed}
your output directory: ${params.output_dir}
"""
/*
====================================================================================================
PROCESS 1 - splitting input bed file with fusionRNAs
====================================================================================================
*/
// channels


// process
process split_fusionRNAs {
	tag "split_fusionRNAs"
	input:
	path(input_bed_handle)

	output:
	path 'fusion*'
	path 'start_time.txt'

	"""
	01_split_fusionRNAs.py -i $input_bed_handle
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""
}

/*
====================================================================================================
THE WORKFLOW
====================================================================================================
*/
workflow {
	split_fusionRNAs(input_bed)
}
