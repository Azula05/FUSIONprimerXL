#!/usr/bin/env nextflow

// Enamble DSL2
nextflow.enable.dsl=2

/*
====================================================================================================
Pipeline: FUSIONprimerXL
Description: A Nextflow pipeline for the design of fusion primers for the amplification 
of fusion genes.
License: MIT
Copyright (c) 2021 Ghent University
====================================================================================================
*/
/*
====================================================================================================
DEFAULT PARAMETERS (can be overwrittien in config file)
====================================================================================================
*/
// Reference and input sequence
params.index_bowtie = "$baseDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "GRCh38_dna"
params.input_seq = "example/path"
// Primer3 settings
params.primer_settings = "$baseDir/assets/primer3plus_settings.txt"
// output directory
params.output_dir = "$baseDir/output/"
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
params.amp_max = 200

params.temp_str_filter = 'on'
params.spec_filter = 'strict'
params.help = false

// Required parameters
input_seq = file(params.input_seq)
index_bowtie = file(params.index_bowtie)


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
	--input_seq		path to input file with fusion sequences (the fusion should be indicated with a *)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)


	Optional pipeline arguments:
	--primer_settings	path to file with primer3plus settings (see primer3 manual)
	--primer3_diff		the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr		the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm		minimum melt temperature of the primers (default: 58)
	--max_tm		maximum melt temperature of the primers(default: 60)
	--opt_tm		optimal melt temperature of the primers(default: 59)
	--diff_tm		maximum difference in melt temperature between the primers(default: 2)
	--min_gc		minimum GC contect of the primers (default: 30)
	--max_gc		maximum GC contect of the primers(default: 80)
	--opt_gc		optimal GC contect of the primers(default: 50)
	--amp_min		minimum amplicon length (default: 50)
	--amp_max		maximum amplicon length (default: 200)
	--temp_str_filter	when set to 'on', secunday structures in the tample are avoided during primer design; when set to 'off', this is not done (default: 'on')
	--spec_filter		when set to 'strict', only 2MM + 3MM are allowed; when set to 'loose', 2MM + 2MM and 2MM + 1MM are also allowed
	--output_dir		path to directory where the output files will be saved


	"""
}

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
if (!file(params.input_seq).exists()) {exit 1, "Input fusion sequence file not found: ${params.input_seq}"}
// --out_dir
if (!file(params.output_dir).exists()) {exit 1, "Output directory not found: ${params.output_dir}"}
// --index_bowtie
if (!file(index_bowtie).exists()) {exit 1, "Index file not found: ${index_bowtie}"}
// --primer_settings
if (!file(params.primer_settings).exists()) {exit 1, "Primer3 settings file not found: ${params.primer_settings}"}
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
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer > 0."}
//	--temp_str_filter
if (params.temp_str_filter != "on" && params.temp_str_filter != 'off'){
	exit 1, "Invalid temp_str_filter filter: ${params.sec_str_filter}. Valid options: 'on','off'."}
//	--spec_filter
if (params.spec_filter != "strict" && params.spec_filter != 'loose'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}
//	checking logic
if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}
if (params.amp_min.toInteger() > params.amp_max.toInteger() ) {exit 1, " amp_min and amp_max: amp_max (${params.amp_max}) should be > amp_min (${params.amp_min})"}

/*
====================================================================================================
RUN INFO
====================================================================================================
*/
log.info """\
==============================================
FUSIONprimerXL pipeline
==============================================
OncoRNALab - Marieke Vromman (DSL1) / Arne Blom (DSL2)
Github -
Docker -
==============================================
your input file: ${params.input_seq}
your output directory: ${params.output_dir}
"""
/*
====================================================================================================
PROCESS 1 - splitting input file
====================================================================================================
*/
// channels
input_seq_handle = Channel.fromPath(params.input_seq)

// process
process split_input{
    input:
    path ind_fusion_file_handle

    output:
	path 'fusion*', emit: ind_fusion_file
	path 'start_time.txt', emit: start_time
	path 'all_fusion.txt', emit: all_fusion

	script:	
	"""
	01_split_input.py -i $input_seq_handle
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""
}

/*
====================================================================================================
THE WORKFLOW
====================================================================================================
*/

workflow {
    // process 1
    split_input(input_seq_handle)
}