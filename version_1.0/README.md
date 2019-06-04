Welcome to UPDio, Version 1.0

# This verison adds mulitsample support.

## Introduction 
UPDio is designed to identify uniparental disomy in probands of trio VCF data.  
	
This directory contains 2 files 
1. README.txt  
	This file
2. UPDio.pl  
	The main calling script

and, three child directories
1. scripts  
	This contains some helper scripts for UPDio to run
2. pre_processing  
	This shall serve as a location to prepare samples for analysis
3. sample_data  
	This contains a pre-processed exome trio ready to be run on UPDio
	
UPDio requires 3 input files: these are the single-sample VCF files corresponding to the proband, mom, and dad samples of a trio  
	# update: or one (child mother father) trio multisample VCF file.

UPDio also recommends including a file of CNV calls for the proband; this is recommended to limit false positives  
	# this step is especially important when detecting smaller UPD events

At the command line, type 'perldoc UPDio.pl' to familiarize yourself with how to run UPDio, then return here.


## Setup
This program is written in Perl and R  
Please ensure that dependencies are installed before attempting to run UPDio; they are all required.  
Dependencies can be downloaded from CPAN and CRAN.
	
	R Dependencies
		quantsmooth
		ggplot2
	Perl Dependencies
		Statistics::R (0.31)
		Path::Class
		Vcf
		Iterator::Simple
		List::MoreUtils
		Math::Round
		Const::Fast
	UPDio was tested using R version 2.14.1

## Pre-processing
UPDio requires two input file format requirements:
	
1. The VCF files must be sorted in numeric chromosome order: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
		We include a sort-vcf script to generate vcfs sorted in this order
2. The VCF files must include homozygous reference (homREF) genotypes, i.e. "0/0" positions
		We include a add_hom_refs_to_sorted_vcfs.pl script that can add these positions
		Note that your VCFs may already include 0/0 positions, if so, you can skip this step.

The pre_processing directory contains instructions for how to prepare input files for UPDio  
Included in this directory is a commands file 'commands_to_preprocess_files.sh' showing the commands used to sort and add homREF positions.  
	# If this doesn't work, email me and I'll make the positions file for you.  
	# SureSelect_v4 is now on github, SureSelect_v5 soon to be uploaded
	
## Running UPDio
First try running UPDio on the example trio that is supplied before running it on your own trio data  
To do this, refer to the file in 'scripts/run_UPDio_example.sh' containing the command to run UPDio on the example trio  
When this script has completed you should be able to observe a UPD event in the output files.

## Output
Output is stored by default in a directory called 'output_dir' but you can specify your own output directory as an option to UPDio

Output files suffixes
1. table  
	a tabulation of informative genotypes by chromosome
2. events_list  
	a print out of informative genotypes found in the proband
3. upd  
	a list of significant UPD events found (these lines can be long; try less -S to read this file)  
	bear in mind that CNVs often masquerade as UPD events 
4. pngs  
	the plot of UPD events
5. log  
	a log file


## Troubleshooting
Q: I'm getting an error message that looks like this: "Can't locate Statistics::R in @INC (@INC contains: ...)"  
A: Your perl dependencies are not installed; install all dependencies before running UPDio

Q: I'm getting an error message that says "cannot find set method in Statistics::R"  
A: You are using an older version of Statistics::R, please upgrade to 0.31

Q: How can I gain access to 1000 MAFs to select common sites in my exome design?  
A: Sorry, this file was too large to load to github! Solutions: https://www.biostars.org/p/6897/ or email me.
