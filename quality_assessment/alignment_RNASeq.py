#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: alignment_RNASeq.py
#
# Date: 12/08/2017
#
##############################################################################
# Import modules
import os.path

from helper_functions import *

# This script takes one fastq files (i.e. single end) per sample(s) as input and aligns
#  input reads against indexed reference genome. It then converts sam file to BAM format,
# sorts and indexed the BAM file.

# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule index_genome:
#    - A single fasta file containing genome of interest is provided as input.
#      Input fasta file is indexed by STAR.
# rule align_to_genome:
#   - For each sample, sequences are aligned in pairs to the STAR-indexed genome.
#    The rule checks if sequence headers for each input pair contain the same header information.
#    If they differ, an error is raised.
# rule sort_sam_to_bam:
#    - For each sample, aligned SAM file is converted to BAM file and sorted.
# rule index_bam:
#    - For each sample, index sorted BAM file.
# rule mapping_stats:
#    - For each sample, calculate mapping statistics for reads.

##############################################################################
# Sample information
##############################################################################
# Bumblebee (Bombus terrestris) workers and queens were purchased from a commercial
# supplier (Agralan).
# Total of 87 samples were sequenced using NextSeq500 (1*75bp)
# All workers were age-controlled to 10 days and were exposed to treatment for
# a period of 4 days. Age of queens is unknown.
# Total RNA was extracted from the head of each individual.
# Library preparation was performed within the TruSeq stranded mRNA library preparation kit.
# Individuals were collected from individual colonies and were randomly assigned to
# either control group or experimental group (consisted of 15 different pesticide treatments)
# Four individuals were collected per treatment.
# Each individual and treatment were assigned unique identifiers
# Example: '2016-Bter-DIN-C47-4-W2-head'
# Explanation:{year_collected}_{species}_{treatment}_{caste}_{colony_number}_{timepoint}_{individual_ID}_{tissue_type}
# species: Bter = Bombus terrestris
# caste: W = Worker; Q=Queen.

##############################################################################
# Prior to use
##############################################################################
# To run alignment_to_coverage.py:
#  1. Download and install the following software:
#   STAR
#   samtools
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
#
#  3. Assign global variables for use within specific rules
#     Please see section below for further information on variable to be assigned
#
#  4. Assignment of wildcards to be used within the rules
#
#  5. Define all the input and outputs of each rule with respect to the above defined wildcards
#
#  6. Each rule willl take assigned input, global variables and generate defined outputs.
#
#  7. Input data should be formatted in the context of defined wildcards: {sample}.fastq
#      For example: 2016-Bter-DIN-C47-4-W2-head.combined.fastq
#
# 8. Make a text.file called 'sample_list.txt' and put in same directory as Snakefile.
#     Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       2014_Bter_P_D_14_260_head
#       2014_Blap_P_D_21_412_thorax
#       2014_Bpas_A_D_20_391_thorax

##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################
# Specify the reference genome file to be indexed
GENOME_DIR    = "./GCF_000214255.1_Bter_1.0_genomic.fna"

# Specify the custom index of Bowtie2-build to use
GENOME_FILE   = "GCF_000214255.1_Bter_1.0_genomic.fna"

# Specify the number of threads
MAX_THREADS   = 20

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('queen_samples.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign path for raw read data for aligning
RAW_DATA          = "./{{samples}}.nextseq_hiseq.combined.fq"

# Output aligned, sorted and indexed data data in these directories, respectively
ALIGNED_DATA      = "raw_temp/01_aligned/{samples}.Aligned.out.sam"

SORTED_DATA       = "raw_temp/02_sorted/{samples}.Aligned.out.sorted.bam"

INDEXED_DATA      = "raw_temp/02_sorted/{samples}.Aligned.out.sorted.bam.bai"

# Output mapping stats here:
STATS_DATA        = "raw_temp/03_mapping_stats/{samples}.stats.txt"


##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project']     = os.path.abspath('../../../../../')
dirs['src']         = os.path.join(dirs['project'], 'bin')

# Create an empty dictionary called 'tools'
tools = {}
tools['STAR']       = os.path.join(dirs['src'], 'STAR')
tools['samtools']   = os.path.join(dirs['src'], 'samtools')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(STATS_DATA, samples=SAMPLES)

# Each rule will specify an intermediate step
# Align raw reads to indexed genome:
rule align_to_genome:
    input:  forward=expand("./{{samples}}.nextseq_hiseq.combined.fq")
    output: ALIGNED_DATA
    run:
        input.reverse=re.split('/|.n', (''.join(input.forward)))[1]
        print(input.reverse)
        #input.reverse2=re.split('com|.f', (''.join(input.reverse)))[0]
        #print(input.reverse2)
        align_list="raw_temp/01_aligned/" + input.reverse + "."
        print(align_list)
        check_files_arent_empty(input)
        shell("{tools[STAR]} --genomeDir ./ \
                             --readFilesIn {input.forward} \
                             --outFileNamePrefix {align_list} \
                             --runThreadN {MAX_THREADS} && [[ -s {output} ]]")

# Convert sam file to bam, sort and output:
rule sort_sam_to_bam:
    input:  ALIGNED_DATA
    output: SORTED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} view -bS {input} \
             | {tools[samtools]} sort - -o {output} && [[ -s {output} ]]")

## Index sorted bam file and output:
rule index_bam:
    input:  SORTED_DATA
    output: INDEXED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} index {input} > {output} && [[ -s {output} ]]")

## Calculate mapping stats:
rule mapping_stats:
    input:  SORTED_DATA
    output: STATS_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} flagstat {input} > {output} && [[ -s {output} ]]")
