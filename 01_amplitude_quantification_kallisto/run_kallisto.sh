#!/bin/env bash
##################################################################################################
##
##  name:run_kallisto.sh
##
## Purpose:
## This script the following commands:
## - Index of input cDNA file for pseudoalignment.
## - Pseudoaligmment of individual single-end FASTQ files against
##   kallisto-indexed cDNA sequences.  
##
## The script takes compressed (gzip) FASTQ files (one per sample) and results in the generation of 
## three files per sample. 
## 1) abundances.h5: HDF5 binary file containing run info, abundance estimates, bootstrap estimates,
##    and transcript length information.
## 2) abundances.tsv: Plaintext format of abundance estimates. Does not contain boostrap estimates. 
## 3) run_info.json: A json file containing information on the run.  
##
##################################################################################################
##================================================================================================
## REMOVE THIS FOR PUBLICATION REPOSITORY?
##================================================================================================
#Filenames are like 
#2016-Bter-CTH-C60-4-W1-head.D711_D506.HM3WHBBXX.s_6.R1.fastq.gz
#2016-Bter-CYH-C03-4-Q-head.D711_D503.HM3WHBBXX.s_6.R1.fastq.gz 

## Create directory on scratch:
mkdir ~/scratch/2018-09-19-kallisto_actualworkers
## Soft Link 
ln -s ~/scratch/2018-09-19-kallisto_actualworkers tmp
##
# Note to self. to create an array the folowing is preferred;
#    fastqs=(input/fastq_all/2016-Bter-${treatment}*-W[1-4]-*gz)   but not using here bc i need ot use them after
##================================================================================================

## The first step involves defining treatments, including control and pesticide treatments:
treatments=(CON CLO IMI)

## Define reference file:
kallisto_ref=tmp/Bter1_cdna

## Generate kallisto-index:
./kallisto index -i ${kallisto_ref} input/reference/Bter1_cdna.fa 

## For quantifying transcript abundance using Kallisto, generate a shell script to 
## run the same kallisto quantification for each sample.  
for treatment in $treatments; do
  echo "running with ${treatment}"
  fastqs=`ls input/fastq_all/2016-Bter-${treatment}[-_]*-W[0-9]-*gz`
  colonies=$(echo $fastqs | cut -d '-' -f 4 | sort -u)
  colonies=($(echo "$colonies"))

  mkdir -p tmp/quant
  for colony in $colonies; do
    echo "running with ${colony}"
    fastqs_for_colony=`echo $fastqs | grep "\-${colony}\-" | tr "\n" " "`
    output=tmp/quant/${treatment}-${colony}
    echo "./kallisto quant -i ${kallisto_ref}  --output-dir=${output} --single --fragment-length=300 --sd=20  --threads=35 ${fastqs_for_colony}" >> ./kallisto_commands.sh
  done
done

## Run shell script to estimate transcript abundance counts for each sample per treatment:  
sh kallisto_commands.sh > kallisto_commands.log  2> kallisto_commands.err
mv tmp/quant results/kallisto_quantifications

