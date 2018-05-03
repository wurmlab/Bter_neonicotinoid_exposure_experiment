#!/bin/env bash

## Take inputs from the command line:
input_fasta=$1
input_gtf=$2
overhang=$3

## Check arguments are provided:
if [ $# -eq 0 ]
  then
    echo "No arguments supplied. Usage: ./make_STAR.sh input.fasta input.gtf overhang"
    echo "For overhang, a value of the maximum read length minus 1 should be provided."

fi

## Create STAR indices
../../../../../bin/STAR --runThreadN 10 \
               --runMode genomeGenerate \
               --genomeDir . \
               --genomeFastaFiles "$input_fasta" \
               --sjdbGTFfile "$input_gtf" \
               --sjdbOverhang "$overhang" 
