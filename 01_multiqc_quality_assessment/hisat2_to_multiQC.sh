############################################################################################
##
## Author: Joe Colgan
## Queen Mary University of London
## All rights reserved
############################################################################################

## Required: MultiQC file for HISAT2 alignments:
module load samtools/1.3.1
module load qualimap/2.2.1 
module load multiqc/0.7

## Link aligned sorted sam file:
ln -s ../2018-09-18-spliced_alignment_worker/tmp/quant/ input_workers
ln -s ../2018-11-30_spliced_alignment_queen-tmp/tmp/quant/ input_queens

## Make results directory:
mkdir results

## For each input sam file, extract mapped reads and output as bam file.
## Calculate summary stats for the new bam file. 
for name in input_workers/*sam; 
do 
new_name="$(echo "$name" | cut -d '/' -f 2 - | cut -d '.' -f 1 - )";  
samtools view -bS "$name" | samtools sort --threads 20 -o input_workers/"$new_name".sorted.bam - ; 
done

## Perform the same for queen data:
for name in input_queens/*sam;
do
new_name="$(echo "$name" | cut -d '/' -f 2 - | cut -d '.' -f 1 - )";
samtools view -bS "$name" | samtools sort --threads 20 -o input_queens/"$new_name".sorted.bam - ;
done

## Run qualimap:
for name in results/*sorted.bam; 
do 
qualimap bamqc -bam "$name"; 
done

## Run multiqc:
multiqc ./results/

## Rename output file:
mv multiqc_report.html supplemental_S3_multiQC_report_all_24_samples.html
