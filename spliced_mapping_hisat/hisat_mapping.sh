# from /data/home/btw928/autoScratch/2017-08-11_bombus_transcriptomics_pesticides_exp_1/results/2018-02-24_multi_analysis/2018-05-14_dexseq_analysis/

## Create hisat2 index:
../../../src/hisat2-2.1.0/hisat2-build ./GCF_000214255.1_Bter_1.0_genomic.fna Bter_1.0

## Create alignment files against the ensembl genome:
for name in linked_worker_data/*fastq; 
do 
new_name="$(echo "$name" | cut -d '/' -f 2 - | cut -d '.' -f 1 - )";  
../../../src/hisat2-2.1.0/hisat2 -x Bter_1.0 -U "$name" -s > "$new_name".sam; 
done

## Extract annotations:
python ../../../src/DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py Bombus_terrestris.Bter_1.0.39.gtf Bombus_terrestris.Bter_1.0.39.gff3

## Count reads aligned over exons:
for name in *.sam; 
do 
python ../../../src/DEXSeq/inst/python_scripts/dexseq_count.py Bombus_terrestris.Bter_1.0.39.gff3 "$name" "$name".fb.txt; 
done
