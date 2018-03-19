#### Pool seq samples ####

#pool hiseq with next seq data

cat 2016-Bter-CON-C61-4-W1-head-54247263.combined.fastq 2016-Bter-CON-C61-4-W1-head-54247263_combined_hiseq.fq > 2016-Bter-CON-C61-4-W1-head-54247263_combined_hiseq_nextseq.fq

#Run fastqc on fastq.gz files
for gz in .*.gz
do fastqc "(2016-Bter-CLO-C02-4-W1-head-54242293.combined.fastq.gz)" "$gz"
done

#create link to genomic index 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_cds_from_genomic.fna.gz 

gunzip GCF_000214255.1_Bter_1.0_cds_from_genomic.fna.gz

ln -s ../../ 

ln -s ../../data/GCF_000214255.1_Bter_1.0_cds_from_genomic.fna . 


######## RUNNING KALLISTO ##########
#create an index for kallisto
kallisto index -i Bter_1.0 GCF_000214255.1_Bter_1.0_rna_from_genomic.fna
#creates index named 'Bter_1.0'

#run quantification e.g.
kallisto quant -i Bter_1.0 -o 2016-Bter-CLO_quant_v2 --single -l 300 -s 20 2016-Bter-CON-C61-4-W1-head-54247263_combined_hiseq_nextseq.fq -t 8

#loop for all files
for name in *.fq
do
new_name="$(echo "$name" | cut -d '.' -f 1,3 - )"
echo "$new_name"
done #cuts filenames

#kallisto quant and rename files
for name in *.fq; do new_name="$(echo "$name" | cut -d '.' -f 1,3 - )";
> kallisto quant -i Bter_1.0 -o "$new_name" --single -l 300 -s 20 "$name" -t 8
> done


####### run kallisto with pseudobam option ###### 
#also convert .sam to .bam
for name in *.fq; 
do new_name="$(echo "$name" | cut -d '.' -f 1,3 - )_pseudo";
kallisto quant -i Bter_1.0 -o "$new_name" --pseudobam --single -l 300 -s 20 "$name" | samtools view -Sb - > "$name.bam" ; done


#check bam files with:
samtools view | grep -c '@NS'
grep -c '@K00'

#create stat file from bam file
#e.g. samtools flagstat 2016-Bter-CON_C06-4-W1-head-54232438_combined_hiseq_nextseq.bam > 2016-Bter-CON_C06-4-W1-head-54232438_combined_hiseq_nextseq.stat

for name in *.fq.bam ; do new_name="$(echo "$name" | cut -d '.' -f 1,1 - ).stat"; echo "$new_name"; 
> samtools flagstat "$name" > "$new_name"
> done



#run qualimap on .bam files
#run multiqc on qualimap files 
