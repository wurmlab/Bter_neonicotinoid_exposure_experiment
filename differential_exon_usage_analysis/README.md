## Differential exon usage analysis:  

Steps involved in the preparation of raw data for differential exon usage analysis.  
Alignment of data was performed using HISAT2("https://ccb.jhu.edu/software/hisat2/index.shtml") as recommended by Mike Love ("https://support.bioconductor.org/p/108576/").

1. Create hisat2 index:
```
hisat2-build ./GCF_000214255.1_Bter_1.0_genomic.fna Bter_1.0
```

2. Create alignment files against the ensembl genome:
```
hisat2 -x Bter_1.0 -U sample.fastq -s > sample.sam
```

3. Create non-redundant exon bins for counting the number of reads mapped:  
```
python dexseq_prepare_annotation.py Bombus_terrestris.Bter_1.0.39.gtf Bombus_terrestris.Bter_1.0.39.dexseq.gff
```

4. Extract counts for number of reads aligned over exons:
```
python dexseq_count.py Bombus_terrestris.Bter_1.0.39.gff3 sample.sam sample.fb.txt 
```

5. Perform DEU analysis using ```dexseq_analysis.Rmd```. 

6. Prepare DEXSeq output for input into ```go_enrichment_analysis.Rmd```.  

DEXSeq outputs a tab-delimited text file with multiple fields, including p-values for each exon per gene within the bumblebee genome.    
For running ```go_enrichment_analysis.Rmd```, we only require two fields, which contain information on gene name and adjusted p-values.  

We use the script ```get_lowest_exon_p_value.sh```, which takes the output of ```dexseq_analysis.Rmd```, extracts fields of interest and formats into the input format for running ```go_enrichment_analysis.Rmd```.  

To run ```get_lowest_exon_p_value.sh```, use:  
```
./get_lowest_exon_p_value.sh results/dexseq_results_all.txt input/dexseq_results_all.input_to_topgo.txt
```
