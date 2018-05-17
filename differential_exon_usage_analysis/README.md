## Differential exon usage analysis:  

Steps involved in the preparation of differential exon usage analysis.  
Steps involved include:
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

5. Perform DEU analysis using ``` ```. 
