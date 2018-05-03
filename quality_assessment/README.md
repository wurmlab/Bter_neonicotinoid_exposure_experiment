## Quality assessment of RNASeq sequences generated from pesticide-exposure studies.  

For quality assessment, hiseq and nextseq runs were combined per sample and aligned against the reference
buff-tailed bumblebee _Bombus terrestris_ reference genome obtained from Ensembl (release .   
Alignment was performed against STAR using the reference GTF available from Ensembl (release 39).  

Reference genome was first indexed using STAR with the following script:  
```
./make_STAR_index.sh input.fasta input.gtf
```
