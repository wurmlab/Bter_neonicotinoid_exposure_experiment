#!/bin/env bash

## Create STAR indices
../../../../../bin/STAR --runThreadN 10 \
               --runMode genomeGenerate \
               --genomeDir . \
               --genomeFastaFiles GCF_000214255.1_Bter_1.0_genomic.fna \
               --sjdbGTFfile GCF_000214255.1_Bter_1.0_genomic.gtf \
               --sjdbOverhang 74
