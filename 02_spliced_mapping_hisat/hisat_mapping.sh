## Load required modules:
module load gcc/7.1.0 

## Create soft links for required software:
ln -s ../../soft/hisat2-2.1.0/hisat2-build
ln -s ../../soft/hisat2-2.1.0/hisat2

## Make temporary directory on scratch and link to working directory:
mkdir ~/scratch/
ln -s ~/scratch/hisat_w tmp

## For the purpose of performing splice-aware aligments, exon locations are required to be extracted
## from a user-provided annotation file (gtf file).  
## obtain exon locations first in 0-based BED-like format
gunzip --stdout input/reference/Bter1.0.40.gtf.gz > tmp/Bter1.0.40.gtf

## Using HISAT2, extract splice sites and exons from the user-defined annotation file:
hisat/hisat2_extract_splice_sites.py tmp/Bter1.0.40.gtf > tmp/Bter_1.0.splice_sites
hisat/hisat2_extract_exons.py        tmp/Bter1.0.40.gtf > tmp/Bter_1.0.exons

## Usun these files, build a reference for alignments:  
gunzip --stdout input/reference/Bter_1.0.fa.gz > tmp/Bter_1.0.fa
./hisat/hisat2-build -p 36 --ss tmp/Bter_1.0.splice_sites --exon tmp/Bter_1.0.exons tmp/Bter_1.0.fa tmp/Bter_1.0 > tmp/Bter_1.0.log  2> tmp/Bter_1.0.err

## Similar to kallisto, define treatments:
treatments=(CON CLO IMI)

## For each sample in each treatment:
## First, extract colony information present within the input name to use as an abbreviated name for temporay and output files. 
## Second, print commands for performing splice-aware alignment for each sample into a shell script.  
for treatment in $treatments; do
  echo "running with ${treatment}"
  fastqs=`ls input/fastq_all/2016-Bter-${treatment}[-_]*-W[0-9]-*gz`
  colonies=$(echo $fastqs | cut -d '-' -f 4 | sort -u)
  colonies=($(echo "$colonies"))

  mkdir -p tmp/quant
  for colony in $colonies; do
    echo "running with ${colony}"
    fastqs_for_colony=`echo $fastqs | grep "\-${colony}\-" | tr '\n' ' ' | sed  -e 's/ i/,i/g'`
    output=tmp/quant/${treatment}-${colony}
    echo "./hisat/hisat2 -x tmp/Bter_1.0 -U ${fastqs_for_colony} --known-splicesite-infile tmp/Bter_1.0.splicesites.txt --threads=35 -S ${output}.sam" >> ./hisat_commands.sh
  done
done

## Run shell script to perform splice-aware alignments for each sample: 
sh hisat_commands.sh > tmp/hisat_commands.sh.log  2> tmp/hisat_commands.sh.err
