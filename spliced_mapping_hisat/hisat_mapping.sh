module load gcc/7.1.0 

ln -s ../../soft/hisat2-2.1.0 hisat

mkdir ~/scratch/hisat_queen
ln -s ~/scratch/hisat_queen tmp

# obtain exon locations first in 0-based BED-like format
gunzip --stdout input/reference/Bter1.0.40.gtf.gz > tmp/Bter1.0.40.gtf

hisat/hisat2_extract_splice_sites.py tmp/Bter1.0.40.gtf > tmp/Bter_1.0.splice_sites
hisat/hisat2_extract_exons.py        tmp/Bter1.0.40.gtf > tmp/Bter_1.0.exons

# build index with these
gunzip --stdout input/reference/Bter_1.0.fa.gz > tmp/Bter_1.0.fa
./hisat/hisat2-build -p 36 --ss tmp/Bter_1.0.splice_sites --exon tmp/Bter_1.0.exons tmp/Bter_1.0.fa tmp/Bter_1.0 > tmp/Bter_1.0.log  2> tmp/Bter_1.0.err



# map each sample
treatments=(CON CLO IMI)

for treatment in $treatments; do
  echo "running with ${treatment}"
  fastqs=`ls input/fastq_all/2016-Bter-${treatment}[-_]*-Q-*gz`
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

sh hisat_commands.sh > hisat_commands.sh.log  2> hisat_commands.sh.err

## Extract annotations:
python ../../../src/DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py Bombus_terrestris.Bter_1.0.39.gtf Bombus_terrestris.Bter_1.0.39.gff3

## Count reads aligned over exons:
for name in *.sam; 
do 
    python ../../../src/DEXSeq/inst/python_scripts/dexseq_count.py Bombus_terrestris.Bter_1.0.39.gff3 "$name" "$name".fb.txt; 
done
