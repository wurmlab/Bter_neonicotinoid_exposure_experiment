#Filenames are like 
#2016-Bter-CTH-C60-4-W1-head.D711_D506.HM3WHBBXX.s_6.R1.fastq.gz
#2016-Bter-CYH-C03-4-Q-head.D711_D503.HM3WHBBXX.s_6.R1.fastq.gz 

treatments=(CON CLO IMI)


###
mkdir ~/scratch/2018-09-19-kallisto_actualworkers
ln -s ~/scratch/2018-09-19-kallisto_actualworkers tmp

# indexing
kallisto_ref=tmp/Bter1_cdna
./kallisto index -i ${kallisto_ref} input/reference/Bter1_cdna.fa 



## running for each
for treatment in $treatments; do
  echo "running with ${treatment}"
  fastqs=`ls input/fastq_all/2016-Bter-${treatment}[-_]*-W[0-9]-*gz`
  colonies=$(echo $fastqs | cut -d '-' -f 4 | sort -u)
  colonies=($(echo "$colonies"))

  mkdir -p tmp/quant
  for colony in $colonies; do
    echo "running with ${colony}"
    fastqs_for_colony=`echo $fastqs | grep "\-${colony}\-" | tr "\n" " "`
    output=tmp/quant/${treatment}-${colony}
    echo "./kallisto quant -i ${kallisto_ref}  --output-dir=${output} --single --fragment-length=300 --sd=20  --threads=35 ${fastqs_for_colony}" >> ./kallisto_commands.sh
  done
done

sh kallisto_commands.sh > kallisto_commands.log  2> kallisto_commands.err

# for W, some issue with colony name for CON (manually renamed), and CYH-C43 has 1 fewer file than others (asked Joe)



# Note to self. to create an array the folowing is preferred;
#    fastqs=(input/fastq_all/2016-Bter-${treatment}*-W[1-4]-*gz)   but not using here bc i need ot use them after

 mkdir results
 mv tmp/quant results/kallisto_quantifications

