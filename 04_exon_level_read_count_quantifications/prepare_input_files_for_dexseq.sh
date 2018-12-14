## Load parallel:
module load parallel

## Link output from splice alignment folder and use as input:
ln -s ../2018-11-30-spliced_alignment_queen-tmp/tmp/quant/ input

## Link folder containing temporary files required for running DEXSeq:
ln -s ../2018-11-30_spliced_alignment_queen-tmp/tmp/ .

## Create soft link for dexseq python scripts:
ln -s ../../soft/DEXSeq_python1.26.0 dexseq_py

## Create soft link to virtual environment for running python scripts:
ln -s ../2018-11-30-spliced_alignment_worker/htseq .

## Activate virtual environment:
source ./htseq/bin/activate

## Using DEXSeq, extract annotations from a user-provided annotation file (.gtf):
python ./dexseq_py/dexseq_prepare_annotation.py tmp/Bter1.0.40.gtf tmp/Bter1.0.40.gtf.gff

## Count reads aligned over exons and output results in a temporary folder:
ls tmp/quant/*.sam | parallel -v python dexseq_py/dexseq_count.py -s reverse tmp/Bter1.0.40.gtf.gff {} {}.counts

## Create output results directory
mkdir results

## Move output count files into new results directory:
mv tmp/quant/*counts results
cd results

## Link DEXSeq prepared annotation files:
ln -s ../tmp/Bter1.0.40.gtf.gff .
