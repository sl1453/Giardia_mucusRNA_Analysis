#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

source activate htseq-clip 

#- Set variables ----------------------------------------------------------------#

bam_dir=/home/mch284/Auto_biting/mRNA/mapped_reads/two_pass/bam_dir
count_dir=/home/mch284/Auto_biting/mRNA/fr_counts_dir
htseq=/home/mch284/.conda/envs/htseq-clip/bin/htseq-count
ref=/home/mch284/genomes/AalbF3/aedes_albopictus_AalbF3.gff3

#- RUN fastqc ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s reverse -t exon -i gene_id ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount

done
