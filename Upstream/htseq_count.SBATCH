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

bam_dir=/home/sl1453/samSortBam
count_dir=/home/sl1453/counting
htseq=/home/sl1453/.conda/envs/htseq-clip/bin/htseq-count
ref=/home/sl1453/genome/GiardiaDB-66_GintestinalisAssemblageAWB2019.gff

#- RUN fastqc ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.sorted.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.sorted.bam`
${htseq} -f bam -r pos -s reverse -t exon -i gene_id --nonunique fraction ${bam_dir}/${base}_Aligned.out.sorted.bam ${ref} > ${count_dir}/${base}_htseqCount

done
