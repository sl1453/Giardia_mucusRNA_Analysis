#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=50G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

#- Load the module for Anaconda-----------------------------------------------#
module load anaconda3


#- Activate my conda virtue environment called htseq-env----------------------#
conda activate htseq-env


#- Set variables ----------------------------------------------------------------#

bam_dir=/home/sl1453/samsortbam
count_dir=/home/sl1453/htseq2.17Counting
ref=/home/sl1453/genome/GiardiaDB-66_GintestinalisAssemblageAWB2019.gff


#- RUN htseq-count ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.sorted.bam)
for file in ${files[@]}
do
    base=$(basename ${file} _Aligned.out.sorted.bam)
    htseq-count -f bam -r pos -s reverse -t exon -i gene_id --nonunique fraction ${bam_dir}/${base}_Aligned.out.sorted.bam ${ref} > ${count_dir}/${base}_htseqCount
done


#- Deactivate conda environment-------------------------------------------------#
conda deactivate

