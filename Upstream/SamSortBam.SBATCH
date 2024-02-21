#!/bin/bash
#SBATCH --job-name=sam2bam --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script converts sam to bam, sorts/indexes the bam files #
#-----------------------------------------------------------------------------#

module load samtools/1.9

#- Set variables ----------------------------------------------------------------#

unsortedBam_dir=/home/sl1453/StarMap/unsortedBam_dir
bam_dir=/home/sl1453/samsortbam

#- RUN samSortBam -------------------------------------------------------------$

files=(${unsortedBam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
  base=$(basename $file _Aligned.out.bam)
  samtools sort -o ${bam_dir}/${base}_Aligned.out.sorted.bam $file
  samtools index ${bam_dir}/${base}_Aligned.out.sorted.bam


done

#- FIN -----------------------------------------------------------------------#
