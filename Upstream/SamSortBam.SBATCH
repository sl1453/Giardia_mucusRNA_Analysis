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

sam_dir=/home/sl1453/StarMap/sam_dir
bam_dir=/home/sl1453/samSortBam

#- RUN fastqc ----------------------------------------------------------------#

files=(${sam_dir}/*.out.bam)
for file in ${files[@]}
do
base=`basename ${file} .out.bam`
samtools sort -o ${bam_dir}/${base}.out.sorted.bam
samtools index ${bam_dir}/${base}.out.sorted.bam

done

#- FIN -----------------------------------------------------------------------#
