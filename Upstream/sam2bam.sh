#!/bin/bash
#SBATCH --job-name=sam2bam --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script converts sam to bam, sorts/indexes the bam files #
#-----------------------------------------------------------------------------#

module load samtools/1.9

#- Set variables ----------------------------------------------------------------#

sam_dir=/home/sm3679/culex_biting/sam_dir
bam_dir=/home/sm3679/culex_biting/bam_dir

#- RUN fastqc ----------------------------------------------------------------#

files=(${sam_dir}/*_Aligned.out.sam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.sam`
samtools view -b -S ${sam_dir}/${base}_Aligned.out.sam | samtools sort -o ${bam_dir}/${base}_Aligned.out.bam
rm -f ${sam_dir}/${base}_Aligned.out.sam
samtools index ${bam_dir}/${base}_Aligned.out.bam

done

#- FIN -----------------------------------------------------------------------#
