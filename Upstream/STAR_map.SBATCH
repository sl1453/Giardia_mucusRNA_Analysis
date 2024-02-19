#!/bin/bash
#SBATCH --job-name=STAR_map --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=14:00:00
#SBATCH --mem=10G

#-----------------------------------------------------------------------------#
# This script maps reads to the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

trim_dir=/home/sl1453/trimmomatic/trim_1

out_dir=/home/sl1453/StarMap/sam_dir

refgen_dir=/home/sl1453/genome/indexed_genome


#- RUN STAR----------------------------------------------------------------#

files=(${trim_dir}/*_1_PE.fastq.gz)
for file in ${files[@]}
do
base=`basename ${file} _L004_1_PE.fastq.gz`
STAR --genomeDir ${refgen_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${trim_dir}/${base}_L004_1_PE.fastq.gz ${trim_dir}/${base}_L004_2_PE.fastq.gz \
        --outFilterMultimapNmax 2 \
        --alignIntronMax 1 \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ${out_dir}/${base}_  
done

#- FIN -----------------------------------------------------------------------#
