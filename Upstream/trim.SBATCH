#!/bin/bash
#SBATCH --job-name=trim --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs trimmomatic to clean up raw reads #
#-----------------------------------------------------------------------------#


#- Set variables --------------------------------------------------------------$

raw_dir=/home/sl1453/raw_untrimmed

trim_dir=/home/sl1453/trimmomatic/trim_1

trim=/home/sl1453/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar

adapter=/home/sl1453/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10

#- RUN Trimmomatic-------------------------------------------------------------$

files=(${raw_dir}/*_R1.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1.fastq.gz`
java -Xmx2G -jar ${trim} PE \
          ${raw_dir}/${base}_R1.fastq.gz \
          ${raw_dir}/${base}_R2.fastq.gz \
          ${trim_dir}/${base}_1_PE.fastq.gz ${trim_dir}/${base}_1_SE.fastq.gz \
          ${trim_dir}/${base}_2_PE.fastq.gz ${trim_dir}/${base}_2_SE.fastq.gz \
          ILLUMINACLIP:${adapter} \
          HEADCROP:15 \
          TRAILING:6 \
          SLIDINGWINDOW:4:15 \
          MINLEN:50
done
