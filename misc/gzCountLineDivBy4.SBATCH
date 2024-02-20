#!/bin/bash
#SBATCH --job-name=line-count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=01:00:00
#SBATCH --mem=1G

#-----------------------------------------------------------------------------#
# This script runs counting reads for FASTQ.gz files #
#-----------------------------------------------------------------------------#

#- Set variables --------------------------------------------------------------$
  
raw_dir=/home/sl1453/trimmomatic/trim_1

counts=/home/sl1453/trimmomatic/trim_1/read_counts

#- RUN LOOP-------------------------------------------------------------$
# Loop through each gz file and calculate the line count divided by 4

files=(${raw_dir}/*PE.fastq.gz)
for file in ${files[@]}
do
    echo "Processing $file"
    count=$(zcat "$file" | wc -l | awk '{print $1/4}')
    echo "$file: $count" >> ${counts}/read_counts.txt

  
done
