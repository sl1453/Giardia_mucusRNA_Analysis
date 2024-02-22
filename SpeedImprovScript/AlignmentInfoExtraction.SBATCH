#!/bin/bash
#SBATCH --job-name=extract-data --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=01:00:00
#SBATCH --mem=5G

#-----------------------------------------------------------------------------#
# This script runs get info from STAR output.out files #
#-----------------------------------------------------------------------------#

#- Set variables --------------------------------------------------------------$
# Set the directory containing your files

DATA_DIR=/home/sl1453/StarMap/sam_dir

OUTPUT_FILE="${DATA_DIR}/extracted_data.txt"
echo "File,Number of Input Reads,Uniquely Mapped Reads Number,Uniquely Mapped Reads %" > "$OUTPUT_FILE"



# Loop through each file and extract data
for file in ${DATA_DIR}/*Log.final.out; do
    num_input_reads=$(grep "Number of input reads" "$file" | awk '{print $6}')
    uniq_mapped_reads_num=$(grep "Uniquely mapped reads number" "$file" | awk '{print $6}')
    uniq_mapped_reads_percent=$(grep "Uniquely mapped reads %" "$file" | awk '{print $6}')
    
    # Write the extracted data to the output file
    echo "${file},${num_input_reads},${uniq_mapped_reads_num},${uniq_mapped_reads_percent}" >> "$OUTPUT_FILE"

done
