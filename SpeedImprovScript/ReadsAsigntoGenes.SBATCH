#!/bin/bash
#SBATCH --job-name=extract-read-counted --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=01:00:00
#SBATCH --mem=5G

#-----------------------------------------------------------------------------#
# This script runs get info from HTseq output.out files #
#-----------------------------------------------------------------------------#

#- Set variables --------------------------------------------------------------$
# Set the directory containing your files

DATA_DIR=/home/sl1453/2.17counting

OUTPUT_FILE="${DATA_DIR}/extracted_data.txt"
echo "File,Reads_Assigned_to_Genes" > "$OUTPUT_FILE"

# Loop through each file and extract data
for file in ${DATA_DIR}/*_htseqCount; do
    ReadsAsigntoGenes=$(head -n -5 "$file" | awk '{total += $2} END {print total}')
    filename=$(basename "$file")
    
    
    # Write the extracted data to the output file
    echo "${filename},${ReadsAsigntoGenes}" >> "$OUTPUT_FILE"

done
