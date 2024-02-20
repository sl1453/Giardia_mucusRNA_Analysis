#!/bin/bash
#SBATCH --job-name=read-extractor
#SBATCH --output=read_extractor_%j.out

#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

# Change to the directory containing your files
#---------------------------------------------------------------------
cd ~/raw_untrimmed

output_file="read_numbers_summary.txt"
echo "File Name,Read Number" > "$output_file"

# Loop through all _stats.txt files in the directory
for file in *_stats.txt; do
    # Extract the read number
    read_number=$(grep "reads" "$file" | awk '{print $4}')

    # Output the file name and the read number to the output file
    echo "${file},${read_number}" >> "$output_file"


done

#----------------------------------------------------------------------#
# Display the contents of the output file                             #
#---------------------------------------------------------------------#
echo "Summary of Read Numbers:"
cat "$output_file"
