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

# Create an output file to store the results
output_file="read_numbers_only.txt"

# Loop through all _stats.txt files in the directory
for file in *_stats.txt; do
    # Extract the read number
    read_number=$(awk '/reads/{print $(NF-1)}' "$file")

    # Output only the read number to the output file
    echo "${read_number}" >> "$output_file"
done

#----------------------------------------------------------------------#
# Display the contents of the output file                             #
#---------------------------------------------------------------------#
echo "Summary of Read Numbers:"
cat "$output_file"

