#!/bin/bash
#SBATCH --job-name=convert_ai --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script uses eutils to get locus tags (other aliases) for current gene ids #
#-----------------------------------------------------------------------------#

geneIDs=~/culex_biting/culex_genome/geneIDlist_segmentai


for line in `cat ${geneIDs}`
do 
efetch -db gene -id ${line} | grep -e "1. LOC" -e "Other Aliases" | paste - - | sed 's/1. //g' | sed 's/Other Aliases: //g' >> output_ai.txt
done
