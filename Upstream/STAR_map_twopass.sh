#!/bin/bash
#SBATCH --job-name=STAR_map_twopass --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script maps reads to the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

trim_dir=/home/sm3679/culex_biting/trim_dir

splice_junctions=/home/sm3679/culex_biting/sam_dir

out_dir=/home/sm3679/culex_biting/sam_dir/twopass

refgen_dir=/home/sm3679/culex_biting/culex_genome/index_genome


#- RUN STAR----------------------------------------------------------------#

files=(${trim_dir}/*_1_PE.fastq.gz)
for file in ${files[@]}
do
base=`basename ${file} _L005_1_PE.fastq.gz`
STAR --genomeDir ${refgen_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${trim_dir}/${base}_L005_1_PE.fastq.gz ${trim_dir}/${base}_L005_2_PE.fastq.gz \
        --sjdbFileChrStartEnd ${splice_junctions}/*_SJ.out.tab \
        --outFileNamePrefix ${out_dir}/${base}_  
done

#- FIN -----------------------------------------------------------------------#
