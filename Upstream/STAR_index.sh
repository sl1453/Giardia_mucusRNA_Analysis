#!/bin/bash
#SBATCH --job-name=STAR_index --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=5G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=/home/sl1453/genome
refgen_index=/home/sl1453/genome/indexed_genome


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/GiardiaDB-66_GintestinalisAssemblageAWB2019_Genome.fasta \
        --sjdbGTFfile ${refgen_dir}/GiardiaDB-66_GintestinalisAssemblageAWB2019.gff \
        --genomeSAindexNbases 11

#- FIN -----------------------------------------------------------------------#
