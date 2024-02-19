#!/bin/bash
#SBATCH --job-name=fastqc --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=5G

#-----------------------------------------------------------------------------#
# This script gives quality control of fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

module load fastqc/0.11.9
fastqc -o /home/sl1453/Qality/fastqc_trim_1 /home/sl1453/trimmomatic/trim_1/*fastq.gz

#- FIN -----------------------------------------------------------------------#
