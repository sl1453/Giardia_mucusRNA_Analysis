#!/bin/bash
#SBATCH --job-name=fastqc --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=1G

#-----------------------------------------------------------------------------#
# This script gives quality control of fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

module load fastqc/0.11.9
fastqc -o /home/sl1453/Qality/fastqc_Luke_trimmed /home/sl1453/raw_Luke_trimmed/*fastq.gz

#- FIN -----------------------------------------------------------------------#
