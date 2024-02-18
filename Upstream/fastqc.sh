#!/bin/bash
#SBATCH --job-name=fastqc --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sl1453@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=1G

#-----------------------------------------------------------------------------#
# This script gives quality control of fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

/home/sm3679/bin/FastQC/fastqc -o /home/sm3679/culex_biting/fastqc_dir /home/sm3679/culex_biting/trim_dir/*PE.fastq.gz

#- FIN -----------------------------------------------------------------------#
