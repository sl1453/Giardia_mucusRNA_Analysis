#!/bin/bash
#SBATCH --job-name=STAR_index --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=/home/sm3679/culex_biting/culex_genome
refgen_index=/home/sm3679/culex_biting/culex_genome/index_genome


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \
        --sjdbGTFfile ${refgen_dir}/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gtf \
        --sjdbOverhang 135 \
        --genomeSAindexNbases 13

#- FIN -----------------------------------------------------------------------#
