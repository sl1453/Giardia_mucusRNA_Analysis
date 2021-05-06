#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

source /home/sm3679/python-env/bin/activate

#- Set variables ----------------------------------------------------------------#

bam_dir=/home/sm3679/culex_biting/bam_dir
count_dir=/home/sm3679/culex_biting/counts_dir
htseq=/home/sm3679/python-env/bin/htseq-count
ref=/home/sm3679/culex_biting/culex_genome/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gtf

#- RUN fastqc ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s yes -t exon -i gene_id ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount

done

#- FIN -----------------------------------------------------------------------#
