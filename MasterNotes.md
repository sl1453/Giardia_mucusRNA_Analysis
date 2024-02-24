# Workflow/pipeline details for the Giardia response to viscoelasticity RNA-seq analysis
### Automation & Software Install (Virtue Env Setup)
Along the way, you might need some automations to quickly extract read numbers or other info from various output files in batches, zipped or unzipped.
Dig in ([here](https://github.com/sl1453/Giardia_mucusRNA_Analysis/tree/main/SpeedImprovScript)) to find scripts that might suit the need for speeding.

If you need more updated softwares than what I used, ([look here](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/SoftwareInstruct/How2SetUp.md))for a beginners instruction on how to create your own environemnt and incorporate updated softwares to sbatch scripts.

*If you are a command line pro, jump to the next section.
## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found ([here](https://docs.google.com/spreadsheets/d/17yXVA9PE-rkG_24pJ2p31cf8QelNo2Ly/edit?rtpof=true#gid=566235505))

### Data Accession
Data was generated, raw reads...

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Raw Reads came in Trimmed and Untrimmed from the sequencing center. Parameter used by Maryland Sequencing Center are: 
Trimmomatic (version 0.33)
parameters: simple clip threshold=7, seed mismatches=2, palindrome threshold=40, minimum sequence length=30, and trailing quality=20.
(TruSeq3-PE-2.fa:2:40:7)

FastQC (v0.11.9) was used for quality control visualization     
Raw_untrimmed:([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/fastqc_untrimmed.SBATCH))    
Sequencing_center_trimmed: ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/fastqc_umd_trimmed.SBATCH))

Most samples have warnings for "per sequence GC content", which I expected because Giardia has high GC contents; and poor "sequence duplication levels" which I expect for RNAseq data. The "per base quality" score always dropps at the end of fragments, and ONLY seen in read 2 sequences.

The QC for the prileminary trimming from Luke has successfully removed the adaptor warnings, however, still show poor "per base GC contents" for the first ~15 bases. I decided to use HEADCROP flag to remove the first 15 bases for the UNTRIMMED RAW READS while keeping the adaptor trimming parameters similar to Luke's. All other flags (TRAILING, SLIDINGWINDOW, and MINLEN) are rather general/default for basic quality of bases and did not result in much difference of trimming. 

Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/trim.SBATCH))

My trimming results in shorter sequences ~130bp, improved the "per base GC contents". The end sequece "per base quality" score dropping was not improved as I only kept the "Trailing" parameter as low as 6, to advoid extremmely short fragemnts. 


### Mapping

Mapping was done using the *_Giardia_ assenblage A WB * reference genome from Giardia DB website ([xu et al. 2019 version](https://giardiadb.org/giardiadb/app/record/dataset/TMPTX_gassAWB2019)). ([GFF file](https://giardiadb.org/giardiadb/app/downloads/Current_Release/GintestinalisAssemblageAWB2019/gff/data/))

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/STAR_index.SBATCH)).

Reads were mapped ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/STAR_map.SBATCH)). Becasue most of VSP genes are duplicated, parameter "outFilterMultimapNmax 2" was chosen;  "alignIntronMax 1" was to account for the fact that _Giardia_ has very few introns (8 cis, 5 trans); "outSAMtype BAM Unsorted". Due to the strand-specific library we used, it is good to add the strand parameter in STAR for stranded alignment "--outSAMstrandField intronMotif".

Output bam files were sorted, and then indexed ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/SamSortBam.SBATCH))

### Generating count matrix with HTSeq (htseq-count)

HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/htseq_count.SBATCH))

Parameters: 
-s reverse                :ths is for "first-stranded" library, using NEB Ultra II Directional kit.
--nonunique fraction      :multimapping reads are fractionally counted (each alignment carrying 1/2 count, since the parameter outFilterMultimapNmax was previously set to 2 in the alignment step). 
or
--nonunique none (default): the read (or read pair) is counted as ambiguous and not counted for any features.

### Alternative Generating count matrix with featurecounts


## Downstream

All downstream analysis done in R using R version 4.0.2

### DESeq
Using DESeq2 (v1.30.1) ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Downstream/DESeq.R))

### GO Enrichment
Using topGO (v2.42.0) ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Downstream/GOenrichment.R))

### KEGG Pathway Enrichment
Using a custom script for KEGG term enrichment in pathways all pulled from KEGG ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Downstream/KEGG.R))

### Obtaining CpipJ IDs from past assembly for gene IDs of the newest assembly
Because we are working with the newest Culex quinquefasciatus assembly from NCBI (GCF_015732765.1), we ran into issues where downstream analysis for GO and KEGG pathway enrichment were hindered since those databases had annotations based on past Culex quinquefasciatus assemblies using locus tags (i.e. the CpipJ IDs as identifiers). Details on how we obtained locus tags is [here](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/misc/GeneID_LocusTag_Conversion.md).
