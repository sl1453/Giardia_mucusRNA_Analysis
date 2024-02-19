# Workflow/pipeline details for the Giardia response to viscoelasticity RNA-seq analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found [here](https://docs.google.com/spreadsheets/d/1y15EVJ7VUNeKWtLNaBUMJ1zZaR_LLv7YgWeYWtGrIpI/edit?usp=sharing)

### Data Accession
Data was generated, raw reads...

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Raw Reads came in Trimmed and Untrimmed from the sequencing center. Parameter used by Maryland Sequencing Center are: 
Trimmomatic (version 0.33)
parameters: simple clip threshold=7, seed mismatches=2, palindrome threshold=40, minimum sequence length=30, and training quality=20.

FastQC (v0.11.9) was used for quality control visualization     
Raw_untrimmed:([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/fastqc.SBATCH))    
Sequencing_center_trimmed: ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/fastqc_umd_trimmed.SBATCH))

Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/sl1453/Giardia_mucusRNA_Analysis/blob/main/Upstream/trim.sh))


Preliminary trimming and fastqc showed a poor "per sequence base content" for the first ~15 bases. We decided to use HEADCROP flag to remove the first 15 bases. All other flags (TRAILING, SLIDINGWINDOW, and MINLEN) are rather general/default for basic quality of bases and did not result in much difference of trimming.

Most samples have warnings or fail for "per sequence GC content" and "sequence duplication levels", which I expect for RNAseq data. An interesting occurence, using the HEADCROP flag to remove the first 15 bases results in all samples failing "per sequence tile quality", while trimming without HEADCROP results only in a warning for this. This is likely because the samples are now being measured over 135b instead of 150b, which I guess makes this check fail for sequence tile quality.

"Overrepresented sequences" all blast to Culex (and related mosquitoes) or else are too repetitive to map to any organism. All the overrepresented sequences from each sample can be found [here](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/misc/OverrepSequences.txt), annotated with whether they could be blasted mapped and to what.

### Mapping

Mapping was done using the *Culex quinquefasciatus* reference genome ([GCF_015732765.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_015732765.1/)) found on NCBI

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Upstream/STAR_index.sh)).

Reads were mapped in a two pass method. The first pass followed typical method with splice junctions from annotations ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Upstream/STAR_map.sh)). The second pass is similar except that it additionally uses the output splice junctions info from the first pass (these would be novel splice junctions) to facilitate mapping ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Upstream/STAR_map_twopass.sh)).

Output sam files were converted to bam and then indexed ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Upstream/sam2bam.sh))

### Generating count matrix with HTSeq (htseq-count)

HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses ([script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/Upstream/htseq_count.sh))

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
