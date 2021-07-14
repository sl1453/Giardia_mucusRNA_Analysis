# Obtaining CpipJ IDs from past assembly for gene IDs of the newest assembly
Because we are working with the newest *Culex quinquefasciatus* assembly from NCBI (GCF_015732765.1), we ran into issues where downstream analysis for GO and KEGG pathway enrichment were hindered since those databases had annotations based on past *Culex quinquefasciatus* assemblies using locus tags (i.e. the CpipJ IDs as identifiers). In order to use these databases we decided to obtain the locus tag aliases for each of the gene IDs that we could. This was done using Entrez E-utilities to pull down alias information for each gene ID from NCBI.

## Make a list of gene IDs based on the gtf file
```
grep -E 'gene_id' GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gtf | cut -d '"' -f 2 | uniq > geneIDlist.txt
```
I decided to split the file into multiple smaller files and ran the script on each of these
```
split -l 2100 geneIDlist.txt geneIDlist_segment
```

## Using efetch to get the alias information 
I ran a script on each of the multiple gene ID files, using efetch to pull alias information from NCBI. The older assembly that NCBI has alias information for is CulPip1.0 (GCA_000209185.1)


[script](https://github.com/srmarzec/Culex_Biting_RNAseq/blob/main/misc/convert_geneID_locusTag.sh)

Because this script queries online, it sometimes fails and retries for some IDs (perhaps timed out while trying for the URL). I ended up pulling the logs and finding a list of the IDs that fail and rerunning those IDs. I did this for each output-log, and then for the resulting output-log to make sure all IDs had an alias if avaiable.

## Combine all the output files from various scripts.
In the end we should have one long list with current gene IDs in column 1 and locus tags/other aliases in column 2.
```
cat output_*.txt >> newID_oldAlias.txt
```
This combined file has 24252 lines which seems right

Next, I want to keep only the aliases in the second column which are old locus tags (i.e. must start with CpipJ_CPIJ). Some loci have multiple aliases (see below)
```
LOC6032085      CpipJ_CPIJ001082, ADA
LOC6032086      CpipJ_CPIJ001083
LOC6032087      CpipJ_CPIJ001084
LOC6032089      CpipJ_CPIJ001086
LOC6032091      CpipJ_CPIJ001088, mahogunin
LOC6032092      CpipJ_CPIJ001089
LOC6032095      CpipJ_CPIJ001092, bitesize
```

I will print the first column as it is and then print everything before the comma (the CpipJ locus tag) in the second column
```
awk '{sub(/,.*$/,"",$2); print $1"\t"$2}' newID_oldAlias.txt > newID_oldLocusTag.txt
```
Get only the loci that have a CpipJ locus tag.
```
grep -E 'CpipJ_CPIJ' newID_oldLocusTag.txt > newID_oldLocusTag_ONLY.txt
```

