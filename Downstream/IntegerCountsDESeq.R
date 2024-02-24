# This script runs DESeq (specifically starting with htseq count files)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install('EnhancedVolcano')
# BiocManager::install("apeglm")
# BiocManager::install("GOplot")

#Load Libraries
library(DESeq2); library(ggplot2); library(EnhancedVolcano); library(GOplot); library(mygene); library(tidyverse); library(pheatmap); library(reshape2); library(viridis)

#library(reshape2); library(viridis); 

# Set working directory to source file location

#Choose directory with htseq-count data
directory<-("../data")

#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- list.files(directory)
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- substr(sampleFiles, 1, 1)#to get conditions I'm pulling the first letter, which is either M (molestus, non-biting) or P (pipiens, biting)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
#dds$condition <- relevel(dds$condition, ref = "ND")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

# Transform normalized counts using the rlog function ("RNASEQ20_Day3_HandsOn.pdf")
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")

# Hierarchial Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld) #retrieve matrix from the rld object
# compute pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# plot the correlation c=values as a heatmap
pheatmap(rld_cor)

# # Look at normalized data
dds_counts <- estimateSizeFactors(dds)
dds_counts <- counts(dds_counts, normalized = TRUE)

# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="condition_P_vs_M", type="apeglm")
res_LFC

# Order the table by smallest p value
resOrdered <- res_LFC[order(res_LFC$pvalue),]
summary(resOrdered)


#Basic MA-plot
plotMA(res, ylim=c(-3,3))
plotMA(res_LFC, ylim=c(-3,3)) # See lots of shrinkage towards the beginning

#If you want to interactivly find genes associated with certain points
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
# Note for selecting fold change cutoff: log2foldchange 0.58 is equal to a 1.5 fold change
EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-2.5, 2.5),
                ylim = c(0,25),
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

# Write out a table of these significant differentially expressed genes
write.csv(dplyr::select(sig_res, gene, log2FoldChange, padj), 
            file="../misc/PvM_LFCshrink_padj.txt", row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res %>% dplyr::select(gene), 
            file="../misc/PvM_DEGids.txt", col.names = F, row.names = F, quote = F)

##################
# Trialing making a heatmap based on the normalized counts

## Extract normalized expression for significant genes, and set the gene column (1) to row names
norm_DEGsig <- dds_counts %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% sig_res$gene) %>% 
  column_to_rownames(var = "gene") 

## Annotate our heatmap (optional)
annotation <- data.frame(row.names = sampleNames,
                          type = sampleConditions)

## Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

## Run pheatmap
pheatmap(norm_DEGsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

##################
# Make output excel sheet to match M.F.P.'s old excel sheets
res_print <- as.data.frame(res_LFC)
res_print$geneID <- row.names(res_print)

head(res_print)

# Find gene names with mygene package
dat <- queryMany(res_print$geneID, scopes="symbol", fields=c("name"))

res_merged <- merge(res_print, dat, by.x="geneID", by.y="query")

keep_cols <- c("geneID", "name", "log2FoldChange", "pvalue", "padj")


# Write out a csv with these data
write.csv(res_merged[keep_cols], 
          file="../misc/Culex_DESeq_results.csv", row.names = F)

