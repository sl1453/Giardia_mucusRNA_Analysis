# This script runs DESeq (specifically starting with htseq count files)

#Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")
BiocManager::install("GOplot")
BiocManager::install("mygene")
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
install.packages("ggrepel")
install.packages("ggplot2")

#Load Libraries
library(DESeq2); library(ggplot2); library(EnhancedVolcano); library(GOplot); library(mygene); library(tidyverse); library(pheatmap); library(reshape2); library(viridis); library(ggrepel); library(viridisLite); library(apeglm)

# Set working directory to source file location
#Choose directory with htseq-count data
setwd("/Users/apple/Documents/Physalia Course  Workshop/dAtA/DESeq2/data/htseq2.17Counting")

###############################################################################
#Round the htseq counts, DO NOT round if your count is normalized (use tools other than DESeq in that case).
htseqFiles <- list.files(pattern = "^(control|methylCell).*htseqCount$")
htseqFiles

rounded_counts_list <- lapply(htseqFiles, function(file) {
  counts <- read.table(file, header = FALSE, sep = "\t", row.names = 1)
  rounded_counts <- round(counts)
  colnames(rounded_counts) <- sub("XSLGU.*", "", file)  # I'm extracting the sample names from the long file names.
  rounded_counts
})


names(rounded_counts_list) <- htseqFiles #Assign names to each rounded data frame to its original htseq file name.

print(names(rounded_counts_list)) #Check the order of the names.

#Check the first 10 rows of the rounded list
first_file_counts <- rounded_counts_list[[htseqFiles[1]]]
head(first_file_counts, 10)

#Check if first 10 rows from the original unrounded files are same
unround_counts <- read.table(htseqFiles[1], header = FALSE, sep = "\t", row.names = 1, nrows = 10)
print(unround_counts)

#combine the rounded list into a single dataframe with samples as columns.
combined_rounded_counts <- do.call(cbind, rounded_counts_list)
colnames(combined_rounded_counts) #verify column names

##############################################################################
#Create the sample table (this could alternatively be made externally and read in)
sampleNames <- sub("XSLGU.*","",htseqFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- substr(htseqFiles, 1, 1)#to get conditions I'm pulling the first letter, which is either c, m, or x (control, methylCell, or xanthanGum).

sampleTable <- data.frame(sampleName = sampleNames,
                          condition = sampleConditions)

str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable
#Check: the column names of combined_rounded_counts should match the sample identifiers (sampleName) in sampleTable.
colnames(combined_rounded_counts)

###############################################################################
# Make the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = combined_rounded_counts,
                              colData = sampleTable,
                              design = ~ condition)
#dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  #directory = directory,
                                  #design = ~ condition)
dds
#Check all of the original htseq files contain same number of genes (5396 in this case)
rounded_counts_list <- lapply(htseqFiles, function(file) {
  counts <- read.table(file, header = FALSE, sep = "\t", row.names = 1)
  if (nrow(counts) != 5396) {
    stop("File ", file, " does not have 5396 rows.")
  }
  round(counts)
})
#Check if all original htseq files has the consistent gene_id name and order
gene_ids_list <- list()

rounded_counts_list <- list()
for (i in seq_along(htseqFiles)) {
  file <- htseqFiles[i]
  counts <- read.table(file, header = FALSE, sep = "\t", row.names = 1)
  rounded_counts <- round(counts)
  gene_ids_list[[i]] <- rownames(rounded_counts)
  rounded_counts_list[[i]] <- rounded_counts
}

all_identical <- all(sapply(gene_ids_list, function(x) all(x == gene_ids_list[[1]])))
print(all_identical)



###############################################################################
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for everything compared to the control. (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "c")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

png("library_sizes_6e7.png", width = 800, height = 600)
barplot(librarySizes, names = names(librarySizes), las = 2, ylim = c(0,6e+07), main = "Barplot of library sizes")
dev.off()

librarySizes <- colSums(counts(dds))
hist(librarySizes, main = "Histogram of Library Sizes", xlab = "Library Size", ylab = "Frequency")

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
#rld <- rlog(dds, blind = TRUE)
# manually generate PCA data to look into other PCA Levels.
sampleNames <- colnames(dds)
rld_mat <- assay(rld)  # rlog-transformed data
pca <- prcomp(t(rld_mat))
pcaData <- as.data.frame(pca$x)
pcaData$SampleName <- sampleNames
pcaData$condition <- dds$condition
pcaData$condition
dds$condition
# For example, plotting PC4 vs PC6
ggplot(pcaData, aes(x = PC1, y = PC4, color = condition)) +
  geom_point() +
  geom_text(aes(label = SampleName), vjust = 1, hjust = 1, check_overlap = TRUE)+
  xlab("PC1") +
  ylab("PC4") +
  ggtitle("PCA - PC1 vs PC4")

# Hierarchial Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld) #retrieve matrix from the rld object
# compute pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# plot the correlation c=values as a heatmap
pheatmap(rld_cor)

sampleNames <- colnames(rld_mat)
sampleNames
# Manually set conditions for each sample
sampleConditions <- c('Condition1', 'Condition2', 'Condition3', 'Condition4','methylCell_1','methylCell_2','methylCell_3','methylCell_4') 

# make name orders
rld_mat_ordered <- rld_mat[, order(sampleConditions)]
rld_cor_ordered <- cor(rld_mat_ordered)
pheatmap(rld_cor_ordered, cluster_rows = FALSE, cluster_cols = FALSE)


# # Look at normalized data
dds_counts <- estimateSizeFactors(dds)
dds_counts <- counts(dds_counts, normalized = TRUE)


# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res_m_vs_c <- results(dds)
res_m_vs_c

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
#resultsNames(dds)
#res_LFC_x_vs_c <- lfcShrink(dds, coef="condition_x_vs_c", type="apeglm")
#res_LFC_x_vs_c

# Run lfcShrink for 'm' vs. 'c' comparison
res_LFC_m_vs_c <- lfcShrink(dds, coef="condition_m_vs_c", type="apeglm")
res_LFC_m_vs_c

# Order the table by smallest p value
#resOrdered_x_vs_c <- res_LFC_x_vs_c[order(res_LFC_x_vs_c$pvalue),]
#summary(resOrdered_x_vs_c)

resOrdered_m_vs_c <- res_m_vs_c[order(res_m_vs_c$pvalue),]
summary(resOrdered_m_vs_c)
?results

summaryOutput <- capture.output(summary(resOrdered_m_vs_c))
writeLines(summaryOutput, "../../output/c_mCell/Summaryresults_m_vs_c.txt")

#Basic MA-plot
plotMA(res_x_vs_c, ylim=c(-3,3), main="MA-Plot for x vs c")
plotMA(res_LFC_x_vs_c, ylim=c(-3,3), main="Shrunk MA-Plot for x vs c")
# See lots of shrinkage towards?
plotMA(res_m_vs_c, ylim=c(-3,3), main="MA-Plot for m vs c")
plotMA(res_LFC_m_vs_c, ylim=c(-3,3), main="Shrunk apeglm MA-Plot for m vs c")


#If you want to interactivly find genes associated with certain points
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
# Note for selecting fold change cutoff: log2foldchange 0.58 is equal to a 1.5 fold change
EnhancedVolcano(res_x_vs_c,
                lab = rownames(res_x_vs_c),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-2.5, 2.5),
                ylim = c(0,10),
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0)

EnhancedVolcano(res_LFC_x_vs_c,
                lab = rownames(res_LFC_x_vs_c),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-2.5, 2.5),
                ylim = c(0,10),
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0)

EnhancedVolcano(res_m_vs_c,
                lab = rownames(res_m_vs_c),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,5),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0)

EnhancedVolcano(res_LFC_m_vs_c,
                lab = rownames(res_LFC_m_vs_c),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,5),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0)
# Can look at which of the result are significant and have a high enough log2 fold change
sig_res_x_vs_c <- res_x_vs_c %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
sig_res_x_vs_c
# Write out a table of these significant differentially expressed genes
write.csv(dplyr::select(sig_res_x_vs_c, gene, log2FoldChange, padj), 
          file="../../output/Gum_vs_C_padj.txt", row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res_x_vs_c %>% dplyr::select(gene), 
            file="../../output/Gum_vs_C_DEGids.txt", col.names = F, row.names = F, quote = F)




sig_res_m_vs_c <- res_m_vs_c %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 0.5)
sig_res_m_vs_c
# Write out a table of these significant differentially expressed genes
write.csv(dplyr::select(sig_res_m_vs_c, gene, log2FoldChange, padj), 
          file="../../output/methylCell_vs_C_padj.txt", row.names = F)
##################
# Trialing making a heatmap based on the normalized counts

## Extract normalized expression for significant genes, and set the gene column (1) to row names
norm_DEGsig_m_vs_c <- dds_counts %>% 
  data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% sig_res_m_vs_c$gene) %>% 
  column_to_rownames(var = "gene")


## Annotate our heatmap (optional)
colnames(norm_DEGsig_m_vs_c)

annotation <- data.frame(row.names = colnames(dds),
                         type = dds$condition)
rownames(annotation)
## Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

## Run pheatmap
pheatmap(norm_DEGsig_m_vs_c, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
################## heatmap All ############
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
rownames(norm_counts_df) <- rownames(dds)
annotation <- data.frame(row.names = colnames(dds), type = dds$condition)
heat_colors <- brewer.pal(6, "YlOrRd")

pheatmap(norm_counts_df, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation_col = annotation,
         border_color = NA, 
         fontsize = 10, 
         scale = "row",
         fontsize_row = 10, 
         height = 20)
##################
# Make output excel sheet to list All genes
res_print <- as.data.frame(resOrdered_m_vs_c)
res_print$geneID <- row.names(res_print)
head(res_print)

# Find gene names with mygene package
library(mygene)
dat <- queryMany(res_print$geneID, scopes="symbol", fields=c("name"))
res_merged <- merge(res_print, dat, by.x="geneID", by.y="query")

keep_cols <- c("geneID", "name", "log2FoldChange", "pvalue", "padj")


# Write out a csv with these data
write.xlsx(res_merged[keep_cols], 
           file="../../output/c_mCell/DESeq2_results_m_vs_c.xlsx", rowNames = F)

