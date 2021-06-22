# KEGG analysis
# The pathway enrichment analysis for this I modified from "https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html"

# BiocManager::install("KEGGREST")

# BiocManager::install("mygene")

# Load libraries

library(KEGGREST); library(tidyverse); library(kableExtra)


# Read in data
gene_list <- read.csv("../misc/PvM_LFCshrink_padj.txt", header = T)

# Read in the list of gene IDs and old locus tags that we found by querying NCBI
locus_tag_list <- read.table("../misc/newID_oldLocusTag_fixed.txt", sep = "\t", header = F)
colnames(locus_tag_list) <- c("gene","locus_tag")

# Merge the data with the locus tag table
gene_list_tag <- left_join(gene_list, locus_tag_list, by = "gene")


# If you want to write out the KEGG ID list and do this online
write.table(gene_list_tag$locus_tag, 
          file="../misc/PvM_keggID.txt", col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways
###



#Trying to automate this process with KEGG pathway enrichment
all_genes_list <- read.csv(file="../misc/Culex_DESeq_results.csv")
head(all_genes_list)

# Merge the data with the locus tag table
colnames(locus_tag_list) <- c("geneID","locus_tag")
all_genes_list_tag <- left_join(all_genes_list, locus_tag_list, by = "geneID")



# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "cqu")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

geneList <- all_genes_list_tag$padj
geneLFClist <- all_genes_list_tag$log2FoldChange
names(geneList) <- sub("cqu:","", all_genes_list_tag$locus_tag) # get rid of the beginning "cqu:" since the gene list we brought from kegg doesn't have this
names(geneLFClist) <- sub("cqu:","", all_genes_list_tag$locus_tag)
head(geneList)
head(geneLFClist)

#genes.by.pathway_40d <- genes.by.pathway[-c(99, 120)]

pathway_pval <- data.frame()

for (pathway in 1:length(genes.by.pathway)){
      pathway.genes <- genes.by.pathway[[pathway]]
      if (!is.na(pathway.genes)){
        list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
        list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
        scores.in.pathway <- geneList[list.genes.in.pathway]
        scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
        if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
        } else{
            p.value <- NA
        }
        new_row <- c(names(genes.by.pathway[pathway]), p.value, length(list.genes.in.pathway), sum(abs(geneLFClist[list.genes.in.pathway])>0.58&geneList[list.genes.in.pathway]<0.05), sum(geneLFClist[list.genes.in.pathway]>0.58&geneList[list.genes.in.pathway]<0.05), sum(geneLFClist[list.genes.in.pathway]< -0.58&geneList[list.genes.in.pathway]<0.05))
        pathway_pval <- rbind(pathway_pval, new_row)
      }
    }

colnames(pathway_pval) <- c("pathwayCode", "pval", "annotated", "DEG", "up", "down")
pathway_pval <- pathway_pval[complete.cases(pathway_pval),]

pathway_pval$pathwayName <- pathways.list[match(pathway_pval$pathwayCode, sub("path:","", names(pathways.list)))]

head(pathway_pval)
pathway_pval$pval <- as.numeric(pathway_pval$pval)

pathway_pval <- pathway_pval[order(pathway_pval$pval),]
head(pathway_pval)

# Write out a csv with these data
write.csv(pathway_pval, 
          file="../output/Culex_keggPathwayEnrichment_UPDown.csv", row.names = F)
