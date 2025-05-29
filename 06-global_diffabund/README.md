Global DESeq2 Analysis of HIV status groups

Load conda environment

```bash
conda activate 2024-HIV_RNASeq
```

```R
# load data to reduce rerun time
load("deseq_global.RData")
```

First comparing H to D

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")

library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
setwd("/home/suzanne/rna_dohmain/06-global_diffabund")
metadata <- read.table("/home/suzanne/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("/home/suzanne/rna_dohmain/homd_map/read_counts.txt", header=T, sep="\t", row.names=1)

# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare HEU to HUU
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "H"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm

# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)

# H is positive, D cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)

# write results to file
<<<<<<< Updated upstream
write.table(resLFC, file="deseq_results_global-HvD.txt", quote=F, sep="\t")
save.image("deseq_global.RData")
=======
write.table(resLFC, file="deseq_results_global-HUUvHEU.txt", quote=F, sep="\t")
save.image("deseq_global-HUUvHEU.RData")
>>>>>>> Stashed changes
```

PCA plot comparing H to D all ADS activity

```R
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_global-HvD.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_global-HvD.pdf")
```

Heatmap of global RNA DE expression 

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")
# install.packages("pheatmap")
library(pheatmap, warn.conflicts = F, quietly = T)
library(RColorBrewer)

df <- as.data.frame(colData(vld)[,c("tooth_health","hiv_status")])
colnames(df) <- c("tooth_health", "hiv_status")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 25)
 
#make a subset of the log transformed counts for just the top 25 varying genes
top25Counts <- assay(vld)[topVarGenes,]
write.csv(top25Counts, file="top25counts.vld.csv", quote=FALSE)

x <- pheatmap(top25Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_HvD.pdf")
system("/home/allie/.iterm2/imgcat heatmap_HvD.pdf")
```

Now I want to do the same thing but compare HIV status categories (the H and D are not that differentially distributed in this dataset)

```R
# filter metadata so that we only compare HEU to HUU
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm

# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)

resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)

# write results to file
write.table(resLFC, file="deseq_results_global-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_global.RData")

df <- as.data.frame(colData(vld)[,c("tooth_health","hiv_status")])
colnames(df) <- c("tooth_health", "hiv_status")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 25)

#make a subset of the log transformed counts for just the top 25 varying genes
top25Counts <- assay(vld)[topVarGenes,]
write.csv(top25Counts, file="top25counts.vld.csv", quote=FALSE)

x <- pheatmap(top25Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_HIvHUU.pdf")
system("/home/allie/.iterm2/imgcat heatmap_HIvHUU.pdf")
```


<<<<<<< Updated upstream
=======
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low10 <- head(sortdf$locus_tag, 10)
# negative top 10
top10 <- tail(sortdf$locus_tag, 10)
# concatenate
labgenes <- c(top10, low10)

# Create volcano plot
p <- EnhancedVolcano(res_ord,
	lab = ifelse(res_ord$locus_tag %in% labgenes, paste(res_ord$Species, res_ord$gene, sep=" "), ""),
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	max.overlaps = Inf,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HvD.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HvD.pdf")
```
# 1. Compare HUU vs HI global abundance
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/global")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("~/rna_dohmain/homd_map/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status-HUUvsHI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)

write.table(resLFC, file="deseq_results_global-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_global-HUUvHI.RData")
>>>>>>> Stashed changes
