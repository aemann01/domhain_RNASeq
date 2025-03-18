# For sanity sake, also going to check with DESeq2


```R
setwd("/home/allie/domhain_RNAseq/05-global_diffabund")
# install packages -- running on pickles
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")

# load required libraries
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# load data
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.table("/home/allie/domhain_RNAseq/03-star_map/02-HOMD_map/read_counts.txt", header=T, sep="\t", row.names = 1)

# format for deseq
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
# only compare HCF to HCD
metadata <- metadata[metadata$aliquot_type == "HCF" | metadata$aliquot_type == "HCD",]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")  
# only keep columns found in metadata
genecounts <- genecounts[, colnames(genecounts) %in% row.names(metadata)]
# reorder by metadata rownames
metadata <- metadata[order(colnames(genecounts)),]
table(colnames(genecounts)==metadata$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = genecounts, colData = metadata, design = ~aliquot_type)
# filter out any genes with fewer than 10 reads total
star_results <- star_results[rowSums(counts(star_results)) >= 10,]
star_results
# set factor level (this determines which direction the comparisions are made -- by default it's by alphabetical order)
star_results$aliquot_type <- factor(star_results$aliquot_type, levels=c("HCD", "HCF"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results)
proc.time() - ptm # this takes about an hour to complete
# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)

# filter out low count genes
resLFC <- lfcShrink(se_star, coef="aliquot_type_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# write results to file
write.table(resLFC, file="deseq_results.txt", quote=F, sep="\t")

# transform for visualizations
vld <- vst(se_star)

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 25)
 
#make a subset of the log transformed counts for just the top 25 varying genes
top25Counts <- assay(vld)[topVarGenes,]
write.csv(top25Counts, file="top25counts.vld.csv", quote=FALSE)
 
#PLOT PCA
#PCA using top 500 varying genes
pdf("pca_hcfvhcd.pdf")
plotPCA(vld, intgroup=c("aliquot_type"), ntop=500) + theme_minimal()
dev.off()

df <- as.data.frame(colData(vld)[,c("hiv_status","aliquot_type")])
colnames(df) <- c("hiv_status", "aliquot_type")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
library(RColorBrewer)
x <- pheatmap(top25Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_deseq.pdf")
```

Ok so a lot less differentiation in these samples as compared to UF -- try with HIV status and HCD?

```R
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
genecounts <- read.table("/home/allie/domhain_RNAseq/03-star_map/02-HOMD_map/read_counts.txt", header=T, sep="\t", row.names = 1)






# combine metadata categories
metadata$hiv_aliquot <- paste(metadata$hiv_status, metadata$aliquot_type, sep="_")
# format for deseq
# only compare HCF to HCD
metadata <- metadata[metadata$hiv_aliquot == "HI_HCD" | metadata$hiv_aliquot == "HUU_HCD",]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")  
# only keep columns found in metadata
genecounts <- genecounts[, colnames(genecounts) %in% row.names(metadata)]
# reorder by metadata rownames
metadata <- metadata[order(colnames(genecounts)),]
table(colnames(genecounts)==metadata$sample_id) # should return all true












# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = genecounts, colData = metadata, design = ~hiv_aliquot)
# filter out any genes with fewer than 10 reads total
star_results <- star_results[rowSums(counts(star_results)) >= 10,]
star_results
# set factor level (this determines which direction the comparisions are made -- by default it's by alphabetical order)
star_results$hiv_aliquot <- factor(star_results$hiv_aliquot, levels=c("HCD", "HCF"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results)
proc.time() - ptm # this takes about an hour to complete
# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)

# filter out low count genes
resLFC <- lfcShrink(se_star, coef="hiv_aliquot_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# write results to file
write.table(resLFC, file="deseq_results.txt", quote=F, sep="\t")

# transform for visualizations
vld <- vst(se_star)

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 25)
 
#make a subset of the log transformed counts for just the top 25 varying genes
top25Counts <- assay(vld)[topVarGenes,]
write.csv(top25Counts, file="top25counts.vld.csv", quote=FALSE)
 
#PLOT PCA
#PCA using top 500 varying genes
pdf("pca_hcfvhcd.pdf")
plotPCA(vld, intgroup=c("hiv_aliquot"), ntop=500) + theme_minimal()
dev.off()

df <- as.data.frame(colData(vld)[,c("hiv_status","hiv_aliquot")])
colnames(df) <- c("hiv_status", "hiv_aliquot")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
library(RColorBrewer)
x <- pheatmap(top25Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_deseq.pdf")
```

Save environment image

```R
save.image("deseq_global.RData")
```