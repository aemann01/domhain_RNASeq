# Differential Gene Expression Analysis Using DESeq2 comparing Dentin cavitated teeth vs healthy teeth (H vs D)

### 1. Load environment

```bash
conda activate 2024-HIV_RNASeq
```

### 2. DESeq2 Analysis of all samples Healthy teeth vs teeth with a Dentin cavity (H v D)

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")

library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
setwd("/home/allie/domhain_RNAseq/07-ads_expression")
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)

# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
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
# [1] "number of genes with adjusted p value lower than 0.05:  595"

# out of 4208 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 537, 13%
# LFC < 0 (down)     : 58, 1.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 2396, 57%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  600"

# out of 4208 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 596, 14%
# LFC < 0 (down)     : 119, 2.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 2497, 59%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HvD.txt", quote=F, sep="\t")
# save.image()
```

### 3. PCA plot comparing H to D all ADS activity


```R
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-HvD.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-HvD.pdf")
# aliquot type?
pdf("pca_pdvpf_ADS-aliquot_type.pdf")
plotPCA(vld, intgroup=c("aliquot_type")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-aliquot_type.pdf")
# hiv group?
pdf("pca_pdvpf_ADS-hiv_status.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-hiv_status.pdf")
# no realy clustering patterns
```

### 4. Volcano plot of DE ADS genes

```R
# install.packages("BiocManager")
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(dplyr)

# add in annotations
homd <- read.table("~/domhain_RNAseq/03-star_map/homd_db/annotations.merge.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 5 genes (based on highest log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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

### 5. HIV Status Volcano Plots H v D

#### 5a. Starting with HIV

```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/allie/domhain_RNAseq/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HI",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
# filter gene count
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
# [1] "number of genes with adjusted p value lower than 0.05:  99"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 74, 5.1%
# LFC < 0 (down)     : 25, 1.7%
# outliers [1]       : 291, 20%
# low counts [2]     : 500, 34%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  91"

# out of 1455 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 158, 11%
# LFC < 0 (down)     : 33, 2.3%
# outliers [1]       : 291, 20%
# low counts [2]     : 416, 29%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HI.HvD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 5 genes (based on log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HI.HvD.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("/home/allie/.iterm2/imgcat volcano-HI.HvD.pdf")

# get object for bubble plot 
dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```

#### 5b. HEU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HEU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
# filter gene count
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
# [1] "number of genes with adjusted p value lower than 0.05:  242"

# out of 1434 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 188, 13%
# LFC < 0 (down)     : 54, 3.8%
# outliers [1]       : 371, 26%
# low counts [2]     : 322, 22%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  223"

# out of 1434 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 277, 19%
# LFC < 0 (down)     : 61, 4.3%
# outliers [1]       : 371, 26%
# low counts [2]     : 281, 20%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HEU.HvD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HEU.HvD.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HEU.HvD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```

#### 5c. HUU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HUU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
# filter gene count
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
# [1] "number of genes with adjusted p value lower than 0.05:  109"

# out of 1207 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 96, 8%
# LFC < 0 (down)     : 13, 1.1%
# outliers [1]       : 348, 29%
# low counts [2]     : 314, 26%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  96"

# out of 1207 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 131, 11%
# LFC < 0 (down)     : 14, 1.2%
# outliers [1]       : 348, 29%
# low counts [2]     : 248, 21%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HUU.HvD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HUU.HvD.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HUU.HvD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```

### 6. Create bubble chart showing average base mean across HIV status groups, split by H vs D

```R
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHHI -> HI
dfHHEU -> HEU
dfHHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "H")

# do the same thing with the disease results
dfDHI -> HI
dfDHEU -> HEU
dfDHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged2<- merged2 %>% mutate(Condition = "D")
# merge them both together!
df <- rbind(merged, merged2)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(Genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$Genus <- factor(df$Genus, levels = species_counts$Genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$Genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$Genus, df$Species, sep = " - ")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Source <- factor(df$Source, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Condition <- factor(df$Condition, levels = c("H", "D"))

pdf("HvD_bubble_plot.pdf", width = 10)
  ggplot(df, aes(x = Source, y = Species_nested, size = average_baseMean, color = Condition)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c("#22A146", "#9B002F")) +
   theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat HvD_bubble_plot.pdf")
```

### 7. Next want to compare healthy teeth but at different oral health status (HCF v HCD)

```R
# reloading data just so it doesn't screw up code down the line
# metadata
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# gene counts
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by aliquot_type
submap <- metadata[metadata$aliquot_type == "HCF" | metadata$aliquot_type == "HCD",]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
subcount <- subcount + 1
submap <- submap[order(colnames(subcount)),]
table(colnames(subcount)==submap$sample_id) # should return all true
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~aliquot_type)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$aliquot_type <- factor(star_results$aliquot_type, levels=c("HCD", "HCF"))
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
# [1] "number of genes with adjusted p value lower than 0.05:  382"

# out of 1627 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 142, 8.7%
# LFC < 0 (down)     : 240, 15%
# outliers [1]       : 0, 0%
# low counts [2]     : 221, 14%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# HCF is positive, HCD cavity negative
resLFC <- lfcShrink(se_star, coef="aliquot_type_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  374"

# out of 1627 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 186, 11%
# LFC < 0 (down)     : 274, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 187, 11%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# write results to file
write.table(resLFC, file="deseq_results_ADS-HCFvHCD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 10)
# negative top 10
top <- tail(sortdf$locus_tag, 10)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HCFvHCD.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("/home/allie/.iterm2/imgcat volcano-HCFvHCD.pdf")
# Interesting! HCD has more upregulated ADS genes -- trying to prevent acidification? Response to lower pH due to cavitated teeth elsewhere in the mouth? Need to split this by HIV status and build another bubble plot
```

### 8. HIV Status Volcano Plots HCF v HCD

#### 8a. Starting with HIV

```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/allie/domhain_RNAseq/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HI",]
submap <- submap[submap$aliquot_type == "HCF" | submap$aliquot_type == "HCD",]
# filter gene count
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~aliquot_type)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$aliquot_type <- factor(star_results$aliquot_type, levels=c("HCD", "HCF"))
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
# [1] "number of genes with adjusted p value lower than 0.05:  312"

# out of 1303 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 209, 16%
# LFC < 0 (down)     : 103, 7.9%
# outliers [1]       : 0, 0%
# low counts [2]     : 152, 12%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="aliquot_type_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  312"

# out of 1303 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 256, 20%
# LFC < 0 (down)     : 113, 8.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 227, 17%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HI.HCFvHCD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 5 genes (based on log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HI.HCFvHCD.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("/home/allie/.iterm2/imgcat volcano-HI.HCFvHCD.pdf")

# get object for bubble plot 
dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```

#### 8b. HEU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HEU",]
# filter out enamel cavity
submap <- submap[submap$aliquot_type == "HCF" | submap$aliquot_type == "HCD",]
# filter gene count
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~aliquot_type)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$aliquot_type <- factor(star_results$aliquot_type, levels=c("HCD", "HCF"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  153"

# out of 1164 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 17, 1.5%
# LFC < 0 (down)     : 136, 12%
# outliers [1]       : 367, 32%
# low counts [2]     : 42, 3.6%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="aliquot_type_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  155"

# out of 1164 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 27, 2.3%
# LFC < 0 (down)     : 151, 13%
# outliers [1]       : 367, 32%
# low counts [2]     : 91, 7.8%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HEU.HCFvHCD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HEU.HCFvHCD.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HEU.HCFvHCD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```

#### 8c. HUU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HUU",]
# filter out enamel cavity
submap <- submap[submap$aliquot_type == "HCF" | submap$aliquot_type == "HCD",]
# filter gene count
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~aliquot_type)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$aliquot_type <- factor(star_results$aliquot_type, levels=c("HCD", "HCF"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  140"

# out of 920 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 121, 13%
# LFC < 0 (down)     : 19, 2.1%
# outliers [1]       : 320, 35%
# low counts [2]     : 105, 11%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="aliquot_type_HCF_vs_HCD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  129"

# out of 920 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 140, 15%
# LFC < 0 (down)     : 31, 3.4%
# outliers [1]       : 320, 35%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HUU.HCFvHCD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HUU.HCFvHCD.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HUU.HCFvHCD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```

### 9. Create bubble chart showing average base mean across HIV status groups, split by HCF vs HCD

```R
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHHI -> HI
dfHHEU -> HEU
dfHHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "HCF")

# do the same thing with the disease results
dfDHI -> HI
dfDHEU -> HEU
dfDHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged2<- merged2 %>% mutate(Condition = "HCD")
# merge them both together!
df <- rbind(merged, merged2)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(Genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$Genus <- factor(df$Genus, levels = species_counts$Genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$Genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$Genus, df$Species, sep = " - ")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Source <- factor(df$Source, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Condition <- factor(df$Condition, levels = c("HCF", "HCD"))

pdf("HCFvHCD_bubble_plot.pdf", width = 10)
  ggplot(df, aes(x = Source, y = Species_nested, size = average_baseMean, color = Condition)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c("#22A146", "#9B002F")) +
   theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat HCFvHCD_bubble_plot.pdf")
```

### 10. Finally, want to compare E vs D among all samples (NOTE! E is negative, D positive in volcano plot)

```R
# reloading data just so it doesn't screw up code down the line
# metadata
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# gene counts
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter data
submap <- metadata[metadata$tooth_health == "E" | metadata$tooth_health == "D",]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
subcount <- subcount + 1
submap <- submap[order(colnames(subcount)),]
table(colnames(subcount)==submap$sample_id) # should return all true
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("E", "D"))

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
# adjusted p-value < 0.05
# LFC > 0 (up)       : 95, 6%
# LFC < 0 (down)     : 341, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 123, 7.8%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# E is negative, D cavity positive -- for some reason can't force factor here
resLFC <- lfcShrink(se_star, coef="tooth_health_D_vs_E", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  437"

# out of 1577 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 158, 10%
# LFC < 0 (down)     : 403, 26%
# outliers [1]       : 0, 0%
# low counts [2]     : 153, 9.7%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_ADS-DvE.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 5 genes (based on log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-DvE.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("/home/allie/.iterm2/imgcat volcano-DvE.pdf")
```

### 11. HIV Status Volcano Plots E v D

#### 11a. Starting with HIV

```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/allie/domhain_RNAseq/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HI",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "E" | submap$tooth_health == "D",]
# filter gene count
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("E", "D"))
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
# [1] "number of genes with adjusted p value lower than 0.05:  89"

# out of 951 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 38, 4%
# LFC < 0 (down)     : 51, 5.4%
# outliers [1]       : 442, 46%
# low counts [2]     : 0, 0%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_D_vs_E", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  89"

# out of 951 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 44, 4.6%
# LFC < 0 (down)     : 68, 7.2%
# outliers [1]       : 442, 46%
# low counts [2]     : 0, 0%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HI.DvE.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 5 genes (based on log fold change) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HI.DvE.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("/home/allie/.iterm2/imgcat volcano-HI.DvE.pdf")

# get object for bubble plot 
dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```

#### 11b. HEU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HEU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "E" | submap$tooth_health == "D",]
# filter gene count
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("E", "D"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  111"

# out of 1318 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 52, 3.9%
# LFC < 0 (down)     : 59, 4.5%
# outliers [1]       : 340, 26%
# low counts [2]     : 278, 21%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

resLFC <- lfcShrink(se_star, coef="tooth_health_D_vs_E", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  84"

# out of 1318 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 76, 5.8%
# LFC < 0 (down)     : 75, 5.7%
# outliers [1]       : 340, 26%
# low counts [2]     : 154, 12%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HEU.DvE.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HEU.DvE.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HEU.DvE.pdf")
dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```

#### 11c. HUU

```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HUU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "E" | submap$tooth_health == "D",]
# filter gene count
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("E", "D"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  74"

# out of 1116 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 18, 1.6%
# LFC < 0 (down)     : 56, 5%
# outliers [1]       : 336, 30%
# low counts [2]     : 216, 19%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_D_vs_E", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  67"

# out of 1116 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 21, 1.9%
# LFC < 0 (down)     : 71, 6.4%
# outliers [1]       : 336, 30%
# low counts [2]     : 150, 13%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HUU.DvE.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)

# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get new dataframe with concatenated genus species and locus id
sigsp <- paste("x", sigloc$Genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$locus_tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# create object for each genus
# remove previously created lists or will mess up
rm(list = grep("^x_", ls(), value=TRUE))
for(genus in names(siglist)){
	assign(genus, unname(siglist[[genus]]))
}

# get an object containing all genes of interest
all_genes <- do.call(c, siglist)
# add as target genes
resdf$target_gene <- rownames(resdf) %in% all_genes
# order by target gene
res_ord <- resdf[order(resdf$target_gene),]

# get a large number of colors from the viridis package to iterate over
# install.packages("viridis")
library(viridis)
keycolors <- viridis(length(siglist))

genus_to_color <- rep(NA, length(unlist(siglist)))
names(genus_to_color) <- unlist(siglist)
# map colors to genus
for (i in seq_along(siglist)){
	species_group <- siglist[[i]]
	color <- keycolors[i]
	genus_to_color[species_group] <- color
}
# add to dataframe
res_ord$color <- genus_to_color[rownames(res_ord)]
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

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
pdf("volcano-HUU.DvE.pdf", width=15, height=10)
p
dev.off()
system("/home/allie/.iterm2/imgcat volcano-HUU.DvE.pdf")
dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```

### 12. Create bubble chart showing average base mean across HIV status groups, split by D vs E

```R
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHHI -> HI
dfHHEU -> HEU
dfHHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "D")

# do the same thing with the E results
dfDHI -> HI
dfDHEU -> HEU
dfDHUU -> HUU 

# subset for clarity
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged2<- merged2 %>% mutate(Condition = "E")
# merge them both together!
df <- rbind(merged, merged2)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(Genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$Genus <- factor(df$Genus, levels = species_counts$Genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$Genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$Genus, df$Species, sep = " - ")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Source <- factor(df$Source, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Condition <- factor(df$Condition, levels = c("D", "E"))

pdf("DvE_bubble_plot.pdf", width = 10)
  ggplot(df, aes(x = Source, y = Species_nested, size = average_baseMean, color = Condition)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c("#22A146", "#9B002F")) +
   theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat DvE_bubble_plot.pdf")
```

Noting very super interesting here, some variation in the degree of ADS activity but nothing worth writing home about (as far as I can see)