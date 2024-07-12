## Differential Gene Expression Analysis Using DESeq2 comparing Dentin cavitated teeth vs healthy teeth (H vs D)

Load environment

```bash
conda activate 2024-HIV_RNASeq
```

DESeq2 Analysis of all samples Healthy teeth vs teeth with a Dentin cavity (H v D)

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

PCA plot comparing H to D all ADS activity

```R
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-HvD.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-HvD.pdf")
```

Volcano plot of DE ADS genes

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

Now I want to make a companion figure for the volcano plot that illustrates the number of significant genes in each species

```R
# treemap plot
# install.packages("treemap")
library(treemap)

# first need to split our target genes by whether they are more highly expressed in health or disease
dfDall <- subset(sortdf, log2FoldChange < 0)
dfHall <- subset(sortdf, log2FoldChange > 0)

# get counts of genus + species
posgroup <- dfpos %>% group_by(Genus, Species) %>% summarise(count = n(), hexcol = first(color), .groups = 'drop')
neggroup <- dfneg %>% group_by(Genus, Species) %>% summarise(count = n(), hexcol = first(color), .groups = 'drop')

# create color map for each dataframe
colors <- setNames(posgroup$hexcol, paste(posgroup$Genus, posgroup$Species, sep="_"))

# create treemap 
pdf("postreemap.pdf")
treemap(posgroup,
        index = c("Genus", "Species"),  # Hierarchical index
        vSize = "count",                # Size of the rectangles
        vColor = "hexcol", 
        type = "color",
        palette = colors,              # Color palette
        border.col = "white",           # Border color of the rectangles
        fontsize.labels = c(15, 10),    # Font size for labels at different levels
        fontcolor.labels = c("black", "black"),  # Font color for labels
        fontface.labels = c(2, 1),      # Font face for labels
        bg.labels = "#CCCCCCDC",      # Background color for labels
        align.labels = list(c("left", "top"), c("center", "center")))
dev.off()
system("/home/allie/.iterm2/imgcat postreemap.pdf")
# treemap for negative group
colors <- setNames(neggroup$hexcol, paste(neggroup$Genus, neggroup$Species, sep="_"))
pdf("negtreemap.pdf")
treemap(neggroup,
        index = c("Genus", "Species"),  # Hierarchical index
        vSize = "count",                # Size of the rectangles
        vColor = "hexcol", 
        type = "color",
        palette = colors,              # Color palette
        border.col = "white",           # Border color of the rectangles
        fontsize.labels = c(15, 10),    # Font size for labels at different levels
        fontcolor.labels = c("black", "black"),  # Font color for labels
        fontface.labels = c(2, 1),      # Font face for labels
        bg.labels = "#CCCCCCDC",      # Background color for labels
        align.labels = list(c("left", "top"), c("center", "center")))
dev.off()
system("/home/allie/.iterm2/imgcat negtreemap.pdf")
```

I'm not a super fan of this visualization -- will use this code to create a bubble chart below

HIV Status Volcano Plots H v D
--------------------------------------

Now filter dataframe by HIV status and plot individual volcano plots for each HIV status group (HvD)

Starting with HIV

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

HEU

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

HUU

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

What about a bubble chart?

```R
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHall -> all
dfHHI -> HI
dfHHEU -> HEU
dfHHUU -> HUU 

# subset for clarity
all <- all %>% select('baseMean', 'Genus', 'Species')
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_all <- all %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_all <- group_all %>% mutate(Source = "All")
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_all, group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "H")

# do the same thing with the disease results
dfDall -> all
dfDHI -> HI
dfDHEU -> HEU
dfDHUU -> HUU 

# subset for clarity
all <- all %>% select('baseMean', 'Genus', 'Species')
HI <- HI %>% select('baseMean', 'Genus', 'Species')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species')

# group and get average base mean
group_all <- all %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_all <- group_all %>% mutate(Source = "All")
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_all, group_HI, group_HEU, group_HUU)
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
df$Source <- factor(df$Source, levels = c("All", "HUU", "HEU", "HI"))
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
    size = "Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c("#22A146", "#9B002F")) +
   theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat HvD_bubble_plot.pdf")
```






