# 1. HvD
## 1.1 HI HvD
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/suzanne/rna_dohmain/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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

write.table(resdf, file="deseq_results_ADS-HI.HvD.txt", quote=F, sep="\t")

############### volcano plot
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
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
system("~/.iterm2/imgcat volcano-HI.HvD.pdf")

# get object for bubble plot 
dfHI <- sortdf
dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```
## 1.2. HEU HvD
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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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

write.table(resdf, file="deseq_results_ADS-HEU.HvD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(combined_df)
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
system("~/.iterm2/imgcat volcano-HEU.HvD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfHEU <- sortdf
dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```
## 1.3. HUU HvD
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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
# write results to file
resdf <- as.data.frame(combined_df)
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

write.table(resdf, file="deseq_results_ADS-HUU.HvD.txt", quote=F, sep="\t")

############### volcano plot
resdf <- as.data.frame(combined_df)
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
system("~/.iterm2/imgcat volcano-HUU.HvD.pdf")

# first need to split our target genes by whether they are more highly expressed in health or disease
dfHUU <- sortdf
dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```
## 1.4. Make bubble plot
```R
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHI -> HI
dfHEU -> HEU
dfHUU -> HUU 

# subset for clarity
HI <- HI %>% select('H', 'Genus', 'Species')
HEU <- HEU %>% select('H', 'Genus', 'Species')
HUU <- HUU %>% select('H', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(H))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(H))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(H))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "H")

# do the same thing with the disease results
dfHI -> HI
dfHEU -> HEU
dfHUU -> HUU 

# subset for clarity
HI <- HI %>% select('D', 'Genus', 'Species')
HEU <- HEU %>% select('D', 'Genus', 'Species')
HUU <- HUU %>% select('D', 'Genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(D))
group_HEU <- HEU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(D))
group_HUU <- HUU %>% group_by(Genus, Species) %>% summarise(average_baseMean = mean(D))

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
system("~/.iterm2/imgcat HvD_bubble_plot.pdf")
```
# 2. Using enamel
## 2.2 HI HvEvD
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/suzanne/rna_dohmain/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
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
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E" | submap$tooth_health == "D",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "E", "H"))

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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-HI.HED.txt", quote=F, sep="\t")

############### volcano plot
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
system("~/.iterm2/imgcat volcano-HI.HvD.pdf")

# get object for bubble plot 
dfHI <- res_ord
dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```
## 2.2. HEU HvEvD
```R
# filter by HIV group
submap <- metadata[metadata$hiv_status == "HEU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E" | submap$tooth_health == "D",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "E", "H"))

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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-HEU.HED.txt", quote=F, sep="\t")

############### volcano plot
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
pdf("volcano-HEU.HvD.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-HEU.HvD.pdf")

# get object for bubble plot 
dfHEU <- res_ord
dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```
## 2.3. HUU HvEvD
```R
submap <- metadata[metadata$hiv_status == "HUU",]
# filter out enamel cavity
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E" | submap$tooth_health == "D",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "E", "H"))

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

baseMeanPerLvl <- sapply( levels(se_star$tooth_health), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$tooth_health== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-HUU.HED.txt", quote=F, sep="\t")

############### volcano plot
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
pdf("volcano-HUU.HvD.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-HUU.HvD.pdf")

# get object for bubble plot 
dfHUU <- res_ord
dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```
## 2.4. Make bubble plot
```R
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHI -> HI
dfHEU -> HEU
dfHUU -> HUU 

# subset for clarity
HI <- HI %>% select('H', 'genus', 'Species')
HEU <- HEU %>% select('H', 'genus', 'Species')
HUU <- HUU %>% select('H', 'genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(H))
group_HEU <- HEU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(H))
group_HUU <- HUU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(H))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "H")

# rename the input dataframes so that you don't have to run all the code above
dfHI -> HI
dfHEU -> HEU
dfHUU -> HUU 

# subset for clarity
HI <- HI %>% select('E', 'genus', 'Species')
HEU <- HEU %>% select('E', 'genus', 'Species')
HUU <- HUU %>% select('E', 'genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(E))
group_HEU <- HEU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(E))
group_HUU <- HUU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(E))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged2 <- merged2 %>% mutate(Condition = "E")

# do the same thing with the disease results
dfHI -> HI
dfHEU -> HEU
dfHUU -> HUU 

# subset for clarity
HI <- HI %>% select('D', 'genus', 'Species')
HEU <- HEU %>% select('D', 'genus', 'Species')
HUU <- HUU %>% select('D', 'genus', 'Species')

# group and get average base mean
group_HI <- HI %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(D))
group_HEU <- HEU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(D))
group_HUU <- HUU %>% group_by(genus, Species) %>% summarise(average_baseMean = mean(D))

# add source column so you can track them back
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged3 <- rbind(group_HI, group_HEU, group_HUU)
# add health indicator
merged3<- merged3 %>% mutate(Condition = "D")
# merge them both together!
df <- rbind(merged, merged2, merged3)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$genus <- factor(df$genus, levels = species_counts$genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$genus, df$Species, sep = " - ")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Source <- factor(df$Source, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Condition <- factor(df$Condition, levels = c("H", "E", "D"))

df_colored <- df %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()


df_colored <- df %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source, na.rm = TRUE),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",
      rank_within_source == max_rank ~ "lowest",
      TRUE ~ "middle"
    ) %>% as.character()  # Ensure output is character
  ) %>%
  ungroup()

df_colored$rank_within_source <- factor(df_colored$rank_within_source, levels=c("highest", "middle", "lowest"))


df_colored <- df_strep %>%
  group_by(Condition, Species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Genus)))

df_colored$rank_within_source <- factor(df_colored$rank_within_source, levels=c("highest", "middle", "lowest"))

pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_colored, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")


df_colored <- df %>%
  group_by(Source, Species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_filt <- filter(df_colored, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")
df_filt <- filter(df_filt, genus == "Streptococcus")


df_filt2 <- df_filt %>%
  group_by(Species, Source) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.flipped.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.flipped.pdf")

df_filt2 <- df_filt %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.hiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.hiv.pdf")


pdf("HED_bubble_plot.healthvhiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.healthvhiv.pdf")

df_colored <- df %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()
df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))
df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Genus)))


pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_strep, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4596C3",
      "middle" = "#4596C3",
      "lowest" = "#4596C3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")
```
# 3. rpoC with enamel
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)
setwd("~/rna_dohmain/07-ads_expression")
# load and clean up rpoC data
map <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
map$aliquot_type <- sub("-", "", map$aliquot_type)
row.names(map) <- map$sample_id
# sequence table
seqtab <- read.table("~/rna_dohmain/11-rpoc_processing/sequence_table.merged.txt", header=T, sep="\t", row.names=1)
tax <- read.table("~/rna_dohmain/11-rpoc_processing/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
notinmeta <- setdiff(row.names(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), row.names(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw
# should both come back as character(0) 
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
# change sp. HMT to oral taxon to match nomenclature differences
filter_sp <- c("Streptococcus_anginosus", "Streptococcus_australis", "Streptococcus_cristatus", "Streptococcus_gordonii", "Streptococcus_sanguinis", "Streptococcus_sp._oral_taxon_056", "Streptococcus_mitis", "Streptococcus_oralis", "Streptococcus_parasanguinis", "Streptococcus_salivarius", "Streptococcus_sp._oral_taxon_066", "Streptococcus_constellatus", "Streptococcus_intermedius", "Actinomyces_naeslundii", "Actinomyces_oris", "Actinomyces_johnsonii", "Actinomyces_sp._oral_taxon_169", "Actinomyces_sp._oral_taxon_170", "Actinomyces_sp._oral_taxon_175", "Actinomyces_viscosus", "Cryptobacterium_curtum", "Cutibacterium_acnes", "Kingella_oralis", "Oribacterium_asaccharolyticum")
# do these actually exist in my phyloseq object?
species_df <- data.frame(QuerySpecies = filter_sp, stringsAsFactors = FALSE)
# Fuzzy match with taxonomy table
tax_levels <- paste0("V", 2:13)
# Find out if these species exist anywhere in the taxonomy and pull asvids
tax_df <- as.data.frame(tax_table(ps.dat)) %>%
  rownames_to_column("ASVID")

# Find matching ASVs for each query species
result_list <- lapply(species_df$QuerySpecies, function(x) {
  # Find rows where the species appears in any taxonomic level
  matches <- tax_df %>%
    filter(if_any(all_of(tax_levels), ~ . == x))
  
  if(nrow(matches) > 0) {
    data.frame(QuerySpecies = x, 
              ASVID = matches$ASVID,
              TaxLevel = apply(matches[, tax_levels], 1, function(row) {
                names(which(row == x))[1]
              }))
  } else {
    NULL
  }
})

# Combine results
matches_df <- bind_rows(result_list)
# Get all unique matching ASV IDs
matching_asvs <- unique(matches_df$ASVID)
# filter by matching ASVIDs
ps.dat.filt <- prune_taxa(matching_asvs, ps.dat)
ps.dat.filt
# head(tax_table(ps.dat.filt))
# first need to get a standardized species column (since the species show up in different levels)
matches_df <- matches_df %>%
  mutate(StandardizedSpecies = QuerySpecies)

# First calculate total number of samples in each comparison group
sample_counts <- sample_data(ps.dat.filt) %>% 
  as_tibble() %>% 
  count(hiv_status, tooth_health, name = "total_samples")
# Conditional mean abundance
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status", "tooth_health")) %>%
  # Calculate metrics
  group_by(Species, hiv_status, tooth_health) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))

df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)

# I want to add in missing species that do not show up in our rpoC data to make the figures comparable
df <- df %>%
  # Ensure all species are included (even missing ones)
  complete(
    Species = filter_sp,
    hiv_status = unique(df$hiv_status),
    tooth_health = unique(df$tooth_health),
    fill = list(conditional_mean = 0, prevalence = 0)
  ) %>%
  # Reapply factor levels to Species
  mutate(Species = factor(Species, levels = filter_sp))

# order by species level
df <- df[order(df$Species),]
# set levels of x axis
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))
# df <- df %>%
#   filter(tooth_health != "E")
# and order of grid
df$tooth_health <- factor(df$tooth_health, levels = c("H", "E", "D"))

df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))

df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Species)))

df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Species)))
df_strep$Species <- factor(df_strep$Species, levels=c("Streptococcus_anginosus", "Streptococcus_australis", "Streptococcus_constellatus", "Streptococcus_cristatus", "Streptococcus_gordonii",
	"Streptococcus_oralis", "Streptococcus_parasanguinis", "Streptococcus_sanguinis", "Streptococcus_sp._oral_taxon_056", "Streptococcus_intermedius", "Streptococcus_mitis",
	 "Streptococcus_salivarius", "Streptococcus_sp._oral_taxon_066"))

pdf("HED_rpoC_bubble_plot.pdf", width = 10)
ggplot(df_strep,
  aes(
    x = hiv_status,
    y = Species,
    size = log10(conditional_mean + 1),
    color = tooth_health
  )
) +
  geom_point(alpha = ifelse(df_strep$conditional_mean == 0, 0.1, 0.7)) +
  scale_size(range = c(0.5, 10), name = "Log10 Conditional Mean") +  
  scale_y_discrete(limits = rev) +  # Reverse y-axis order (top-to-bottom)
  scale_color_manual(values = c("#22A146", "#22A146", "#22A146")) +
  labs(
    x = "Source",
    y = "Species",
    color = "Tooth Health"
  ) +
  facet_grid(. ~ tooth_health, switch = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Improve x-axis label readability
  )
dev.off()
system("~/.iterm2/imgcat HED_rpoC_bubble_plot.pdf")


df_colored <- df %>%
  group_by(Species, hiv_status) %>%
  mutate(
    rank_within_source = rank(-conditional_mean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()
df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))
df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Species)))
df_strep$Species <- factor(df_strep$Species, levels=c("Streptococcus_anginosus", "Streptococcus_australis", "Streptococcus_constellatus", "Streptococcus_cristatus", "Streptococcus_gordonii",
	"Streptococcus_oralis", "Streptococcus_parasanguinis", "Streptococcus_sanguinis", "Streptococcus_sp._HMT_056", "Streptococcus_intermedius", "Streptococcus_mitis",
	 "Streptococcus_salivarius", "Streptococcus_sp._HMT_066"))

pdf("HED_rpoC_bubble_plot.pdf", width = 10)
ggplot(df_strep,
  aes(
    x = hiv_status,
    y = Species,
    size = log10(conditional_mean + 1),
    color = rank_color
  )
) +
  geom_point(alpha = ifelse(df_colored$conditional_mean == 0, 0.1, 0.7)) +
  scale_size(range = c(0.5, 10), name = "Log10 Conditional Mean") +  
  scale_y_discrete(limits = rev) +  # Reverse y-axis order (top-to-bottom)
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
  labs(
    x = "Source",
    y = "Species",
    color = ""
  ) +
  facet_grid(. ~ tooth_health, switch = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Improve x-axis label readability
  )
dev.off()
system("~/.iterm2/imgcat HED_rpoC_bubble_plot.pdf")
```
# 4. rpoC long pcoAs
```R
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
setwd("~/long_oral/diversity")
load("~/long_oral/master_phyloseq.RData")
# set up some color palettes for major groupings
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
sample_data(ps.dat)$hiv_status <- factor(sample_data(ps.dat)$hiv_status, levels=c("HI", "HEU", "HUU"))
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~hiv_status)
# capscale plot by HIV status group
pdf("./bdiv_cap.hiv_status.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="hiv_status") + 
    theme_minimal() + 
    scale_color_manual(values=hivCols)
dev.off()

ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~aliquot_type)
pdf("./bdiv_cap.aliquot_type.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="aliquot_type") + 
    theme_minimal() +
    scale_color_manual(values=healthCols)
dev.off()
```
# 5. HvD
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
setwd("~/rna_dohmain/07-ads_expression")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
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
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(dplyr)

# add in annotations
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
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


colormap <- setNames(res_ord$color, res_ord$Genus)

res_ord$color2 <- ifelse(res_ord$Genus == "Actinomyces", "#F26F5F",
ifelse(res_ord$Genus == "Cryptobacterium", "#791A92",
ifelse(res_ord$Genus == "Cutibacterium", "#F4E6BF",
ifelse(res_ord$Genus == "Oribacterium", "#235EC7",
ifelse(res_ord$Genus == "Peptostreptococcaceae", "#F1F046",
ifelse(res_ord$Genus == "Streptococcus", "#53B36D",
"#808080"))))))
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color2)] <- NA
# change NA genus to grey
res_ord$color2[is.na(res_ord$color2)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_ord$color2, res_ord$Genus)
# Create volcano plot
p <- EnhancedVolcano(res_ord,
	lab = ifelse(res_ord$locus_tag == "bleep", paste(res_ord$Species, res_ord$gene, sep=" "), ""),
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap2 ,
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
system("~/.iterm2/imgcat volcano-HvD.pdf")
```
# 6. HIV status compare
## 6.2 H HIvHEUvHUU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/suzanne/rna_dohmain/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-H.hiv.txt", quote=F, sep="\t")

############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
filter(sigloc, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")

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
pdf("volcano-H.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-H.HUUvHI.pdf")

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(res_ord)[1:5] <- paste0(colnames(res_ord)[1:5], "_HUUvHI")
head(res_ord)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$locus_tag <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$locus_tag<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(res_ord), as.data.frame(resHEUvHI), by = "locus_tag")
dfH <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "locus_tag")
```
## 6.2. E hiv status
```R
# filter by HIV group
submap <- metadata[metadata$tooth_health == "E",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-E.hiv.txt", quote=F, sep="\t")

############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
filter(sigloc, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")

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
pdf("volcano-E.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-E.HUUvHI.pdf")

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(res_ord)[1:5] <- paste0(colnames(res_ord)[1:5], "_HUUvHI")
head(res_ord)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$locus_tag <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$locus_tag<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(res_ord), as.data.frame(resHEUvHI), by = "locus_tag")
dfE <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "locus_tag")
```
## 6.3. HUU HvEvD
```R
# filter by HIV group
submap <- metadata[metadata$tooth_health == "D",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
resdf <- as.data.frame(combined_df)
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
resdf$genus <- resdf$Genus

write.table(resdf, file="deseq_results_ADS-D.hiv.txt", quote=F, sep="\t")

############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
filter(sigloc, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")

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
pdf("volcano-D.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-D.HUUvHI.pdf")

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(res_ord)[1:5] <- paste0(colnames(res_ord)[1:5], "_HUUvHI")
head(res_ord)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$locus_tag <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$locus_tag<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(res_ord), as.data.frame(resHEUvHI), by = "locus_tag")
dfD <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "locus_tag")
```
## 6.4. Make bubble plot
```R
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HUU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HUU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HUU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )
group_E <- E %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )

# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged <- rbind(group_H, group_E, group_D)
# add health indicator
merged <- merged %>% mutate(Condition = "HUU")

# rename the input dataframes so that you don't have to run all the code above
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HEU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HEU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HEU', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
group_E <- E %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged2 <- rbind(group_H, group_E, group_D)
# add health indicator
merged2 <- merged2 %>% mutate(Condition = "HEU")

# do the same thing with the disease results
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HI', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HI', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HI', 'genus', 'Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )
group_E <- E  %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(genus, Species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )

# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged3 <- rbind(group_H, group_E, group_D)
# add health indicator
merged3 <- merged3 %>% mutate(Condition = "HI")

# merge them both together!
df <- rbind(merged, merged2, merged3)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$genus <- factor(df$genus, levels = species_counts$genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$genus, df$Species, sep = " - ")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Condition <- factor(df$Condition, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Source <- factor(df$Source, levels = c("H", "E", "D"))

# order within HIV status
df_colored <- df %>%
  group_by(Species, Source) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source, na.rm = TRUE),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",
      rank_within_source == max_rank ~ "lowest",
      TRUE ~ "middle"
    ) %>% as.character()  # Ensure output is character
  ) %>%
  ungroup()

df_colored$rank_within_source <- factor(df_colored$rank_within_source, levels=c("highest", "middle", "lowest"))
df_filt <- filter(df_colored, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")
df_filt <- filter(df_filt, genus == "Streptococcus")

df_colored <- df_filt %>%
  group_by(Source, Species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_colored$upregulation <- ifelse(df_colored$any_significant_HUUvHI == "TRUE", )
pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_colored, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")








df_colored <- df %>%
  group_by(Source, Species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_filt <- filter(df_colored, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")
df_filt <- filter(df_filt, genus == "Streptococcus")


df_filt2 <- df_filt %>%
  group_by(Species, Source) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.flipped.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.flipped.pdf")

df_filt2 <- df_filt %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.hiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.hiv.pdf")


pdf("HED_bubble_plot.healthvhiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.healthvhiv.pdf")

df_colored <- df %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()
df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))
df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Genus)))


pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_strep, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4596C3",
      "middle" = "#4596C3",
      "lowest" = "#4596C3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")
```
# 7. HvD American
```sh
awk -F'\t' 'BEGIN { OFS=FS }
NR == 1 {
  print
  next
}
{
  if ($2 == "PD") $2 = "D"
  print
}' map.txt > temp

awk -F'\t' 'BEGIN { OFS=FS }
NR == 1 {
  print
  next
}
{
  if ($2 == "PE") $2 = "E"
  print
}' temp > temp2

awk -F'\t' 'BEGIN { OFS=FS }
NR == 1 {
  print
  next
}
{
  if ($2 == "PF") $2 = "H"
  print
}' temp2 > map.txt
sed -i 's/tooth_type/tooth_health/' map.txt

```
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
setwd("~/ads_plaque/10-vince")
metadata <- read.table("~/ads_plaque/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("ads.read_counts.txt", header=T, sep="\t", row.names = 1)

# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "X", replacement = "UF")  
genecounts <- genecounts[1:(length(genecounts)-1)]

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
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(dplyr)

# add in annotations
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
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


colormap <- setNames(res_ord$color, res_ord$Genus)

res_ord$color2 <- ifelse(res_ord$Genus == "Actinomyces", "#F26F5F",
ifelse(res_ord$Genus == "Cryptobacterium", "#791A92",
ifelse(res_ord$Genus == "Cutibacterium", "#F4E6BF",
ifelse(res_ord$Genus == "Oribacterium", "#235EC7",
ifelse(res_ord$Genus == "Peptostreptococcaceae", "#F1F046",
ifelse(res_ord$Genus == "Streptococcus", "#53B36D",
"#808080"))))))
# if no color, remove genus label
res_ord$Genus[is.na(res_ord$color2)] <- NA
# change NA genus to grey
res_ord$color2[is.na(res_ord$color2)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_ord$color2, res_ord$Genus)
# Create volcano plot
p <- EnhancedVolcano(res_ord,
	lab = ifelse(res_ord$locus_tag == "bleep", paste(res_ord$Species, res_ord$gene, sep=" "), ""),
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap2 ,
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
pdf("volcano-HvD.USA.pdf", width=15, height=10)
p
dev.off()
system("~/.iterm2/imgcat volcano-HvD.USA.pdf")
```
# 6. HIV status compare callopsed by species
```py
import pandas as pd
import numpy as np  # For handling NaN values

# read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/07-ads_expression/arcGene_read_counts.cleaned.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.rename(columns={"Unnamed: 0": "Geneid"})
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = '/home/suzanne/rna_dohmain/homd_map/annotations.arc.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, gene_base, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus_Species

# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / 1)
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads
# select only columns that start with arcA, arcB, or arcC
genes_of_interest = [col for col in percentage_reads_by_sample.columns 
                     if col.startswith(('arcA', 'arcB', 'arcC'))]

# Subset the dataframe to those columns
subset_df = percentage_reads_by_sample[genes_of_interest]
subset_df['row_sum'] = subset_df.sum(axis=1)

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('species_reads.arcABC.txt', sep="\t", index=True, header = True)
```
## 6.2 H HIvHEUvHUU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(dplyr)
library(EnhancedVolcano)

setwd("/home/suzanne/rna_dohmain/07-ads_expression")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.csv("./species_reads.arcABC.txt", header=T, sep="\t", row.names=1)
genecounts <- t(genecounts)
# fix sample names in gene counts so they match the metadata
# colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL


############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- combined_df %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(combined_df)[1:5] <- paste0(colnames(combined_df)[1:5], "_HUUvHI")
combined_df$species <- rownames(combined_df)
head(combined_df)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$species <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$species<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(combined_df), as.data.frame(resHEUvHI), by = "species")
dfH <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "species")
```
## 6.2. E hiv status
```R
# filter by HIV group
submap <- metadata[metadata$tooth_health == "E",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL


############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- combined_df %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(combined_df)[1:5] <- paste0(colnames(combined_df)[1:5], "_HUUvHI")
combined_df$species <- rownames(combined_df)
head(combined_df)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$species <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$species<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(combined_df), as.data.frame(resHEUvHI), by = "species")
dfE <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "species")
```
## 6.3. HUU HvEvD
```R
submap <- metadata[metadata$tooth_health == "D",]
# filter out enamel cavity
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU" | submap$hiv_status == "HI",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

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
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
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

baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL


############### volcano plot
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- combined_df %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# get object for bubble plot 
resHUUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")


# get heu v huu
submap <- submap[submap$hiv_status == "HUU" | submap$hiv_status == "HEU",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c( "HEU", "HUU"))

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
resHUUvHEU <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")

colnames(combined_df)[1:5] <- paste0(colnames(combined_df)[1:5], "_HUUvHI")
combined_df$species <- rownames(combined_df)
head(combined_df)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$species <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$species<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(combined_df), as.data.frame(resHEUvHI), by = "species")
dfD <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "species")
```
## 6.4. Make bubble plot
```R
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HUU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HUU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HUU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )
group_E <- E %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HUU, na.rm = TRUE)
  )

# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged <- rbind(group_H, group_E, group_D)
# add health indicator
merged <- merged %>% mutate(Condition = "HUU")

# rename the input dataframes so that you don't have to run all the code above
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HEU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HEU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HEU', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
group_E <- E %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HEU, na.rm = TRUE)
  )
# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged2 <- rbind(group_H, group_E, group_D)
# add health indicator
merged2 <- merged2 %>% mutate(Condition = "HEU")

# do the same thing with the disease results
dfH -> H
dfE -> E
dfD -> D 

# subset for clarity
H <- H %>% select('HI', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
E <- E %>% select('HI', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
D <- D %>% select('HI', 'species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_H <- H %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )
group_E <- E  %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )
group_D <- D %>%
  group_by(species) %>%
  summarise(
    direction_HUUvHI = case_when(
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) &
        any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HUUvHEU = case_when(
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) &
        any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "both",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE) ~ "up",
      any(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),
    
    direction_HEUvHI = case_when(
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) &
        any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "both",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE) ~ "up",
      any(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE) ~ "down",
      TRUE ~ "not_significant"
    ),

    average_baseMean = mean(HI, na.rm = TRUE)
  )

# add source column so you can track them back
group_H <- group_H %>% mutate(Source = "H")
group_E <- group_E %>% mutate(Source = "E")
group_D <- group_D %>% mutate(Source = "D")

# merge all into one
merged3 <- rbind(group_H, group_E, group_D)
# add health indicator
merged3 <- merged3 %>% mutate(Condition = "HI")

# merge them both together!
df <- rbind(merged, merged2, merged3)

# now that the data is formatted properly, can make a bubble plot
df$species <- as.factor(df$species)
# Order dataframe based on Genus factor levels
df <- df[order(df$species), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$species)
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Condition <- factor(df$Condition, levels = c("HUU", "HEU", "HI"))
# and order of grid
df$Source <- factor(df$Source, levels = c("H", "E", "D"))

# order within HIV status
df_colored <- df %>%
  group_by(species, Source) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source, na.rm = TRUE),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",
      rank_within_source == max_rank ~ "lowest",
      TRUE ~ "middle"
    ) %>% as.character()  # Ensure output is character
  ) %>%
  ungroup()

df_colored$rank_within_source <- factor(df_colored$rank_within_source, levels=c("highest", "middle", "lowest"))
df_filt <- filter(df_colored, species == "Streptococcus_oralis" | species =="Streptococcus_cristatus" | species =="Streptococcus_gordonii" | species == "Streptococcus_parasanguinis" | species =="Streptococcus_sanguinis"| species == "Streptococcus_mitis" | species == "Streptococcus_sinensis")
# df_strep <- df_colored %>% filter(grepl("^Strep", as.character(species)))

df_colored <- df_filt %>%
  group_by(Source, species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_colored$upregulation <- ifelse(df_colored$any_significant_HUUvHI == "TRUE", )
df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))
pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_colored, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")








df_colored <- df %>%
  group_by(Source, Species) %>%
  mutate(
    rank_within_condition = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_condition, na.rm = TRUE),
    rank_color = case_when(
      rank_within_condition == 1 ~ "highest",
      rank_within_condition == max_rank ~ "lowest",
      TRUE ~ "middle"
    )
  ) %>%
  ungroup()

df_filt <- filter(df_colored, Species == "oralis" | Species =="cristatus" | Species =="gordonii" | Species == "parasanguinis" | Species =="sanguinis"| Species == "mitis" | Species == "sinensis")
df_filt <- filter(df_filt, genus == "Streptococcus")


df_filt2 <- df_filt %>%
  group_by(Species, Source) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.flipped.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.flipped.pdf")

df_filt2 <- df_filt %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()

pdf("HED_bubble_plot.hiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Condition, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Source, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.hiv.pdf")


pdf("HED_bubble_plot.healthvhiv.pdf", width = 10)
  ggplot(df_filt2, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4391BD",
      "middle" = "#AEC3DF",
      "lowest" = "#EEEAF3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.healthvhiv.pdf")

df_colored <- df %>%
  group_by(Species, Condition) %>%
  mutate(
    rank_within_source = rank(-average_baseMean, ties.method = "first"),
    max_rank = max(rank_within_source),
    rank_color = case_when(
      rank_within_source == 1 ~ "highest",          # highest
      rank_within_source == max_rank ~ "lowest",  # lowest
      TRUE ~ "middle"                                  # middle
    )
  ) %>%
  ungroup()
df_colored$rank_color <- factor(df_colored$rank_color, levels=c("highest", "middle", "lowest"))
df_strep <- df_colored %>%
  filter(grepl("^Strep", as.character(Genus)))


pdf("HED_bubble_plot.pdf", width = 10)
  ggplot(df_strep, aes(x = Source, y = Species_nested, size = average_baseMean, color = rank_color)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Species",
    size = "Average Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c(
      "highest" = "#4596C3",
      "middle" = "#4596C3",
      "lowest" = "#4596C3"
    )) +
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HED_bubble_plot.pdf")
```