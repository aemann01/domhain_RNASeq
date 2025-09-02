# 1. oralis
## 1.1 H HUUvHI
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
res_filt <- filter(res_ord, Genus_Species == "Streptococcus_sanguinis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-H-sanguinis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-H-sanguinis.HUUvHI.pdf")

res_filt <- filter(res_ord, Genus_Species == "Streptococcus_oralis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-H-oralis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-H-oralis.HUUvHI.pdf")
```
# 2. D HUUvHI
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
res_filt <- filter(res_ord, Genus_Species == "Streptococcus_sanguinis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-D-sanguinis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-D-sanguinis.HUUvHI.pdf")

res_filt <- filter(res_ord, Genus_Species == "Streptococcus_oralis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-D-oralis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-D-oralis.HUUvHI.pdf")
```
# 2. E HUUvHI
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
res_filt <- filter(res_ord, Genus_Species == "Streptococcus_sanguinis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-E-sanguinis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-E-sanguinis.HUUvHI.pdf")

res_filt <- filter(res_ord, Genus_Species == "Streptococcus_oralis")
res_filt$color2 <- ifelse(grepl("^arcA", as.character(res_filt$gene)), "#009E73",
                   ifelse(grepl("^arcB", as.character(res_filt$gene)), "#0072B2",
                   ifelse(grepl("^arcC", as.character(res_filt$gene)), "#E69F00",
                          "#808080")))
# if no color, remove genus label
res_filt$gene[is.na(res_filt$Genus)] <- NA
# change NA genus to grey
res_filt$color2[is.na(res_filt$Genus)] <- "#808080"

# get key value pairs for plotting
colormap2 <- setNames(res_filt$color2, res_filt$gene)

# Create volcano plot
p <- EnhancedVolcano(res_filt,
  lab = ifelse(res_filt$locus_tag %in% "bleep", paste(res_filt$Species, res_filt$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap2,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_filt) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_filt) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-E-oralis.HUUvHI.pdf", width=15, height=10)
p
dev.off()
# print to current window
system("~/.iterm2/imgcat volcano-E-oralis.HUUvHI.pdf")
```