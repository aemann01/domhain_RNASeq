# 1. Run DESeq HI vs HUU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/03-global-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
# submap <- submap[submap$sample_id %in% sample_list, ]
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
[1] "number of genes with adjusted p value lower than 0.05:  210857"
# out of 51006 with nonzero total read count
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 96045, 1.5%
# LFC < 0 (down)     : 114812, 1.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 5564940, 84%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  210857"# out of 51006 with nonzero total read count
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 123978, 1.9%
# LFC < 0 (down)     : 147149, 2.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 5564940, 84%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results-HIvHUU.RData")
```
Valcona Plot
```R
load("deseq_results-HIvHUU.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$tag
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
sigsp <- paste("x", sigloc$genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# remove any lists (i.e., genera) with fewer than three hits
siglist <- Filter(function(x) length(x) >=3, siglist)

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
res_ord <- res_ord %>% filter(gene != "none")
# get a large number of colors from the viridis package to iterate over
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
res_ord$genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$genus)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$tag, 10)
# negative top 10
top <- tail(sortdf$tag, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$species, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHUU.global.pdf", width=30, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.global.pdf")
```
# 2. Run DESeq HEU vs HUU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/03-global-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HEU", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  208374"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 59755, 0.91%
# LFC < 0 (down)     : 148619, 2.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 5617690, 85%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# HUU is positive, HEU negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  208374"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 90515, 1.4%
# LFC < 0 (down)     : 185339, 2.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 5617690, 85%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.table(resLFC, file="deseq_results-HEUvHUU.txt", quote=F, sep="\t")
save.image("deseq_results-HEUvHUU.RData")
```
Valcona Plot
```R
load("deseq_results-HEUvHUU.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$tag
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
sigsp <- paste("x", sigloc$genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# remove any lists (i.e., genera) with fewer than three hits
siglist <- Filter(function(x) length(x) >=3, siglist)

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
res_ord <- res_ord %>% filter(gene != "none")
# get a large number of colors from the viridis package to iterate over
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
res_ord$genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$genus)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$tag, 10)
# negative top 10
top <- tail(sortdf$tag, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$species, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HEUvHUU.global.pdf", width=30, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.global.pdf")
```
# 3. Run DESeq HI vs HEU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/03-global-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
# submap <- submap[submap$sample_id %in% sample_list, ]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  196687"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 130968, 2%
# LFC < 0 (down)     : 65719, 1%
# outliers [1]       : 0, 0%
# low counts [2]     : 5482201, 83%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# HEU is positive, HI negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  196687"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 158949, 2.4%
# LFC < 0 (down)     : 95800, 1.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 5482201, 83%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results-HIvHEU.txt", quote=F, sep="\t")
save.image("deseq_results-HIvHEU.RData")
```
Valcona Plot
```R
load("deseq_results-HIvHEU.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$tag
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
sigsp <- paste("x", sigloc$genus, sep="_")
sigdf <- as.data.frame(cbind(sigloc$tag, sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# remove any lists (i.e., genera) with fewer than three hits
siglist <- Filter(function(x) length(x) >=3, siglist)

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
res_ord <- res_ord %>% filter(gene != "none")
# get a large number of colors from the viridis package to iterate over
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
res_ord$genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$genus)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$tag, 10)
# negative top 10
top <- tail(sortdf$tag, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$species, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHEU.global.pdf", width=30, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHEU.global.pdf")
```
# 4. Proportion of aligned reads that belong to the red complex
```sh
python3 species_distro.py > species_distro.log
python3 species_distro_average.py > species_distro_average.log
```
Make plots
```R
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(ggpubr)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/03-global-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
counts <- read.csv("species_reads.txt", header=T,sep = "\t")
counts$sample <- gsub(x = counts$sample, pattern = "\\.red", replacement = "") 
counts$sample <- gsub(x = counts$sample, pattern = "\\.", replacement = "-") 
data <- melt(counts)
data <- left_join(data, metadata, by = join_by(sample ==  sample_id))
data$genus <- gsub(x = data$variable, pattern = "_.*", replacement = "")

pdf("rna_reads.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=sample, y=value,fill=genus),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./rna_reads.pdf")

# red complex plot
sub_data <- data[data$variable == "Porphyromonas_gingivalis" | data$variable == "Treponema_denticola" | data$variable == "Tannerella_forsythia",]
pdf("red_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./red_rna_reads.pdf")
mean(sub_data$value)

# Trep denticola
sub_data <- data[data$variable == "Treponema_denticola",]
pdf("tdent_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./tdent_rna_reads.pdf")
pdf("tdent_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./tdent_rna_reads.hist.pdf")
mean(sub_data$value)
# Porphyromonas_gingivalis
sub_data <- data[data$variable == "Porphyromonas_gingivalis",]
pdf("pging_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./pging_rna_reads.pdf")
pdf("pging_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pging_rna_reads.hist.pdf")
mean(sub_data$value)
# Trep denticola
sub_data <- data[data$variable == "Tannerella_forsythia",]
pdf("tfors_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./tfors_rna_reads.pdf")
pdf("tfors_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./tfors_rna_reads.hist.pdf")
mean(sub_data$value)

# compare across the thuree groups
sub_data <- data[data$variable == "Porphyromonas_gingivalis" | data$variable == "Treponema_denticola" | data$variable == "Tannerella_forsythia",]
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
pdf("red_rna_prop.pdf")
ggplot(sub_data, aes(x=factor(hiv_status, levels=hiv_stat),y=value))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./red_rna_prop.pdf")