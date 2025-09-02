# 0. ADS pathway
## 0.1 Run deseq
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(treeio)
setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
# sum togther ads genes
# Create a new column with the part before the underscore
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
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
#get baseman
# basemean_norm=rowMeans(counts(se_star, normalized=TRUE))


combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$SEQ_ID
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")


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
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
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

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HUUvHI")
head(resdf)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$SEQ_ID <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$SEQ_ID<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(resdf), as.data.frame(resHEUvHI), by = "SEQ_ID")
dfHIV <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "SEQ_ID")


#now do it for tooth health
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
# sum togther ads genes
# Create a new column with the part before the underscore
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
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
resdf <- as.data.frame(combined_df)
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$SEQ_ID
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resEvD <- lfcShrink(se_star, coef="tooth_health_E_vs_D", type="apeglm")


# get heu v huu
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E",]
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
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c( "E", "H"))

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
resHvE <- lfcShrink(se_star, coef="tooth_health_H_vs_E", type="apeglm")

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HvD")
head(resdf)

colnames(resEvD) <- paste0(colnames(resEvD), "_EvD")
resEvD$SEQ_ID <- rownames(resEvD)
head(resEvD)

colnames(resHvE) <- paste0(colnames(resHvE), "_HvE")
resHvE$SEQ_ID<- rownames(resHvE)
head(resHvE)

combined <- left_join(as.data.frame(resdf), as.data.frame(resEvD), by = "SEQ_ID")
dftooth <- left_join(as.data.frame(combined), as.data.frame(resHvE), by = "SEQ_ID")

dfHIV2 <- dfHIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
dftooth2 <- dftooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(dfHIV, dftooth, by = c("SEQ_ID","Species", "genus", "Genus", "Genus_Species"))
write.table(joined_df, file="deseq_results_ADS-all.txt", quote=F, sep="\t")
save.image("tree_data.RData")
```
## 0.2 Make heatmap for strep and leptos of interest
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_salivarius/align/salivarius.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_anginosus" | Genus_Species == "Streptococcus_australis" | Genus_Species == "Streptococcus_cristatus" | Genus_Species == "Streptococcus_constellatus" | Genus_Species == "Streptococcus_mitis" | Genus_Species == "Streptococcus_sanguinis" | Genus_Species == "Streptococcus_parasanguinis" | Genus_Species == "Streptococcus_salivarius" | Genus_Species == "Streptococcus_sinesnis" | Genus_Species == "Streptococcus_gordonii" | Genus_Species == "Streptococcus_intermedius" | Genus_Species == "Leptotrichia_sp._HMT_212" | Genus_Species == "Leptotrichia_sp._HMT_215")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_anginosus" | Genus_Species == "Streptococcus_australis" | Genus_Species == "Streptococcus_cristatus" | Genus_Species == "Streptococcus_constellatus" | Genus_Species == "Streptococcus_mitis" | Genus_Species == "Streptococcus_sanguinis" | Genus_Species == "Streptococcus_parasanguinis" | Genus_Species == "Streptococcus_salivarius" | Genus_Species == "Streptococcus_sinesnis" | Genus_Species == "Streptococcus_gordonii" | Genus_Species == "Streptococcus_intermedius" | Genus_Species == "Leptotrichia_sp._HMT_212" | Genus_Species == "Leptotrichia_sp._HMT_215")

filt_HIV2 <- filt_HIV %>% select('Genus_Species','HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('Genus_Species','H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("SEQ_ID", "Genus_Species"))

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(Genus_Species, SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .03
first_tip <- 0.367
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = .5,
         offset = 0.33,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.salivarius.pdf", height = 12, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = .5,
         offset = 0.445,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.salivarius.pdf")
```
## 0.3 Run deseq for denovo
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(treeio)
setwd("~/rna_dohmain/12-denovo_analyses")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id

# read in gene counts file
genecounts <- read.table("~/rna_dohmain/12-denovo_analyses/read_counts.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# remove empty last column (if you run this more than once it will start removing actual samples, make sure your dim after is expected)
genecounts <- genecounts[, -ncol(genecounts)]
dim(genecounts)
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
# sum togther ads genes
# Create a new column with the part before the underscore
locus2trans <- read.table("~/rna_dohmain/12-denovo_analyses/locustag_transcript.txt", header=F, sep="\t")
locus2trans$V1 <- as.character(locus2trans$V1)
locus2trans$V2 <- as.character(locus2trans$V2)
subcount_filt <- subcount[rownames(subcount) %in% locus2trans$V2, ]
subcount_filt$V1 <- locus2trans$V1[match(rownames(subcount_filt), locus2trans$V2)]
summed_subcount <- subcount_filt %>%
  group_by(V1) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()
rownames(summed_subcount) <- summed_subcount$V1
summed_subcount$V1 <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
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
#get baseman
# basemean_norm=rowMeans(counts(se_star, normalized=TRUE))


combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

# write results to file
homd <- read.csv("~/rna_dohmain/12-denovo_analyses/annotations.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
homd2 <- distinct(homd)
ann <- homd2[homd2$seqid %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$seqid
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resHEUvHI <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")


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
locus2trans <- read.table("~/rna_dohmain/12-denovo_analyses/locustag_transcript.txt", header=F, sep="\t")
locus2trans$V1 <- as.character(locus2trans$V1)
locus2trans$V2 <- as.character(locus2trans$V2)
subcount_filt <- subcount[rownames(subcount) %in% locus2trans$V2, ]
subcount_filt$V1 <- locus2trans$V1[match(rownames(subcount_filt), locus2trans$V2)]
summed_subcount <- subcount_filt %>%
  group_by(V1) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()
rownames(summed_subcount) <- summed_subcount$V1
summed_subcount$V1 <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
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

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HUUvHI")
head(resdf)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$seqid <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$seqid<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(resdf), as.data.frame(resHEUvHI), by = "seqid")
dfHIV <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "seqid")


#now do it for tooth health
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
# sum togther ads genes
# Create a new column with the part before the underscore
locus2trans <- read.table("~/rna_dohmain/12-denovo_analyses/locustag_transcript.txt", header=F, sep="\t")
locus2trans$V1 <- as.character(locus2trans$V1)
locus2trans$V2 <- as.character(locus2trans$V2)
subcount_filt <- subcount[rownames(subcount) %in% locus2trans$V2, ]
subcount_filt$V1 <- locus2trans$V1[match(rownames(subcount_filt), locus2trans$V2)]
summed_subcount <- subcount_filt %>%
  group_by(V1) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()
rownames(summed_subcount) <- summed_subcount$V1
summed_subcount$V1 <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
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
homd <- read.csv("~/rna_dohmain/12-denovo_analyses/annotations.txt", header=T, sep="\t", quote="")
resdf <- as.data.frame(combined_df)
homd2 <- distinct(homd)
ann <- homd2[homd2$seqid %in% rownames(resdf),]
# reorder
rownames(ann) <- ann$seqid
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf$genus <- resdf$Genus

# get object for bubble plot 
resEvD <- lfcShrink(se_star, coef="tooth_health_E_vs_D", type="apeglm")


# get heu v huu
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "E",]
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
locus2trans <- read.table("~/rna_dohmain/12-denovo_analyses/locustag_transcript.txt", header=F, sep="\t")
locus2trans$V1 <- as.character(locus2trans$V1)
locus2trans$V2 <- as.character(locus2trans$V2)
subcount_filt <- subcount[rownames(subcount) %in% locus2trans$V2, ]
subcount_filt$V1 <- locus2trans$V1[match(rownames(subcount_filt), locus2trans$V2)]
summed_subcount <- subcount_filt %>%
  group_by(V1) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()
rownames(summed_subcount) <- summed_subcount$V1
summed_subcount$V1 <- NULL
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~tooth_health)
# star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$tooth_health <- factor(star_results$tooth_health, levels=c( "E", "H"))

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
resHvE <- lfcShrink(se_star, coef="tooth_health_H_vs_E", type="apeglm")

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HvD")
head(resdf)

colnames(resEvD) <- paste0(colnames(resEvD), "_EvD")
resEvD$seqid <- rownames(resEvD)
head(resEvD)

colnames(resHvE) <- paste0(colnames(resHvE), "_HvE")
resHvE$seqid<- rownames(resHvE)
head(resHvE)

combined <- left_join(as.data.frame(resdf), as.data.frame(resEvD), by = "seqid")
dftooth <- left_join(as.data.frame(combined), as.data.frame(resHvE), by = "seqid")

dfHIV2 <- dfHIV %>% select('HUU', 'HEU', 'HI', 'seqid','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
dftooth2 <- dftooth %>% select('H', 'E', 'D', 'seqid','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(dfHIV, dftooth, by = c("seqid","taxonomy"))
write.table(joined_df, file="deseq_results_ADS-denovo.txt", quote=F, sep="\t")
save.image("denovo_data.RData")
```
# 1. All strep
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep && cd strep
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn  > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
mv all.fnn all.fna

# cluster using panaroo
cd ~/rna_dohmain/07-ads_expression/vince_trees/panaroo_test
mkdir gff_files
cd gff_files
awk '{print $1}' ../../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | while read line; do grep $line ../../../../homd_map/annotations.merge.txt; done > ads_annotations.tsv
grep -Ew 'Streptococcus' ads_annotations.tsv \
| awk '{print $3}' \
| sort -u \
| parallel -j 190 'cp ../../../homd_map/gff_files/{}.gff gff_files/{}.gff'
cd ../
panaroo -i ./gff_files/*gff -o results --clean-mode strict -t 190 -a core --aligner mafft -c 0.6 --core_threshold 1.00
numb_refs=$(ls gff_files/*gff | wc -l)

mkdir core_genes
grep -c ">" ./results/aligned_gene_sequences/*fas | grep ":$numb_refs" | sed 's/.\/results\/aligned_gene_sequences\///' | sed 's/:.*//' | while read line; do cp ./results/aligned_gene_sequences/$line ./core_genes/$line; done
cd core_genes
ls ./*fas | sed 's/.\/results\/aligned_gene_sequences\///' | while read line; do grep -n ">" ./$line | sed 's/;.*//' |sort |uniq | wc -l; done | grep "$numb_refs" -c
ls *fas | wc -l
# test for recombination
parallel -j 190 'Phi -o -f {} > {.}.rec' ::: *aln.fas

grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb1

grep "Max" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb2

grep "NSS" *.rec | grep permutations | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb3

cat <(awk '{print $1}' recomb1) <(awk '{print $1}' recomb2) <(awk '{print $1}' recomb3)| sed 's/rec.*/fas/' | sort | uniq -c | grep -w 3 | awk '{print $2}' | while read line; do mv $line $line.rd; done
sed -i 's/;.*//' *fas
sed -i 's/_R_//' *fas

# check alignments are all equal length
ls *fas | while read line; do samtools faidx $line; done
ls *fai | while read line; do awk '{print $2}' $line | sort | uniq | wc -l; done
# combine alignments
python3 ./combine_core.py
grep ">" core_genome.align.fa -c
samtools faidx core_genome.align.fa
# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep.core -safe

# correct annotations file
paste -d '\t' annotations.merge.txt <(awk '{print $7}' annotations.merge.txt | sed 's/Genus/Corrected/') > annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3626.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3616.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3623.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3617.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3609.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3610.1" { $8 = "Streptococcus_oralis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF3651.1" { $8 = "Streptococcus_infantis" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF4397.1" { $8 = "Streptococcus_cristatus" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF6886.1" { $8 = "Streptococcus_cristatus" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF4414.1" { $8 = "Streptococcus_gordonii" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF9611.1" { $8 = "Streptococcus_sp._HMT_056" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} $3 == "SEQF1706.1" { $8 = "Streptococcus_constellatus" } 1' annotations.updated.txt > temp && mv temp annotations.updated.txt
awk 'BEGIN{OFS=FS="\t"} {gsub(/_subsp..*/, "", $8)}1' annotations.updated.txt > temp && mv temp annotations.updated.txt

```
## 1.2 Add in the heatmpa
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)
setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/panaroo_0.6/core_genes/strep.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus == "Streptococcus")
filt_tooth <- filter(dftooth, Genus == "Streptococcus")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID', 'Genus_Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID', 'Genus_Species','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("SEQ_ID", "Genus_Species"))

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID, Genus_Species) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)
group_SEQ <- group_SEQ %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID)                  
  )

increased_by <- .005
first_tip <- 1.3
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(aes(label=species2), align = TRUE, size = 3.5, offset = 0.0) +
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 7, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 7, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.01,
         offset = 0.069,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 3,
         color = "black")

pdf("tree.strep_ref.pdf", height = 135, width =80)
  gheatmap(p1,
         group_SEQ3,
         width = 0.01,
         offset = 0.0825,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 3,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.strep_ref.pdf")


p1 <- gheatmap(p,
         log10(group_SEQ2),
         width = 0.01,
         offset = 0.069,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 3,
         color = "black")

pdf("tree.strep_ref.log10.pdf", height = 135, width =80)
  gheatmap(p1,
         log10(group_SEQ3),
         width = 0.01,
         offset = 0.0825,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 3,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.strep_ref.log10.pdf")
```
## 1.3 Create rpoC strep tree
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_rpoC && cd strep_rpoC
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep "RNA polymerase subunit beta'" -A1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa
fasttree -nt c0.align.fa > rpoC.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep.core -safe
```
# 2. Strep oralis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_oralis && cd strep_oralis
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'oralis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

# alt program to test for recomdinariotion
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix oralis.core -safe
```
## 2.1 Color tree based on ads genes
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HUUvHI")
head(resdf)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$locus_tag <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$locus_tag<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(resdf), as.data.frame(resHEUvHI), by = "locus_tag")
dfHIV <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "locus_tag")
dffilt <- filter(dfHIV, Genus_Species == "Streptococcus_oralis")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_oralis/align/oralis.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
dffilt2 <- dffilt %>% select('HUU', 'HEU', 'HI', 'genus', 'Species', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_hiv <- dffilt2 %>%
  group_by(SEQ_ID) %>%
  summarise(
    # Direction calls
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
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE)
  )
group_hiv <- as.data.frame(group_hiv)
# row.names(group_hiv) <- group_hiv$SEQ_ID
# groupHI<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HI','up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHI$condition <- 'HI'
# colnames(groupHI)[2] <- "average_baseMean"
# groupHEU<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HEU', 'up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHEU$condition <- 'HEU'
# colnames(groupHEU)[2] <- "average_baseMean"
# groupHUU<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HUU', 'up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHUU$condition <- 'HUU'
# colnames(groupHUU)[2] <- "average_baseMean"

# merged <- rbind(groupHI, groupHEU, groupHUU)
# merged$all_HIvHUU <- ifelse(merged$down_HUUvHI >2, "all_HI", ifelse(merged$up_HUUvHI >2, "all_HUU", NA))
# merged$all_HEUvHUU <- ifelse(merged$down_HUUvHEU >2, "all_HEU", ifelse(merged$up_HUUvHEU >2, "all_HUU", NA))
# merged$all_HEUvHI <- ifelse(merged$down_HEUvHI >2, "all_HI", ifelse(merged$up_HEUvHI >2, "all_HEU", NA))


group_hiv$all_HIvHUU <- ifelse(group_hiv$down_HUUvHI >2, "all_HI", ifelse(group_hiv$up_HUUvHI >2, "all_HUU", "not_significant"))
group_hiv$all_HEUvHUU <- ifelse(group_hiv$down_HUUvHEU >2, "all_HEU", ifelse(group_hiv$up_HUUvHEU >2, "all_HUU", "not_significant"))
group_hiv$all_HEUvHI <- ifelse(group_hiv$down_HEUvHI >2, "all_HI", ifelse(group_hiv$up_HEUvHI >2, "all_HEU", "not_significant"))

all((group_hiv$down_HUUvHI > 0) + (group_hiv$up_HUUvHI > 0) <= 1)
group_hiv$one_HIvHUU <- ifelse(group_hiv$down_HUUvHI >0, "HI", ifelse(group_hiv$up_HUUvHI >0, "HUU", "not_significant"))
all((group_hiv$down_HUUvHEU > 0) + (group_hiv$up_HUUvHEU > 0) <= 1)
group_hiv$one_HEUvHUU <- ifelse(group_hiv$down_HUUvHEU >0, "HEU", ifelse(group_hiv$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_hiv$down_HEUvHI > 0) + (group_hiv$up_HEUvHI > 0) <= 1)
group_hiv$one_HEUvHI <- ifelse(group_hiv$down_HEUvHI >0, "HI", ifelse(group_hiv$up_HEUvHI >0, "HEU", "not_significant"))


group_hiv2 <- group_hiv %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_hiv2) <- group_hiv2$SEQ_ID
group_hiv2 <- rename(group_hiv2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )

group_hiv2$SEQ_ID <- NULL
pruned_tree <- keep.tip(tree.root, group_hiv$SEQ_ID)

p <- ggtree(pruned_tree) %<+% group_hiv +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = all_HIvHUU, x = 0.061), size = 4, shape = 16) + 
  scale_color_manual(
    name = "HI vs HUU",
    values = c("all_HI" = "#8213A0", "all_HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = all_HEUvHUU, x = 0.063), size = 4, shape = 16) + 
  scale_color_manual(
    name = "HEU vs HUU",
    values = c("all_HEU" = "#FA78FA", "all_HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = all_HEUvHI, x = 0.065), size = 4, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("all_HEU" = "#FA78FA", "all_HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = 0.067), size = 4, shape = 15) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHUU, x = 0.069), size = 4, shape = 15) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = 0.071), size = 4, shape = 15) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

pdf("tree.oralis.pdf", height = 20, width =15)
gheatmap(p,
         group_hiv2,
         width = 0.5,
         offset = 0.019,
         colnames_position = "top",
         # colnames_angle = 90,
         font.size = 8,
         color = "black")+
	scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.oralis.pdf")
```
## 2.2 Color tree based on ads pathway and tooth health
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_oralis/align/oralis.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_oralis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_oralis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .0075
first_tip <- 0.136
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 8, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 8, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 8, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 8, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 8, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 8, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.25,
         offset = 0.072,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 5,
         color = "black")

pdf("tree.oralis.pdf", height = 15, width =12)
  gheatmap(p1,
         group_SEQ3,
         width = 0.25,
         offset = 0.1,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.oralis.pdf")
```
## 2.3 See sample distro for specific lineage (SEQF6415.1)
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
# sum togther ads genes
# Create a new column with the part before the underscore
subcount$SEQF_group <- sub("^([^_]+)_.*", "\\1", rownames(subcount))
summed_subcount <- subcount %>%
  group_by(SEQF_group) %>%
  summarise(across(where(is.numeric), sum))
summed_subcount <- as.data.frame(summed_subcount)
row.names(summed_subcount) <- summed_subcount$SEQF_group
summed_subcount$SEQF_group <- NULL

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = summed_subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU", "HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
#get relative abundance of dna
load("../../rpoc/ps.RData")
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Streptococcus_oralis")
data <- psmelt(actino) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, visit_num, tooth_health, Abundance, V8, OTU)
red_dna$nucl <- "dna"
red_dna <- red_dna %>%
  rename(sample = Sample)
red_dna <- red_dna %>%
  rename(species = V8)
red_dna <- red_dna %>%
  rename(value = Abundance)
red_dna <- red_dna %>%
  rename(ASV = OTU)
red_dna[red_dna$sample == "DM00008V1PQ16-2", ]
#get relative abundance of rna
# red_rna <- as.data.frame(counts(se_star, normalized = TRUE)["SEQF6415.1", ])
# write results to file
homd <- read.csv("~/rna_dohmain/homd_map/annotations.updated.txt", header=T, sep="\t", quote="")
strep_orals_ids <- unique((filter(homd, Corrected_Species == "Streptococcus_oralis")$SEQ_ID))
norm_counts <- as.data.frame(counts(se_star, normalized = TRUE))
strep_norm_counts_long <- norm_counts[rownames(norm_counts) %in% strep_orals_ids, ]
strep_norm_counts_long$SEQ_ID <- rownames(strep_norm_counts_long)
strep_norm_counts_long <- melt(strep_norm_counts_long, id.vars = "SEQ_ID", variable.name = "sample", value.name = "value")

red_rna <- strep_norm_counts_long %>%
  mutate(SEQ_ID = ifelse(SEQ_ID == "SEQF6415.1", SEQ_ID, "other_oralis")) %>%
  group_by(sample, SEQ_ID) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  arrange(sample, SEQ_ID)
red_rna$nucl <- "rna"

map$sample <- row.names(map)
meta <- as.data.frame(as.matrix(map)) %>% dplyr::select(sample, hiv_status, tooth_health, visit_num)
red_rna <- left_join(meta, red_rna, by = "sample")
red_rna$species <- "Streptococcus_oralis"
red_rna$ASV <- "NaN"
red_rna <- red_rna[, c("sample", "hiv_status", "visit_num", "tooth_health","value", "species", "ASV" "nucl")]

#combine
combined_df <- rbind(red_rna, red_dna)
# plot
combined_df$hiv_status <- factor(combined_df$hiv_status, levels = c("HUU", "HEU", "HI"))
combined_df$tooth_health <- factor(combined_df$tooth_health, levels = c("H", "E", "D"))

rna_df <- combined_df[combined_df$nucl == "rna", ]
dna_df <- combined_df[combined_df$nucl == "dna", ]

rna_df <- red_rna %>%
  mutate(hiv_status = factor(hiv_status, levels = c("HUU", "HEU", "HI"))) %>%
  arrange(hiv_status, desc(value)) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

pdf("SEQF6415.RNA.hiv.pdf", width =10, height =10)
ggplot(rna_df, aes(x = sample, y = value, fill = SEQ_ID)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1)
  ) +
  labs(
    title = "Streptococcus oralis RNA: SEQF6415.1 vs Other IDs (by HIV status)",
    x = "Sample",
    y = "Normalized Counts",
    fill = "SEQ_ID"
  ) +
  scale_y_continuous(labels = scales::comma)
dev.off()
system("~/.iterm2/imgcat ./SEQF6415.RNA.hiv.pdf")

rna_df2 <- rna_df %>%
  mutate(tooth_health = factor(tooth_health, levels = c("H", "E", "D"))) %>%
  arrange(tooth_health, desc(value)) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

pdf("SEQF6415.RNA.tooth.pdf", width =8, height =10)
ggplot(rna_df2, aes(x = sample, y = value, fill = SEQ_ID)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1)
  ) +
  labs(
    title = "Streptococcus oralis RNA: SEQF6415.1 vs Other IDs (by Tooth Health)",
    x = "Sample",
    y = "Normalized Counts",
    fill = "SEQ_ID"
  ) +
  scale_y_continuous(labels = scales::comma)
dev.off()
system("~/.iterm2/imgcat ./SEQF6415.RNA.tooth.pdf")
```
## 2.4 Make tree including denovo: tree based on ADS
```sh
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_oralis_denovo && cd strep_oralis_denovo
# get reference ADS
grep Streptococcus -A 1 ../../ads_operons.fna | grep -w 'oralis' -A 1 > ads_operons.filt.fna
sed -i 's/|.*//' ads_operons.filt.fna
sed -i '/^--$/d' ads_operons.filt.fna
# get denovo ADS
grep Streptococcus_oralis ~/rna_dohmain/12-denovo_analyses/annotations.txt | awk '{print $1}' | sed 's/_R_//' | while read line; do grep -A 1 $line  ~/rna_dohmain/12-denovo_analyses/arcABC_operons.fna; done > ads_denovo.filt.fna
sed -i 's/|.*//' ads_denovo.filt.fna
sed -i '/^--$/d' ads_denovo.filt.fna

cat ads_denovo.filt.fna ads_operons.filt.fna > ads_operans.comb.fna
mafft  --thread -1 --adjustdirectionaccurately \
  ads_operans.comb.fna > arcABC_operons.align.fna

trimal -in arcABC_operons.align.fna \
  -out arcABC_operons.trim.fna \
  -htmlout arcABC_operons.trim.html \
  -gt 0.5 \
  -resoverlap 0.5 \
  -seqoverlap 50

iqtree -s arcABC_operons.trim.fna -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix oralis.denovo -safe
```
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_oralis_denovo/oralis.denovo.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tree.root$tip.label <- gsub("^_R_", "", tree.root$tip.label)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_oralis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_oralis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

#now repeat for denovo data
load("~/rna_dohmain/12-denovo_analyses/denovo_data.RData")

filt_HIV <- dfHIV[grepl("Streptococcus_oralis", dfHIV$taxonomy), ]
filt_tooth <- dftooth[grepl("Streptococcus_oralis", dftooth$taxonomy), ]

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'seqid','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'seqid','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "seqid")
joined_df <- joined_df %>%
  rename(SEQ_ID = seqid)
group_SEQ_2 <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ_2 <- as.data.frame(group_SEQ_2)
group_SEQ_2$SEQ_ID <- gsub("^_R_", "",group_SEQ_2$SEQ_ID)

all((group_SEQ_2$down_HUUvHI > 0) + (group_SEQ_2$up_HUUvHI > 0) <= 1)
group_SEQ_2$one_HIvHUU <- ifelse(group_SEQ_2$down_HUUvHI >0, "HI", ifelse(group_SEQ_2$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HUUvHEU > 0) + (group_SEQ_2$up_HUUvHEU > 0) <= 1)
group_SEQ_2$one_HEUvHUU <- ifelse(group_SEQ_2$down_HUUvHEU >0, "HEU", ifelse(group_SEQ_2$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HEUvHI > 0) + (group_SEQ_2$up_HEUvHI > 0) <= 1)
group_SEQ_2$one_HEUvHI <- ifelse(group_SEQ_2$down_HEUvHI >0, "HI", ifelse(group_SEQ_2$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ_2$down_HvD > 0) + (group_SEQ_2$up_HvD > 0) <= 1)
group_SEQ_2$one_HvD <- ifelse(group_SEQ_2$down_HvD >0, "D", ifelse(group_SEQ_2$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ_2$down_HvE > 0) + (group_SEQ_2$up_HvE > 0) <= 1)
group_SEQ_2$one_HvE <- ifelse(group_SEQ_2$down_HvE >0, "E", ifelse(group_SEQ_2$up_HvE >0, "H", "not_significant"))
all((group_SEQ_2$down_EvD > 0) + (group_SEQ_2$up_EvD > 0) <= 1)
group_SEQ_2$one_EvD <- ifelse(group_SEQ_2$down_EvD >0, "D", ifelse(group_SEQ_2$up_EvD >0, "E", "not_significant"))

group_SEQ_21 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ_21) <- group_SEQ_21$SEQ_ID
group_SEQ_21 <- rename(group_SEQ_21,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ_21$SEQ_ID <- NULL

group_SEQ_22 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ_22) <- group_SEQ_22$SEQ_ID
group_SEQ_22 <- rename(group_SEQ_22,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ_22$SEQ_ID <- NULL


# combine
group1 <- rbind(group_SEQ_21, group_SEQ2)
group2 <- rbind(group_SEQ_22, group_SEQ3)
group_overall <- rbind(group_SEQ_2, group_SEQ)

pruned_tree <- drop.tip(tree.root, "SEQF1998.1")


increased_by <- .008
first_tip <- .303
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_overall +
  geom_tiplab(align = TRUE, size = 3.5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group1,
         width = 0.1,
         offset = 0.14,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.7,
         font.size = 5,
         color = "black")

pdf("tree.denovo_oralis.pruned.pdf", height = 18, width =20)
  gheatmap(p1,
         group2,
         width = 0.1,
         offset = 0.14,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.denovo_oralis.pruned.pdf")

p1 <- gheatmap(p,
         log10(group1),
         width = 0.1,
         offset = 0.405,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.7,
         font.size = 5,
         color = "black")

pdf("tree.denovo_oralis.pruned_log.pdf", height = 18, width =20)
  gheatmap(p1,
         log10(group2),
         width = 0.1,
         offset = 0.475,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "log10(Average basemean)")
dev.off()
system("~/.iterm2/imgcat tree.denovo_oralis.pruned_log.pdf")

# group2[row.names(group2) == "SEQF6415.1", ]
```
# 3. Strep sanguinis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_sang && cd strep_sang
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ../ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'sanguinis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab



# # try using mcl blast
gffread ALL_genomes.gff -g ALL_genomes.fna -y proteins.faa
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/proteins.faa > ../ALL_genomes.faa
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/.[123456789]_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.faa' > ads.faa
grep Streptococcus_sanguinis ~/rna_dohmain/homd_map/annotations.merge.txt | awk '{print $3}' | sort | uniq | parallel -j 100 'grep -A 1 {} ./ads.faa' > all.faa
# format fasta headers 
grep -v "^>" ads.faa > seqs
grep "^>" ads.faa | sed 's/ .*//' | sed 's/>//' | awk '{print $0 "|" $0}' | sed 's/_.*|/|/' | sed 's/^/>/' > headers
paste -d '\n' headers seqs > homd_for_mcl.faa
# make blast database
makeblastdb -in homd_for_mcl.faa -dbtype prot -out homd_for_mcl.blast.db
# now run blast 
cat homd_for_mcl.faa | parallel --block 100k --recstart '>' --pipe blastp -db homd_for_mcl.blast.db -evalue 1e-5 -outfmt 6 -query - > blast.out





# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix sang.core -safe
```
## 3.1 Color tree
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
# reload data to filter samples of interest
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter by HIV group
# submap <- metadata[metadata$tooth_health == "H",]
# filter out enamel cavity
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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

colnames(resdf)[1:5] <- paste0(colnames(resdf)[1:5], "_HUUvHI")
head(resdf)
colnames(resHEUvHI) <- paste0(colnames(resHEUvHI), "_HEUvHI")
resHEUvHI$locus_tag <- rownames(resHEUvHI)
head(resHEUvHI)
colnames(resHUUvHEU) <- paste0(colnames(resHUUvHEU), "_HUUvHEU")
resHUUvHEU$locus_tag<- rownames(resHUUvHEU)
head(resHUUvHEU)

combined <- left_join(as.data.frame(resdf), as.data.frame(resHEUvHI), by = "locus_tag")
dfHIV <- left_join(as.data.frame(combined), as.data.frame(resHUUvHEU), by = "locus_tag")
dffilt <- filter(dfHIV, Genus_Species == "Streptococcus_sanguinis")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_sang/align/sang.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
dffilt2 <- dffilt %>% select('HUU', 'HEU', 'HI', 'genus', 'Species', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')

# group and get average base mean
group_hiv <- dffilt2 %>%
  group_by(SEQ_ID) %>%
  summarise(
    # Direction calls
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
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE)
  )
group_hiv <- as.data.frame(group_hiv)
# row.names(group_hiv) <- group_hiv$SEQ_ID
# groupHI<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HI','up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHI$condition <- 'HI'
# colnames(groupHI)[2] <- "average_baseMean"
# groupHEU<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HEU', 'up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHEU$condition <- 'HEU'
# colnames(groupHEU)[2] <- "average_baseMean"
# groupHUU<- select(group_hiv, c('SEQ_ID', 'average_baseMean_HUU', 'up_HUUvHI', 'down_HUUvHI', 'up_HUUvHEU', 'down_HUUvHEU', 'up_HEUvHI','down_HEUvHI')) 
# groupHUU$condition <- 'HUU'
# colnames(groupHUU)[2] <- "average_baseMean"

# merged <- rbind(groupHI, groupHEU, groupHUU)
# merged$all_HIvHUU <- ifelse(merged$down_HUUvHI >2, "all_HI", ifelse(merged$up_HUUvHI >2, "all_HUU", NA))
# merged$all_HEUvHUU <- ifelse(merged$down_HUUvHEU >2, "all_HEU", ifelse(merged$up_HUUvHEU >2, "all_HUU", NA))
# merged$all_HEUvHI <- ifelse(merged$down_HEUvHI >2, "all_HI", ifelse(merged$up_HEUvHI >2, "all_HEU", NA))


group_hiv$all_HIvHUU <- ifelse(group_hiv$down_HUUvHI >2, "all_HI", ifelse(group_hiv$up_HUUvHI >2, "all_HUU", "not_significant"))
group_hiv$all_HEUvHUU <- ifelse(group_hiv$down_HUUvHEU >2, "all_HEU", ifelse(group_hiv$up_HUUvHEU >2, "all_HUU", "not_significant"))
group_hiv$all_HEUvHI <- ifelse(group_hiv$down_HEUvHI >2, "all_HI", ifelse(group_hiv$up_HEUvHI >2, "all_HEU", "not_significant"))


all((group_hiv$down_HUUvHI > 0) + (group_hiv$up_HUUvHI > 0) <= 1)
group_hiv <- group_hiv %>%
  mutate(one_HIvHUU = case_when(
    down_HUUvHI > 0 & up_HUUvHI > 0 ~ "both",
    down_HUUvHI > 0 ~ "HI",
    up_HUUvHI > 0 ~ "HUU",
    TRUE ~ "not_significant"
  ))
all((group_hiv$down_HUUvHEU > 0) + (group_hiv$up_HUUvHEU > 0) <= 1)
group_hiv <- group_hiv %>%
  mutate(one_HEUvHUU = case_when(
    down_HUUvHEU > 0 & up_HUUvHEU > 0 ~ "both",
    down_HUUvHEU > 0 ~ "HI",
    up_HUUvHI > 0 ~ "HUU",
    TRUE ~ "not_significant"
  ))
all((group_hiv$down_HEUvHI > 0) + (group_hiv$up_HEUvHI > 0) <= 1)
group_hiv <- group_hiv %>%
  mutate(one_HEUvHI = case_when(
    down_HEUvHI > 0 & up_HEUvHI > 0 ~ "both",
    down_HEUvHI > 0 ~ "HI",
    up_HUUvHI > 0 ~ "HUU",
    TRUE ~ "not_significant"
  ))


group_hiv2 <- group_hiv %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_hiv2) <- group_hiv2$SEQ_ID
group_hiv2 <- rename(group_hiv2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )

group_hiv2$SEQ_ID <- NULL
pruned_tree <- keep.tip(tree.root, group_hiv$SEQ_ID)

p <- ggtree(pruned_tree) %<+% group_hiv +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HIvHUU, x = 0.065), size = 4, shape = 15) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F", "both" = "green")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHUU, x = 0.067), size = 4, shape = 15) +
  scale_color_manual(
    name = "HEU vs HUU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F", "both" = "green")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = 0.069), size = 4, shape = 15) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F", "both" = "green")
  ) +
  theme(legend.position = "top")

pdf("tree.sang.pdf", height = 20, width =15)
gheatmap(p,
         group_hiv2,
         width = 0.5,
         offset = 0.015,
         colnames_position = "top",
         # colnames_angle = 90,
         font.size = 8,
         color = "black")+
	scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.sang.pdf")
```
## 3.2 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_sang/align/sang.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_sanguinis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_sanguinis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .0018
first_tip <- 0.0495
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.2,
         offset = 0.0203,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.8,
         font.size = 5,
         color = "black")

pdf("tree.sang.pdf", height = 21, width =12)
  gheatmap(p1,
         group_SEQ3,
         width = 0.2,
         offset = 0.0282,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.sang.pdf")
# dffilt2[dffilt2$SEQ_ID == "SEQF5728.1", ]
```
# 4. Strep parasanguinis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_para && cd strep_para
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'parasanguinis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.01)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix para.core -safe
```
## 4.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_para/align/para.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_parasanguinis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_parasanguinis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .0018
first_tip <- 0.0465
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.2,
         offset = 0.0203,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.8,
         font.size = 5,
         color = "black")

pdf("tree.para.pdf", height = 21, width =12)
  gheatmap(p1,
         group_SEQ3,
         width = 0.2,
         offset = 0.0282,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.para.pdf")
```
# 5. Strep gordonii
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_gordonii && cd strep_gordonii
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'gordonii' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix gordonii.core -safe
```
## 5.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_gordonii/align/gordonii.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_gordonii")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_gordonii")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .025
first_tip <- 0.85
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.1,
         offset = 0.294,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.8,
         font.size = 5,
         color = "black")
p2 <- p1 + new_scale_fill()

pdf("tree.gordonii.pdf", height = 21, width =12)
  gheatmap(p1,
         group_SEQ3,
         width = 0.1,
         offset = 0.37,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.gordonii.pdf")
```
# 6. Strep cristatus
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_cristatus && cd strep_cristatus
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'cristatus' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix cristatus.core -safe
```
## 6.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_cristatus/align/cristatus.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_cristatus")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_cristatus")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .02
first_tip <- 0.513
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.18,
         offset = 0.22,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")
p2 <- p1 + new_scale_fill()

pdf("tree.cristatus.pdf", height = 14, width =12)
  gheatmap(p1,
         group_SEQ3,
         width = 0.18,
         offset = 0.3,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.cristatus.pdf")
```
# 7. Strep mitis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_mitis && cd strep_mitis
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'mitis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix mitis.core -safe
```
## 7.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_mitis/align/mitis.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_mitis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_mitis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .035
first_tip <- 0.513
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.4,
         offset = 0.42,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.1,
         font.size = 5,
         color = "black")

pdf("tree.mitis.pdf", height = 5, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = 0.4,
         offset = 0.55,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.mitis.pdf")
```
# 8. Streptococcus anginosus
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_anginosus && cd strep_anginosus
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'anginosus' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix anginosus.core -safe
```
## 8.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_anginosus/align/anginosus.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_anginosus")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_anginosus")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .05
first_tip <- 0.76
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 0.38,
         offset = 0.585,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.anginosus.pdf", height = 15, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = 0.38,
         offset = 0.775,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.anginosus.pdf")
```
# 9. Streptococcus australis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_australis && cd strep_australis
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'australis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix australis.core -safe
```
## 9.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_australis/align/australis.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_australis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_australis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .03
first_tip <- 0.233
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 2.5,
         offset = 0.32,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.australis.pdf", height = 12, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = 2.5,
         offset = 0.485,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.australis.pdf")
```
# 10. Streptococcus constellatus
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_constellatus && cd strep_constellatus
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'constellatus' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix constellatus.core -safe
```
## 10.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_constellatus/align/constellatus.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_constellatus")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_constellatus")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .015
first_tip <- 0.1
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 4,
         offset = 0.155,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.constellatus.pdf", height = 12, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = 4,
         offset = 0.24,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.constellatus.pdf")
```
# 11. Streptococcus intermedius
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_intermedius && cd strep_intermedius
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'intermedius' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.01)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix intermedius.core -safe
```
## 11.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_intermedius/align/intermedius.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_intermedius")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_intermedius")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .03
first_tip <- 0.367
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(pruned_tree) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = .5,
         offset = 0.33,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.intermedius.pdf", height = 12, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = .5,
         offset = 0.445,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.intermedius.pdf")
```
# 12. Streptococcus salivarius
only one species had annotated ADS
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_salivarius && cd strep_salivarius
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'salivarius' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
# ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
# grep "Normal" *.rec | while IFS=: read -r file var pval; do
#     # Remove leading/trailing whitespace
#     pval=$(echo "$pval" | xargs)

#     # Check for missing values
#     if [[ "$pval" == "--" ]]; then
#         echo "$file:$var:        $pval  # Not significant (missing)"
#     else
#         # Convert scientific notation to float and compare with 0.05
#         if awk "BEGIN {exit !($pval < 0.01)}"; then
#             echo "$file:$var:        $pval  # Not significant (p > 0.05)"
#         fi
#     fi
# done > recomb
# awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix salivarius.core -safe
```
## 12.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_salivarius/align/salivarius.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_salivarius")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_salivarius")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

# pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .027
first_tip <- 0.3
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(tree.root) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = .5,
         offset = 0.29,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.salivarius.pdf", height = 2, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = .5,
         offset = 0.38,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.salivarius.pdf")
```
# 13. Streptococcus sinensis
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_sinensis && cd strep_sinensis
# awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
awk '{print $1}' ../../arcGene_read_counts.cleaned.txt | sed 's/_.*//' | sort | uniq | parallel -j 100 'grep -A 1 {} ../ALL_genomes.ffn' > ads.fnn
grep Streptococcus -A 1 ads.fnn | grep -w 'sinensis' -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ~/rna_dohmain/11-perio/06-phylogenies/single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.001)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb
awk '{print $1}' recomb | sed 's/rec.*/align.fa/' | while read line; do mv $line $line.rd; done

## alt recombination method
# ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
# grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
# cat recomb | while read line; do mv $line $line.rd; done

# combine single copy core
# combine alignments
python3 ~/rna_dohmain/11-perio/06-phylogenies/combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix sinensis.core -safe
```
also only has one genome with ads activity
## 13.1 Color tree based on ads pathway
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_sinensis/align/core_genome.fast.tre")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Streptococcus_sinensis")
filt_tooth <- filter(dftooth, Genus_Species == "Streptococcus_sinensis")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

# pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .009
first_tip <- 0.065
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(tree.root) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = 2,
         offset = 0.095,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.sinensis.pdf", height = 2, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = 2,
         offset = 0.14,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.sinensis.pdf")
```
# 14.1 Leptos
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/lepto_212/align/lepto.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus_Species == "Leptotrichia_sp._HMT_212" | Genus_Species == "Leptotrichia_sp._HMT_215")
filt_tooth <- filter(dftooth, Genus_Species == "Leptotrichia_sp._HMT_212" | Genus_Species == "Leptotrichia_sp._HMT_215" )

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = "SEQ_ID")

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

# pruned_tree <- keep.tip(tree.root, group_SEQ$SEQ_ID)

increased_by <- .027
first_tip <- 0.3
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(tree.root) %<+% group_SEQ +
  geom_tiplab(align = TRUE, size = 5, offset = 0.0009) + 
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ2,
         width = .5,
         offset = 0.29,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.3,
         font.size = 5,
         color = "black")

pdf("tree.lepto.pdf", height = 2, width =8)
  gheatmap(p1,
         group_SEQ3,
         width = .5,
         offset = 0.38,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.lepto.pdf")
```
# 15. Denovo ADS tree for strep
```sh
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir denovo_strep && cd denovo_strep

grep -E 'Streptococcus' ~/rna_dohmain/12-denovo_analyses/annotations.txt | awk '{print $1}' | sed 's/_R_//' | while read line; do grep -A 1 $line  ~/rna_dohmain/12-denovo_analyses/arcABC_operons.fna; done > ads_denovo.filt.fna
sed -i 's/|.*//' ads_denovo.filt.fna
sed -i '/^--$/d' ads_denovo.filt.fna
grep -c ">" ads_denovo.filt.fna
mafft  --thread -1 --adjustdirectionaccurately \
  ads_denovo.filt.fna  > arcABC_operons.align.fna

trimal -in arcABC_operons.align.fna \
  -out arcABC_operons.trim.fna \
  -htmlout arcABC_operons.trim.html \
  -gt 0.5 \
  -resoverlap 0.5 \
  -seqoverlap 50

iqtree -s arcABC_operons.trim.fna -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep.denovo -safe
```
## 15.1 Make tree
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/denovo_strep/strep.denovo.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tree.root$tip.label <- gsub("^_R_", "", tree.root$tip.label)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#now repeat for denovo data
load("~/rna_dohmain/12-denovo_analyses/denovo_data.RData")

filt_HIV <- dfHIV[grepl("Streptococcu", dfHIV$taxonomy), ]
filt_tooth <- dftooth[grepl("Streptococcus", dftooth$taxonomy), ]

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'seqid','padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI', 'taxonomy')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'seqid','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD', 'taxonomy')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("seqid", "taxonomy"))
joined_df <- joined_df %>%
  rename(SEQ_ID = seqid)
group_SEQ_2 <- joined_df %>%
  group_by(SEQ_ID, taxonomy) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ_2 <- as.data.frame(group_SEQ_2)
group_SEQ_2$SEQ_ID <- gsub("^_R_", "",group_SEQ_2$SEQ_ID)

all((group_SEQ_2$down_HUUvHI > 0) + (group_SEQ_2$up_HUUvHI > 0) <= 1)
group_SEQ_2$one_HIvHUU <- ifelse(group_SEQ_2$down_HUUvHI >0, "HI", ifelse(group_SEQ_2$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HUUvHEU > 0) + (group_SEQ_2$up_HUUvHEU > 0) <= 1)
group_SEQ_2$one_HEUvHUU <- ifelse(group_SEQ_2$down_HUUvHEU >0, "HEU", ifelse(group_SEQ_2$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HEUvHI > 0) + (group_SEQ_2$up_HEUvHI > 0) <= 1)
group_SEQ_2$one_HEUvHI <- ifelse(group_SEQ_2$down_HEUvHI >0, "HI", ifelse(group_SEQ_2$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ_2$down_HvD > 0) + (group_SEQ_2$up_HvD > 0) <= 1)
group_SEQ_2$one_HvD <- ifelse(group_SEQ_2$down_HvD >0, "D", ifelse(group_SEQ_2$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ_2$down_HvE > 0) + (group_SEQ_2$up_HvE > 0) <= 1)
group_SEQ_2$one_HvE <- ifelse(group_SEQ_2$down_HvE >0, "E", ifelse(group_SEQ_2$up_HvE >0, "H", "not_significant"))
all((group_SEQ_2$down_EvD > 0) + (group_SEQ_2$up_EvD > 0) <= 1)
group_SEQ_2$one_EvD <- ifelse(group_SEQ_2$down_EvD >0, "D", ifelse(group_SEQ_2$up_EvD >0, "E", "not_significant"))

group_SEQ_21 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ_21) <- group_SEQ_21$SEQ_ID
group_SEQ_21 <- rename(group_SEQ_21,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ_21$SEQ_ID <- NULL

group_SEQ_22 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ_22) <- group_SEQ_22$SEQ_ID
group_SEQ_22 <- rename(group_SEQ_22,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ_22$SEQ_ID <- NULL


# pruned_tree <- drop.tip(tree.root, "SEQF1998.1")
group_SEQ_2$clean_taxonomy <- sapply(group_SEQ_2$taxonomy, function(tax) {
  # If it starts with "Streptococcus_sp.", leave it unchanged
  if (grepl("^Streptococcus_sp\\.", tax)) {
    tax
  } else {
    # Otherwise split by "_" and keep first two parts if possible
    parts <- strsplit(tax, "_")[[1]]
    if (length(parts) >= 2) {
      paste(parts[1:2], collapse = "_")
    } else {
      parts[1]
    }
  }
})

increased_by <- .01
first_tip <- .665
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by



p <- ggtree(tree.root) %<+% group_SEQ_2 +
  geom_tiplab(aes(color = clean_taxonomy, stroke =2), align = TRUE, size = 3.5, offset = 0.0009) +
  scale_color_manual(
    name = "Taxonomy",
    values = c(
      "Streptococcus_anginosus" = "#000000",
      "Streptococcus_sanguinis" = "#E69F00",
      "Streptococcus_parasanguinis" = "#56B4E9",
      "Streptococcus_gordonii" = "#009E73",
      "Streptococcus_oralis" = "#8C5109",
      "Streptococcus_sp._A12" = "#0072B2",
      "Streptococcus_cristatus" = "#D55E00",
      "Streptococcus_pneumoniae" = "#CC79A7",
      "Streptococcus_sp._NPS_308" = "#999999"
    )
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(legend.position = "top")

p1 <- gheatmap(p,
         group_SEQ_21,
         width = 0.07,
         offset = 0.143,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.7,
         font.size = 5,
         color = "black")

pdf("tree.denovo_strep.pdf", height = 26, width =20)
  gheatmap(p1,
         group_SEQ_22,
         width = 0.07,
         offset = 0.19,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.denovo_strep.pdf")

# creat output file for Vince
summary_df <- group_SEQ_2 %>%
  group_by(clean_taxonomy, SEQ_ID) %>%
  summarise(
    count_HUU = sum(one_HEUvHUU == "HUU", na.rm = TRUE) + sum(one_HIvHUU == "HUU", na.rm = TRUE),
    count_HEU = sum(one_HEUvHUU == "HEU", na.rm = TRUE) + sum(one_HEUvHI == "HEU", na.rm = TRUE),
    count_HI = sum(one_HIvHUU == "HI", na.rm = TRUE) + sum(one_HEUvHI == "HI", na.rm = TRUE),
    count_H = sum(one_HvE == "H", na.rm = TRUE) + sum(one_HvD == "H", na.rm = TRUE),
    count_E = sum(one_HvE == "E", na.rm = TRUE) + sum(one_EvD == "E", na.rm = TRUE),
    count_D = sum(one_HvD == "D", na.rm = TRUE) + sum(one_EvD == "D", na.rm = TRUE)
    ) %>%
  group_by(clean_taxonomy) %>%
  summarise(
    number_of_lineages = n_distinct(SEQ_ID),
    HUU = sum(count_HUU),
    HEU  = sum(count_HEU),
    HI  = sum(count_HI),
    H     = sum(count_H),
    E     = sum(count_E),
    D     = sum(count_D),
    .groups = "drop"
  )
summary_df <- summary_df %>%
  bind_rows(
    summary_df %>%
      summarise(
        clean_taxonomy = "all",
        number_of_lineages = sum(number_of_lineages, na.rm = TRUE),
        HUU = sum(HUU, na.rm = TRUE),
        HEU = sum(HEU, na.rm = TRUE),
        HI = sum(HI, na.rm = TRUE),
        H = sum(H, na.rm = TRUE),
        E = sum(E, na.rm = TRUE),
        D = sum(D, na.rm = TRUE)
      )
  )
print(summary_df)
```
# 16. Denovo + ref tree
```sh
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_denovo && cd strep_denovo
# get reference ADS
grep Streptococcus -A 1 ../../ads_operons.fna > ads_operons.filt.fna
sed -i 's/|.*//' ads_operons.filt.fna
sed -i '/^--$/d' ads_operons.filt.fna
# get denovo ADS
grep Streptococcus ~/rna_dohmain/12-denovo_analyses/annotations.txt | awk '{print $1}' | sed 's/_R_//' | while read line; do grep -A 1 $line  ~/rna_dohmain/12-denovo_analyses/arcABC_operons.fna; done > ads_denovo.filt.fna
sed -i 's/|.*//' ads_denovo.filt.fna
sed -i '/^--$/d' ads_denovo.filt.fna

cat ads_denovo.filt.fna ads_operons.filt.fna > ads_operans.comb.fna
mafft  --thread -1 --adjustdirectionaccurately \
  ads_operans.comb.fna > arcABC_operons.align.fna

trimal -in arcABC_operons.align.fna \
  -out arcABC_operons.trim.fna \
  -htmlout arcABC_operons.trim.html \
  -gt 0.5 \
  -resoverlap 0.5 \
  -seqoverlap 50

iqtree -s arcABC_operons.trim.fna -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep.denovo_ref -safe
```
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")

#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_denovo/strep.denovo_ref.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tree.root$tip.label <- gsub("^_R_", "", tree.root$tip.label)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))

#Collapse and average by arc genes
filt_HIV <- filter(dfHIV, Genus == "Streptococcus")
filt_tooth <- filter(dftooth, Genus == "Streptococcus")

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'SEQ_ID', 'Genus_Species', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'SEQ_ID','Genus_Species','padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("SEQ_ID", "Genus_Species"))

# group and get average base mean
group_SEQ <- joined_df %>%
  group_by(SEQ_ID, Genus_Species) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ <- as.data.frame(group_SEQ)

all((group_SEQ$down_HUUvHI > 0) + (group_SEQ$up_HUUvHI > 0) <= 1)
group_SEQ$one_HIvHUU <- ifelse(group_SEQ$down_HUUvHI >0, "HI", ifelse(group_SEQ$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ$down_HUUvHEU > 0) + (group_SEQ$up_HUUvHEU > 0) <= 1)
group_SEQ$one_HEUvHUU <- ifelse(group_SEQ$down_HUUvHEU >0, "HEU", ifelse(group_SEQ$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ$down_HEUvHI > 0) + (group_SEQ$up_HEUvHI > 0) <= 1)
group_SEQ$one_HEUvHI <- ifelse(group_SEQ$down_HEUvHI >0, "HI", ifelse(group_SEQ$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ$down_HvD > 0) + (group_SEQ$up_HvD > 0) <= 1)
group_SEQ$one_HvD <- ifelse(group_SEQ$down_HvD >0, "D", ifelse(group_SEQ$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ$down_HvE > 0) + (group_SEQ$up_HvE > 0) <= 1)
group_SEQ$one_HvE <- ifelse(group_SEQ$down_HvE >0, "E", ifelse(group_SEQ$up_HvE >0, "H", "not_significant"))
all((group_SEQ$down_EvD > 0) + (group_SEQ$up_EvD > 0) <= 1)
group_SEQ$one_EvD <- ifelse(group_SEQ$down_EvD >0, "D", ifelse(group_SEQ$up_EvD >0, "E", "not_significant"))

group_SEQ2 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ2) <- group_SEQ2$SEQ_ID
group_SEQ2 <- rename(group_SEQ2,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ2$SEQ_ID <- NULL

group_SEQ3 <- group_SEQ %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ3) <- group_SEQ3$SEQ_ID
group_SEQ3 <- rename(group_SEQ3,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ3$SEQ_ID <- NULL

#now repeat for denovo data
load("~/rna_dohmain/12-denovo_analyses/denovo_data.RData")

filt_HIV <- dfHIV[grepl("Streptococcus", dfHIV$taxonomy), ]
filt_tooth <- dftooth[grepl("Streptococcus", dftooth$taxonomy), ]

filt_HIV2 <- filt_HIV %>% select('HUU', 'HEU', 'HI', 'seqid', 'taxonomy', 'padj_HUUvHI', 'padj_HUUvHEU', 'padj_HEUvHI', 'log2FoldChange_HUUvHI', 'log2FoldChange_HUUvHEU', 'log2FoldChange_HEUvHI')
filt_tooth2 <- filt_tooth %>% select('H', 'E', 'D', 'seqid', 'taxonomy', 'padj_HvD', 'padj_HvE', 'padj_EvD', 'log2FoldChange_HvD', 'log2FoldChange_HvE', 'log2FoldChange_EvD')

joined_df <- left_join(filt_HIV2, filt_tooth2, by = c("seqid", "taxonomy"))
joined_df <- joined_df %>%
  rename(SEQ_ID = seqid)
group_SEQ_2 <- joined_df %>%
  group_by(SEQ_ID, taxonomy) %>%
  summarise(
    up_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI > 0, na.rm = TRUE),
    down_HUUvHI = sum(padj_HUUvHI < 0.05 & log2FoldChange_HUUvHI < 0, na.rm = TRUE),

    up_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU > 0, na.rm = TRUE),
    down_HUUvHEU = sum(padj_HUUvHEU < 0.05 & log2FoldChange_HUUvHEU < 0, na.rm = TRUE),

    up_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI > 0, na.rm = TRUE),
    down_HEUvHI = sum(padj_HEUvHI < 0.05 & log2FoldChange_HEUvHI < 0, na.rm = TRUE),

    up_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD > 0, na.rm = TRUE),
    down_HvD = sum(padj_HvD < 0.05 & log2FoldChange_HvD < 0, na.rm = TRUE),

    up_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE > 0, na.rm = TRUE),
    down_HvE = sum(padj_HvE < 0.05 & log2FoldChange_HvE < 0, na.rm = TRUE),

    up_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD > 0, na.rm = TRUE),
    down_EvD = sum(padj_EvD < 0.05 & log2FoldChange_EvD < 0, na.rm = TRUE),
    # Average expression
    average_baseMean_HUU = mean(HUU, na.rm = TRUE),
    average_baseMean_HEU = mean(HEU, na.rm = TRUE),
    average_baseMean_HI = mean(HI, na.rm = TRUE),
    average_baseMean_H = mean(H, na.rm = TRUE),
    average_baseMean_E = mean(E, na.rm = TRUE),
    average_baseMean_D = mean(D, na.rm = TRUE),
  )
group_SEQ_2 <- as.data.frame(group_SEQ_2)
group_SEQ_2$SEQ_ID <- gsub("^_R_", "",group_SEQ_2$SEQ_ID)

all((group_SEQ_2$down_HUUvHI > 0) + (group_SEQ_2$up_HUUvHI > 0) <= 1)
group_SEQ_2$one_HIvHUU <- ifelse(group_SEQ_2$down_HUUvHI >0, "HI", ifelse(group_SEQ_2$up_HUUvHI >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HUUvHEU > 0) + (group_SEQ_2$up_HUUvHEU > 0) <= 1)
group_SEQ_2$one_HEUvHUU <- ifelse(group_SEQ_2$down_HUUvHEU >0, "HEU", ifelse(group_SEQ_2$up_HUUvHEU >0, "HUU", "not_significant"))
all((group_SEQ_2$down_HEUvHI > 0) + (group_SEQ_2$up_HEUvHI > 0) <= 1)
group_SEQ_2$one_HEUvHI <- ifelse(group_SEQ_2$down_HEUvHI >0, "HI", ifelse(group_SEQ_2$up_HEUvHI >0, "HEU", "not_significant"))

all((group_SEQ_2$down_HvD > 0) + (group_SEQ_2$up_HvD > 0) <= 1)
group_SEQ_2$one_HvD <- ifelse(group_SEQ_2$down_HvD >0, "D", ifelse(group_SEQ_2$up_HUUvHI >0, "H", "not_significant"))
all((group_SEQ_2$down_HvE > 0) + (group_SEQ_2$up_HvE > 0) <= 1)
group_SEQ_2$one_HvE <- ifelse(group_SEQ_2$down_HvE >0, "E", ifelse(group_SEQ_2$up_HvE >0, "H", "not_significant"))
all((group_SEQ_2$down_EvD > 0) + (group_SEQ_2$up_EvD > 0) <= 1)
group_SEQ_2$one_EvD <- ifelse(group_SEQ_2$down_EvD >0, "D", ifelse(group_SEQ_2$up_EvD >0, "E", "not_significant"))

group_SEQ_21 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_HUU, average_baseMean_HEU, average_baseMean_HI)
row.names(group_SEQ_21) <- group_SEQ_21$SEQ_ID
group_SEQ_21 <- rename(group_SEQ_21,
    HUU = average_baseMean_HUU,
    HEU = average_baseMean_HEU,
    HI = average_baseMean_HI
  )
group_SEQ_21$SEQ_ID <- NULL

group_SEQ_22 <- group_SEQ_2 %>%
  filter(SEQ_ID %in% tree.root$tip.label) %>% select(SEQ_ID, average_baseMean_H, average_baseMean_E, average_baseMean_D)
row.names(group_SEQ_22) <- group_SEQ_22$SEQ_ID
group_SEQ_22 <- rename(group_SEQ_22,
    H = average_baseMean_H,
    E = average_baseMean_E,
    D = average_baseMean_D
  )
group_SEQ_22$SEQ_ID <- NULL

group_SEQ_2$Genus_Species <- sapply(group_SEQ_2$taxonomy, function(tax) {
  # If it starts with "Streptococcus_sp.", leave it unchanged
  if (grepl("^Streptococcus_sp\\.", tax)) {
    tax
  } else {
    # Otherwise split by "_" and keep first two parts if possible
    parts <- strsplit(tax, "_")[[1]]
    if (length(parts) >= 2) {
      paste(parts[1:2], collapse = "_")
    } else {
      parts[1]
    }
  }
})
group_SEQ_2$taxonomy <- NULL

# combine
group1 <- rbind(group_SEQ_21, group_SEQ2)
group2 <- rbind(group_SEQ_22, group_SEQ3)
group_overall <- rbind(group_SEQ_2, group_SEQ)
group_overall <- group_overall %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID)                  
  )
# pruned_tree <- drop.tip(tree.root, "SEQF1998.1")


increased_by <- .008
first_tip <- 2.11
second_tip <- first_tip + increased_by
third_tip <- second_tip +increased_by
fourth_tip <- third_tip +increased_by
fifth_tip <- fourth_tip +increased_by
sixth_tip <- fifth_tip +increased_by

p <- ggtree(tree.root) %<+% group_overall +
  geom_tiplab(aes(label=species2), align = TRUE, size = 3.5, offset = 0.0009) +
  # scale_color_manual(
  #   name = "Taxonomy",
  #   values = c(
  #     "Streptococcus_anginosus" = "#e6194b",
  #     "Streptococcus_oralis" = "#3cb44b",
  #     "Streptococcus_sanguinis" = "#ffe119",
  #     "Streptococcus_cristatus" = "#4363d8",
  #     "Streptococcus_parasanguinis" = "#f58231",
  #     "Streptococcus_sp._NPS_308" = "#911eb4",
  #     "Streptococcus_sp._A12" = "#46f0f0",
  #     "Streptococcus_gordonii" = "#f032e6",
  #     "Streptococcus_pneumoniae" = "#bcf60c",
  #     "Streptococcus_intermedius" = "#fabebe",
  #     "Streptococcus_sp._HMT_056" = "#008080",
  #     "Streptococcus_sp._HMT_066" = "#e6beff",
  #     "Streptococcus_australis" = "#9a6324",
  #     "Streptococcus_constellatus" = "#fffac8",
  #     "Streptococcus_sinensis" = "#800000",
  #     "Streptococcus_infantis" = "#aaffc3",
  #     "Streptococcus_mutans" = "#808000",
  #     "Streptococcus_mitis" = "#ffd8b1",
  #     "Streptococcus_lactarius" = "#000075",
  #     "Streptococcus_downei" = "#808080",
  #     "Streptococcus_salivarius" = "#000000"
  #   )
  # ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHUU, x = first_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HUU vs HEU",
    values = c("HEU" = "#FA78FA", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HIvHUU, x = second_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HUU",
    values = c("HI" = "#8213A0", "HUU" = "#40A0FA", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HEUvHI, x = third_tip), size = 6, shape = 16) +
  scale_color_manual(
    name = "HI vs HEU",
    values = c("HEU" = "#FA78FA", "HI" = "#8213A0", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvE, x = fourth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs E",
    values = c("H" = "#24B45A", "E" = "#F0F032", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_HvD, x = fifth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "H vs D",
    values = c("D" = "#AA0A3B", "H" = "#24B45A", "not_significant" = "#7F7F7F")
  ) +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = one_EvD, x = sixth_tip), size = 6, shape = 15) +
  scale_color_manual(
    name = "E vs D",
    values = c("E" = "#F0F032", "D" = "#AA0A3B", "not_significant" = "#7F7F7F")
  ) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16),
   legend.title = element_text(size = 18, face = "bold"),
   legend.key.size = unit(1.2, "cm"),
   legend.spacing = unit(0.6, "cm")
  )

p1 <- gheatmap(p,
         group1,
         width = 0.02,
         offset = 0.135,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.7,
         font.size = 5,
         color = "black")

pdf("tree.denovo_streps.pdf", height = 120, width =80)
  gheatmap(p1,
         group2,
         width = 0.02,
         offset = 0.180,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.denovo_streps.pdf")







p1 <- gheatmap(p,
         log10(group1),
         width = 0.02,
         offset = 0.13,
         colnames_position = "top",
         colnames_angle = 45,
         colnames_offset_y =.7,
         font.size = 5,
         color = "black")

pdf("tree.denovo_streps.log.pdf", height = 120, width =80)
  gheatmap(p1,
         log10(group2),
         width = 0.02,
         offset = 0.178,
         colnames_position = "top",
         # colnames_angle = 45,
         # colnames_offset_y =.8,
         font.size = 5,
         color = "black") +
  scale_fill_gradient(low = "white", high = "black", name = "Average basemean")
dev.off()
system("~/.iterm2/imgcat tree.denovo_streps.log.pdf")

# group2[row.names(group2) == "SEQF6415.1", ]
```
# 17. Do tree for every strep even ones without ADS
```sh
# strart making phylogenies
cd ~/rna_dohmain/07-ads_expression/vince_trees
mkdir strep_all && cd strep_all
mkdir gff_files
grep -Ew 'Streptococcus' ../../../homd_map/annotations.merge.txt \
| awk '{print $3}' \
| sort -u \
| parallel -j 190 'cp ../../../homd_map/gff_files/{}.gff gff_files/{}.gff'
# run panaroo
panaroo -i ./gff_files/*gff -o results --clean-mode strict -t 190 -a core --aligner mafft -c 0.5 --core_threshold 1.00
numb_refs=$(ls gff_files/*gff | wc -l)

mkdir core_genes
grep -c ">" ./results/aligned_gene_sequences/*fas | grep ":$numb_refs" | sed 's/.\/results\/aligned_gene_sequences\///' | sed 's/:.*//' | while read line; do cp ./results/aligned_gene_sequences/$line ./core_genes/$line; done
cd core_genes
ls ./*fas | sed 's/.\/results\/aligned_gene_sequences\///' | while read line; do grep -n ">" ./$line | sed 's/;.*//' |sort |uniq | wc -l; done | grep "$numb_refs" -c
ls *fas | wc -l
# test for recombination
parallel -j 190 'Phi -o -f {} > {.}.rec' ::: *aln.fas

grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb1

grep "Max" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb2

grep "NSS" *.rec | grep permutations | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb3

cat <(awk '{print $1}' recomb1) <(awk '{print $1}' recomb2) <(awk '{print $1}' recomb3)| sed 's/rec.*/fas/' | sort | uniq -c | grep -w 3 | awk '{print $2}' | while read line; do mv $line $line.rd; done
sed -i 's/;.*//' *fas
sed -i 's/_R_//' *fas

# check alignments are all equal length
ls *fas | while read line; do samtools faidx $line; done
ls *fai | while read line; do awk '{print $2}' $line | sort | uniq | wc -l; done

# combine alignments
python3 ../../panaroo_0.6/core_genes/combine_core.py
grep ">" core_genome.align.fa -c
samtools faidx core_genome.align.fa
# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix strep_all.core -safe
```
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)

setwd("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees")
load("tree_data.RData")
homd <- read.csv("~/rna_dohmain/homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
#make tree
tree <- read.tree("/home/suzanne/rna_dohmain/07-ads_expression/vince_trees/strep_all/core_genes/strep_all.core.treefile")
# tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tree.root$tip.label <- gsub("^_R_", "", tree.root$tip.label)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# get the tip labels that have ADS for Strep
filt_HIV <- filter(dfHIV, Genus == "Streptococcus")
group_SEQ <- filt_HIV %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID)                  
  )

# make annotation file 
resdf <- as.data.frame(tree.root$tip.label)
names(resdf)[1] <- "SEQ_ID"
row.names(resdf) <- resdf$SEQ_ID
homd$locus_tag <- NULL
homd$gene <- NULL
homd$gene_base <- NULL
homd2 <- distinct(homd)
ann <- homd2[homd2$SEQ_ID %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$SEQ_ID
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
resdf <- resdf[ , -1]
resdf <- resdf %>%
  mutate(
    species_part = str_replace(Genus_Species, "^[^_]+_", ""),
    species2 = paste0(species_part, "_", SEQ_ID),
    species2 = str_remove(species2, "subsp\\.?_oralis_?"))
resdf$SEQ_ID

# check to see which have ADS
resdf <- resdf %>%
  mutate(has_ADS = SEQ_ID %in% filt_HIV$SEQ_ID)
# add boostrapping
bs_values <- as.numeric(tree.root$node.label)
p <- ggtree(tree.root) %<+% resdf +
  geom_tiplab(aes(label=species2, color = has_ADS), align = TRUE, size = 3, offset = 0.0) +
  scale_color_manual(values = c("TRUE" = "#F8766D", "FALSE" = "gray"))

pdf("tree.strep_all.pdf", height = 135, width = 80)
ggtree(tree.root) %<+% resdf +
  geom_tiplab(aes(label = species2, color = has_ADS), align = FALSE, size = 2, offset = 0.0) +
  scale_color_manual(values = c("TRUE" = "#F8766D", "FALSE" = "gray")) +
  geom_text2(aes(subset = !isTip & !is.na(as.numeric(label)), label = label), hjust = -0.3, size = 3)
dev.off()
system("~/.iterm2/imgcat tree.strep_all.pdf")
```
## 17.1 Filter tree based on completeness of genome
```sh
cd ~/rna_dohmain/homd_map
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

grep -w "Streptococcus" SEQID_info.csv | awk -F "," '{print $1, $3, $4, $8}' | sed 's/"//g' | sed 's/htt.*GCA/GCA/' | sed 's/_AS.*//' | sed 's/ /\t/g' > strep_homd
awk '{print $4}' strep_homd | while read line; do grep -wm1 $line assembly_summary.txt | awk -F "\t" '{print $12}'; done | sort | uniq -c 
awk '{print $3}' ~/rna_dohmain/07-ads_expression/vince_trees/panaroo_test/gff_files/ads_annotations.tsv | sort | uniq | while read line; do grep -wm 1 $line strep_homd; done > ads_homd
awk '{print $4}' ads_homd | while read line; do grep -wm1 $line assembly_summary.txt | awk -F "\t" '{print $12}'; done | sort | uniq -c 

#calculate n50
cd ~/rna_dohmain/homd_map/fna_files
wget -q -O - https://www.homd.org/ftp/genomes/PROKKA/V10.1/fna/ | \
grep -oP 'SEQ[^"]+\.fna' | \
parallel -j 190 wget https://www.homd.org/ftp/genomes/PROKKA/V10.1/fna/{}

```
# 18. Sort by N50
```sh
cd ~/rna_dohmain/homd_map/fna_files
grep -Ew 'Streptococcus' ../annotations.merge.txt | awk '{print $3}' | sort -u | while read line; do
  echo -e "$line\t$(n50 "$line.fna")"
done > strep_n50.txt
sort -k2,2nr strep_n50.txt > strep_n50_sorted.txt
tail strep_n50_sorted.txt
awk '{ sum += $2; count++ } END { if (count > 0) print sum / count; else print 0 }' strep_n50_sorted.txt
# 688076
grep -Ew 'Streptococcus' ../ads_homd | awk '{print $1}' | sort -u | while read line; do
  echo -e "$line\t$(n50 "$line.fna")"
done > strep_ads_n50.txt
sort -k2,2nr strep_ads_n50.txt > strep_ads_n50_sorted.txt
awk '{ sum += $2; count++ } END { if (count > 0) print sum / count; else print 0 }' strep_ads_n50_sorted.txt
```
