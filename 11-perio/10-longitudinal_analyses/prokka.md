# 1. DESEQ for V1
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../../homd_map/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "1",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  456134"
# out of 2297042 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 237446, 10%
# LFC < 0 (down)     : 218688, 9.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 756754, 33%
# (mean count < 2)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  438252"
# out of 2297042 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 293901, 13%
# LFC < 0 (down)     : 303646, 13%
# outliers [1]       : 0, 0%
# low counts [2]     : 623482, 27%
# (mean count < 2)

# https://support.bioconductor.org/p/63567/
baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

write.table(resLFC, file="deseq_results_V1-HIvHUU.prokka.txt", quote=F, sep="\t")
save.image("deseq_results_V1-HIvHUU.prokka.RData")
```
# 2. DESEQ for V2
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../../homd_map/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "2",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  456134"
# out of 2297042 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 237446, 10%
# LFC < 0 (down)     : 218688, 9.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 756754, 33%
# (mean count < 2)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  130752"
# out of 713739 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 79770, 11%
# LFC < 0 (down)     : 99675, 14%
# outliers [1]       : 0, 0%
# low counts [2]     : 179886, 25%
# (mean count < 2)

# https://support.bioconductor.org/p/63567/
baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL

write.table(combined_df, file="deseq_results_V2-HIvHUU.prokka.txt", quote=F, sep="\t")
save.image("deseq_results_V2-HIvHUU.prokka.RData")
```
# 3. DESEQ for V3
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../../homd_map/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "3",]
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

# https://support.bioconductor.org/p/63567/
baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )

# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.05:  456134"
# out of 2297042 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 237446, 10%
# LFC < 0 (down)     : 218688, 9.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 756754, 33%
# (mean count < 2)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  130752"
# out of 713739 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 79770, 11%
# LFC < 0 (down)     : 99675, 14%
# outliers [1]       : 0, 0%
# low counts [2]     : 179886, 25%
# (mean count < 2)

# get basemean by group
# https://support.bioconductor.org/p/63567/
baseMeanPerLvl <- sapply( levels(se_star$hiv_status), function(lvl) rowMeans( counts(se_star,normalized=TRUE)[,se_star$hiv_status== lvl] ) )
resLFC_df <- as.data.frame(resLFC)
baseMeanPerLvl_df <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl_df$gene_id <- rownames(as.data.frame(baseMeanPerLvl_df))
resLFC_df$gene_id <- rownames(resLFC_df)
combined_df <- merge(resLFC_df, baseMeanPerLvl_df, by = "gene_id", all.x = TRUE)
rownames(combined_df) <- combined_df$gene_id
combined_df$gene_id <- NULL
write.table(combined_df, file="deseq_results_V3-HIvHUU.prokka.txt", quote=F, sep="\t")
save.image("deseq_results_V3-HIvHUU.prokka.RData")
```
# 4. P. gingivalis trees
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
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "kgp", "fimA", "rpoC")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/prokka/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
for (virulence2 in virulence_genes_pg) {
resLFC <- read.csv("./deseq_results_V1-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
# add in annotations
homd <- read.csv("../../homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
homd$species <- paste0(homd$Genus,"_", homd$Species)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V1_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

# V2
resLFC <- read.csv("./deseq_results_V2-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V2_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

# V3
resLFC <- read.csv("./deseq_results_V3-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V3_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  SEQ_ID = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("SEQ_ID", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(SEQ_ID, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_wrap(~type) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  facet_wrap(~type, ncol =6,  labeller = labeller(type = facet_labels)) +
  # ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  # geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = padj), color = "black") +
  # scale_fill_distiller(type = "seq", palette = "Greys", direction = -1, 
  #                      guide = guide_colorbar(title = "pvalue", title.position = "top")) +
  # facet_wrap(~type, ncol =9,  labeller = labeller(type = facet_labels)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)
head(combined_metrics_complete)
pdf(paste0("./pging_tree_", virulence2, "_prokka.pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_tree_", virulence2 , "_prokka.pdf", sep=""))
}
```
# 5. T. denticola tree
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
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_td <- c("oppA","flaB1", "flaB3", "flaB3", "fliE", "fliF", "fliG", "fliM", "fliY", "flgC", "flgE", "motB", "cheX", "cheY", "hbpA", "troA")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/prokka/t_dent/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
for (virulence2 in virulence_genes_td) {
resLFC <- read.csv("./deseq_results_V1-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
# add in annotations
homd <- read.csv("../../homd_map/annotations.merge.txt", header=T, sep="\t", quote="")
homd$species <- paste0(homd$Genus,"_", homd$Species)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Treponema_denticola") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V1_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

# V2
resLFC <- read.csv("./deseq_results_V2-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Treponema_denticola") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V2_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

# V3
resLFC <- read.csv("./deseq_results_V3-HIvHUU.prokka.txt", header=T, sep="\t", row.names=1)
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

res_ord <- resdf %>% filter(gene_base != "none")
#combine species and gene
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene_base)
res_ord$gene_name <- gsub(x = res_ord$gene_base, pattern = "_.", replacement = "") 
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene_base ==virulence2 & species == "Treponema_denticola") %>%
  group_by(SEQ_ID, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(SEQ_ID, species, avg_log2FoldChange_factor, avg_padj) %>%
  mutate(metric = "avg_log2FoldChange_factor") %>%
  rename(value = avg_log2FoldChange_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)

sig_average_baseMean_V3_virus <- average_pg %>%
  select(SEQ_ID, species, avg_baseMean_factor, avg_padj) %>%
  mutate(metric = "avg_baseMean_factor") %>%
  rename(value = avg_baseMean_factor,
         padj = avg_padj) %>%
  select(SEQ_ID, species, metric, value, padj)


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  SEQ_ID = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("SEQ_ID", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(SEQ_ID, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_wrap(~type) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  facet_wrap(~type, ncol =6,  labeller = labeller(type = facet_labels)) +
  # ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  # geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = padj), color = "black") +
  # scale_fill_distiller(type = "seq", palette = "Greys", direction = -1, 
  #                      guide = guide_colorbar(title = "pvalue", title.position = "top")) +
  # facet_wrap(~type, ncol =9,  labeller = labeller(type = facet_labels)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)
head(combined_metrics_complete)
pdf(paste0("./tdent_tree_", virulence2, "_prokka.pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./tdent_tree_", virulence2 , "_prokka.pdf", sep=""))
}
```

