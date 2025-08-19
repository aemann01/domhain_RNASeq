# 2, Run DESeq HUU vs HEU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-t-denticola")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./tdent_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
ptm <- proc.time()
#set HUU as reference
star_results$hiv_status <- relevel(star_results$hiv_status, ref = "HUU")
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.05:  5392"
# out of 18242 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5138, 28%
# LFC < 0 (down)     : 254, 1.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 6012, 33%
# (mean count < 2)
```
```R
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
vld <- varianceStabilizingTransformation(se_star)
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 50)
#make a subset of the log transformed counts for just the top 25 varying genes
topCounts <- assay(vld)[topVarGenes,]
df <- as.data.frame(colData(vld)[,c("hiv_status","tooth_health")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(hiv_status = c("#8213A0","#FA78FA","#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HEU", "HUU")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V8=="Treponema_denticola")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "tooth_health", "rel_abundance")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$rel_abundance <- as.numeric(merged_df$rel_abundance)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("tdent.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.test.tdent.pdf")
system("~/.iterm2/imgcat ./heatmap.test.tdent.pdf")
```
By top gene
```R
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("tdent.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
group_all <- long_df %>% group_by(Var2, Var1) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = Var1,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$Var1
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.test.tdent.pdf")
system("~/.iterm2/imgcat ./heatmap.test.tdent.pdf")
```