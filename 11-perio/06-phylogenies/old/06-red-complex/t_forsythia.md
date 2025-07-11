# Find which samples have more than rel abudannce 1% of T. forsythia
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
setwd("~/rna_dohmain/11-perio/06-red-complex")
#get relative abundance of t. forsythia
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
rel <- microbiome::transform(glom, "compositional")
actino <- subset_taxa(rel, V8=="Tannerella_forsythia")
data <- psmelt(actino) # create dataframe from phyloseq object
samples <- data %>% filter(Abundance > 0.00)
sample_list <- samples$Sample
samples <- data$Sample<- factor(data$Sample,levels=unique(data$Sample))
```
# 2, Run DESeq HI vs HEU for samples that have T. forsythia present
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../03-t-forsythia/tforsythia_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
submap <- submap[submap$sample_id %in% sample_list, ]
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
# [1] "number of genes with adjusted p value lower than 0.05:  924"

# out of 10318 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 344, 3.3%
# LFC < 0 (down)     : 580, 5.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 3801, 37%
# (mean count < 3)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  5875"

# out of 26924 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2826, 10%
# LFC < 0 (down)     : 5137, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 13558, 50%
write.table(resLFC, file="deseq_results_tforsythia-HEUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_tforsythia-HEUvHI.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("../03-t-forsythia/tforsythia.annotations.txt", header=T, sep="\t", quote="")
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
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$SEQ_ID,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 10)
# negative top 10
top <- tail(sortdf$locus_tag, 10)
# concatenate
labgenes <- c(top, low)


res_sub <- res_ord %>% filter(gene == "susB" | gene == "kly" | gene == "eno" )
| gene == "hagA" | gene == "fimA")

res_ord %>% filter(gene == "eno")

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
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
pdf("volcano-HEUvHI.tforsythia.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHI.tforsythia.pdf")
#Create volcano plot
overall_plot <- EnhancedVolcano(res_sub,
	lab = res_sub$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	# colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	max.overlaps = Inf,
	pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HEUvHI.tforsythia.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHI.tforsythia.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HEUvHI.tforsythia.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HEUvHI.tforsythia.pdf")
```
Make dendogram of rpoC activity
```R
#Get top varying genes
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
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
mat_colors <- list(hiv_status = c("#8213A0","#FA78FA"))
names(mat_colors$hiv_status) <- c("HI", "HEU")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V8=="Tannerella_forsythia")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "tooth_health", "rel_abundance")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$rel_abundance <- as.numeric(merged_df$rel_abundance)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("tforsythia.annotations.txt", header=T, sep="\t", quote="")
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
save_pheatmap_pdf(x, "heatmap.HEU_v_HI.tforsythia.pdf")
system("~/.iterm2/imgcat ./heatmap.HEU_v_HI.tforsythia.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("../03-t-forsythia/tforsythia.annotations.txt", header=T, sep="\t", quote="")
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
save_pheatmap_pdf(x, "heatmap.HEU_v_HI.tforsythia.pdf")
system("~/.iterm2/imgcat ./heatmap.HEU_v_HI.tforsythia.pdf")
```