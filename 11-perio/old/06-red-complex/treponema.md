# 1. Get reads from Treponema
```sh
grep "Treponema" ~/rna_dohmain/homd_map/annotations.merge.txt | grep -v none > trep.annotations.txt
cat <(head -n 1 ~/rna_dohmain/homd_map/annotations.merge.txt) trep.annotations.txt > temp
mv temp trep.annotations.txt
parallel -a <(awk '{print $1}' trep.annotations.txt | sed '1 i\Geneid') -j 7 -k "grep -w '{}' ../../homd_map/read_counts.txt" > trep_counts.txt
```
# 2. Look at rpoC distro
```R
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(microbiome)
setwd("~/rna_dohmain/11-perio/06-red-complex")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V7=="Treponema")
data <- psmelt(actino) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("treponema.dna.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./treponema.dna.pdf")


data %>% group_by(hiv_status) %>% summarise(average = mean(c(Abundance)))
aggregate(Abundance ~ hiv_status, data = data, FUN = mode)

pdf("treponema.hist.pdf", width =15, heigh =10)
ggplot(data, aes(x = Abundance)) +
  geom_histogram(binwidth = 0.001, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./treponema.hist.pdf")

#look at species implicated in perio
actino <- subset_taxa(rel, V8=="Treponema_denticola" | V8=="Treponema_parvum"|V8=="Treponema_medium" | V8=="Treponema_maltophilum" | V8=="Treponema_lecithinolyticum" | V8=="Treponema_amylovorum" | V8=="Treponema_socranskii")
data <- psmelt(actino) # create dataframe from phyloseq object
# data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("treponema_prio.dna.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./treponema_prio.dna.pdf")

#HUU precent with Treponema
sub_data <- data %>% group_by(Sample, hiv_status, V7) %>% summarise(Abundance = sum(c(Abundance)))
abundance_summary <- sub_data %>%
  group_by(hiv_status) %>%
  summarise(total_samples = n(),  # Total samples in each group
            samples_above_0 = sum(Abundance > 0),  # Samples with Abundance > 0
            percent_above_0 = (samples_above_0 / total_samples) * 100)  # Percentage
abundance_summary
abundance_summary <- sub_data %>%
  group_by(hiv_status) %>%
  summarise(total_samples = n(),  # Total samples in each group
            samples_above_0 = sum(Abundance > 0.01),  # Samples with Abundance > 0
            percent_above_0 = (samples_above_0 / total_samples) * 100)  # Percentage
abundance_summary
```
# Look at total RNA distro by species
```R
library(tidyverse)
library(reshape2)
library(ggolot2)
read_counts <- t(read.csv("~/rna_dohmain/09-urease/09-global-distro/species_reads.txt", header=T, row.names=1,  sep = "\t"))
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)

long_df <- melt(read_counts)
long_df$genus <- gsub(x = long_df$Var1, pattern = "_.*", replacement = "") 
sub_df <- long_df[long_df$genus == "Treponema",]
sub_df <- left_join(sub_df, metadata, by = join_by(Var2 ==  sample_id))
subcount <- sub_df[, colnames(sub_df) %in% c("Var1","Var2","value","hiv_status","tooth_health")]

# plot
pdf("treponema.rna.pdf", width =15, heigh =10)
ggplot(sub_df)+
  geom_bar(aes(x=Var2, y=value,fill=Var1),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./treponema.rna.pdf")
```
Look at relationship between RNA and DNA
```R
library(bcdstats)
library(ggstatsplot)
library(ggside)
sub_data <- data[, colnames(data) %in% c("Sample","V8","hiv_status","Abundance","tooth_health")]
names(sub_data) <- c("sample", "abundance","tooth_health", "hiv_status","species")
names(subcount) <- c("species", "sample","abundance", "tooth_health","hiv_status")
sub_data$nucl <- "DNA"
sub_data$abundance <- sub_data$abundance*100
subcount$nucl <- "RNA"

combined_df_columns <- rbind(sub_data, subcount)
pdf("dnavrna.trep.pdf")
ggplot(combined_df_columns,aes(x = sample,y = abundance)) + 
    geom_bar(aes(fill = nucl),stat = "identity",position = "dodge") + 
     facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./dnavrna.trep.pdf")

#run a correlation
agg_dna <- sub_data %>% group_by(sample, hiv_status) %>% summarise(abundance_dna = sum(c(abundance)))
agg_rna <- subcount %>% group_by(sample, hiv_status) %>% summarise(abundance_dna = sum(c(abundance)))
merged_nucs <- inner_join(agg_dna, agg_rna, by = c("sample" = "sample", "hiv_status" = "hiv_status"))
names(merged_nucs) <- c("sample", "hiv_status","abundance_dna", "abundance_rna")
test <- cor(merged_nucs$abundance_dna, merged_nucs$abundance_rna, method = "spearman")
adjust.corr(test)
hiv_colors <- c("#FA78FA","#8213A0","#40A0FA")
pdf("dnavrna.trep.pdf")
ggplot(merged_nucs, aes(x = abundance_dna, y = abundance_rna)) +
  geom_point(aes(color = hiv_status, size = 1))+
  scale_color_manual(values = hiv_colors) +
  geom_xsidedensity(
    aes(
      y    = ..count../sum(..count..),
    ),
    alpha    = 0.5,
    size     = 1,
    fill="cyan4"
  ) +
  scale_xsidey_continuous(minor_breaks = NULL)+
  scale_xfill_manual(values = "cyan4")+
  geom_ysidedensity(
    aes(
      x    = ..count../sum(..count..),
    ),
    alpha    = 0.5,
    size     = 1,
    fill="orange"
  )+
  stat_cor(method="spearman")+
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
  theme_bw()
dev.off()
system("~/.iterm2/imgcat ./dnavrna.trep.pdf")

```
# 3. Run DESeq HUU vs HI
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  16111"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 831, 1.2%
# LFC < 0 (down)     : 15280, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 32803, 48%
# (mean count < 1)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  16111"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1560, 2.3%
# LFC < 0 (down)     : 16783, 25%
# outliers [1]       : 0, 0%
# low counts [2]     : 32803, 48%
# (mean count < 1)
write.table(resLFC, file="deseq_results_trep-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_trep-HUUvHI.RData")
```
Valcona Plot
```R
load("deseq_results_trep-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-HUUvHI.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUUvHI.trep.pdf")

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
pdf("volcano-HUUvHI.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUUvHI.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HUUvHI.trep.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HUUvHI.trep.pdf")
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
df <- as.data.frame(colData(vld)[,c("hiv_status","tooth_health", "cd4_group")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(hiv_status = c("#8213A0","#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HUU")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "tooth_health","cd4_group", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HUUvHI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUUvHI.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HUU_v_HI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUU_v_HI.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HUU_v_HI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUU_v_HI.trep.pdf")
```
# 3, Run DESeq HUU vs HEU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  15300"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 511, 0.75%
# LFC < 0 (down)     : 14789, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 31734, 47%
# (mean count < 1)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  15300"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1270, 1.9%
# LFC < 0 (down)     : 16403, 24%
# outliers [1]       : 0, 0%
# low counts [2]     : 31734, 47%
# (mean count < 1)
write.table(resLFC, file="deseq_results_trep-HUUvHEU.txt", quote=F, sep="\t")
save.image("deseq_results_trep-HUUvHEU.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-HUUvHEU.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUUvHEU.trep.pdf")

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
pdf("volcano-HUUvHEU.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUUvHEU.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HUUvHEU.trep.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HUUvHEU.trep.pdf")
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
mat_colors <- list(hiv_status = c("#FA78FA","#40A0FA"))
names(mat_colors$hiv_status) <- c("HEU", "HUU")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "tooth_health", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HUUvHEU.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUUvHEU.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HUU_v_HEU.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUU_v_HEU.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HUU_v_HEU.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HUU_v_HEU.trep.pdf")
```
# 3, Run DESeq HEU vs HI
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HEU",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  3234"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1769, 2.6%
# LFC < 0 (down)     : 1465, 2.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 30198, 44%
# (mean count < 1)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  3234"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2235, 3.3%
# LFC < 0 (down)     : 2468, 3.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 30198, 44%
# (mean count < 1)
write.table(resLFC, file="deseq_results_trep-HEUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_trep-HEUvHI.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-HEUvHI.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHI.trep.pdf")

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
pdf("volcano-HEUvHI.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHI.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HEUvHI.trep.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HEUvHI.trep.pdf")
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
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "tooth_health", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HEUvHI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEUvHI.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HEUvHI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEUvHI.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HEUvHI.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEUvHI.trep.pdf")
```




























































































# 4. Run DESeq HI: low v normal CD4 counts
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI",]
submap <- submap[submap$cd4_group == "low" | submap$cd4_group == "normal",]

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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~cd4_group)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$cd4_group <- factor(star_results$cd4_group, levels=c("low", "normal"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  5153"
# out of 32252 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3739, 12%
# LFC < 0 (down)     : 1414, 4.4%
# outliers [1]       : 7122, 22%
# low counts [2]     : 0, 0%
# (mean count < 1)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="cd4_group_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  5153"
# out of 32252 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5552, 17%
# LFC < 0 (down)     : 1551, 4.8%
# outliers [1]       : 7122, 22%
# low counts [2]     : 0, 0%
# (mean count < 1)
write.table(resLFC, file="deseq_results_trep-HI-lowvnormal.txt", quote=F, sep="\t")
save.image("deseq_results_trep-HI-lowvnormal.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-HI-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HI-lowvnormal.trep.pdf")

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
pdf("volcano-HI-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HI-lowvnormal.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HI_.trep.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HUUvHEU.trep.pdf")
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
df <- as.data.frame(colData(vld)[,c("tooth_health", "cd4_group")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(cd4_group = c("orange","cyan4"))
names(mat_colors$cd4_group) <- c("low", "normal")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","tooth_health", "cd4_group", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HI.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HI.lowvhigh.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HI.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HI.lowvhigh.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HI.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HI.lowvhigh.trep.pdf")
```
# 4. Run DESeq HEU: low v normal CD4 counts
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU",]
submap <- submap[submap$cd4_group == "low" | submap$cd4_group == "normal",]

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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~cd4_group)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$cd4_group <- factor(star_results$cd4_group, levels=c("low", "normal"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  3709"
# out of 34179 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3653, 11%
# LFC < 0 (down)     : 56, 0.16%
# outliers [1]       : 561, 1.6%
# low counts [2]     : 21524, 63%
# (mean count < 12)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="cd4_group_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  3279"
# out of 34179 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5634, 16%
# LFC < 0 (down)     : 117, 0.34%
# outliers [1]       : 561, 1.6%
# low counts [2]     : 19654, 58%
# (mean count < 7)
write.table(resLFC, file="deseq_results_trep-HEU-lowvnormal.txt", quote=F, sep="\t")
save.image("deseq_results_trep-HEU-lowvnormal.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-HEU-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEU-lowvnormal.trep.pdf")

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
pdf("volcano-HEU-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEU-lowvnormal.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HEU_lowvnormal.trep.pdf")
plotPCA(vld, intgroup=c("cd4_group")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HEU_lowvnormal.trep.pdf")
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
df <- as.data.frame(colData(vld)[,c("tooth_health", "cd4_group")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(cd4_group = c("orange","cyan4"))
names(mat_colors$cd4_group) <- c("low", "normal")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","tooth_health", "cd4_group", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HEU.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEU.lowvhigh.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HI.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HI.lowvhigh.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HEU.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEU.lowvhigh.trep.pdf")
```
# 6. Run DESeq overall low v normal CD4 counts
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./trep_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$cd4_group == "low" | metadata$cd4_group == "normal",]

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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~cd4_group)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
#set normal as reference
star_results$cd4_group <- factor(star_results$cd4_group, levels=c("low", "normal"))
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
# [1] "number of genes with adjusted p value lower than 0.05:  2129"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 960, 1.4%
# LFC < 0 (down)     : 1169, 1.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 30153, 44%
# (mean count < 1)

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="cd4_group_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  2129"
# out of 68201 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1550, 2.3%
# LFC < 0 (down)     : 1643, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 30153, 44%
# (mean count < 1)
write.table(resLFC, file="deseq_results_trep-lowvnormal.txt", quote=F, sep="\t")
save.image("deseq_results_trep-lowvnormal.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
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
sigsp <- paste("x", sigloc$Species, sep="_")
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
res_ord$Species[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Species)

#combine species and gene
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
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


res_sub <- res_ord %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene == "arcA")

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
pdf("volcano-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-lowvnormal.trep.pdf")

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
pdf("volcano-lowvnormal.trep.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-lowvnormal.trep.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_lowvnormal.trep.pdf")
plotPCA(vld, intgroup=c("cd4_group")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_lowvnormal.trep.pdf")
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
df <- as.data.frame(colData(vld)[,c("tooth_health", "cd4_group")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(cd4_group = c("orange","cyan4"))
names(mat_colors$cd4_group) <- c("low", "normal")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[6])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V7=="Treponema")
abund <- t(as.data.frame(otu_table(sub)))
merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","tooth_health", "cd4_group", "log10")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$log10 <- as.numeric(merged_df$log10)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.lowvhigh.trep.pdf")
```
# Top genes by average base mean
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
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
save_pheatmap_pdf(x, "heatmap.HI.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HI.lowvhigh.trep.pdf")
```
# Top virulence factors
```R
library(reshape2)
library(tidyr)

gene_counts <- assay(vld)
#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("trep.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$Species, "_", rows2$Gene)
rows2$gene <- paste0(rows2$Gene)
row <- rows2[, c("seq", "name", "gene")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
long_df <- left_join(row, long_df,by = join_by(name ==  Var1))
long_sub <- long_df %>% filter(gene == "oppA" | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA")
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "seq")])
long_sub <- as.data.frame(long_sub[, !(colnames(long_sub) %in% "gene")])
group_all <- long_sub %>% group_by(Var2, name) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = name,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 50)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$name
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.HEU.lowvhigh.trep.pdf")
system("~/.iterm2/imgcat ./heatmap.HEU.lowvhigh.trep.pdf")
```