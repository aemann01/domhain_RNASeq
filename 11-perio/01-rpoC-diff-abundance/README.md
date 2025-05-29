# 1. Subset annotation file
```sh
grep rpoC ~/rna_dohmain/homd_map/annotations.merge.txt > rpoC.annotations.txt
cat <(head -n 1 ~/rna_dohmain/homd_map/annotations.merge.txt) rpoC.annotations.txt > temp
mv temp rpoC.annotations.txt
```
# 2. Compare HUU to HI
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)


#load data
setwd("~/rna_dohmain/11-perio/01-rpoC-diff-abundance")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("~/rna_dohmain/09-urease/09-global-distro/rpoC_counts.txt", header=T, sep="\t", row.names=1)
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
# [1] "number of genes with adjusted p value lower than 0.05:  324"

# out of 8642 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 255, 3%
# LFC < 0 (down)     : 69, 0.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 5656, 65%

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  337"
# out of 8642 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 292, 3.4%
# LFC < 0 (down)     : 115, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 5993, 69%
write.table(resLFC, file="deseq_results_rpoC-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_rpoC-HUUvHI.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("rpoC.annotations.txt", header=T, sep="\t", quote="")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 1,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HUUvHI.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUUvHI.pdf")
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HUUvHI_rpoC.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HUUvHI_rpoC.pdf")
```
Make dendogram of rpoC activity
```R
#Get top varying genes
library(pheatmap)
library(RColorBrewer)
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
mat_colors <- list(hiv_status = c("#8213A0","#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HUU")
mat_colors$tooth_health <- c("#F35E5A","#25AE2B" "#4F87FF")
names(mat_colors$tooth_health) <- c("D", "E", "H")
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("rpoC.annotations.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  locus_tag))
rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name
x <- pheatmap(topCounts, annotation_col = df, annotation_colors = mat_colors,  color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "test.pdf")
system("~/.iterm2/imgcat ./test.pdf")
```
# 2. rpoB
Subset data
```sh
grep rpoB ~/rna_dohmain/homd_map/annotations.merge.txt > rpoB.annotations.txt
cat <(head -n 1 ~/rna_dohmain/homd_map/annotations.merge.txt) rpoB.annotations.txt > temp
mv temp rpoB.annotations.txt
parallel -a <(awk '{print $1}' rpoB.annotations.txt) -j 7 -k "grep '{}' ../../homd_map/read_counts.txt"> rpoB_counts.txt
cat <(head -n 1 ../../homd_map/read_counts.txt) rpoB_counts.txt > temp
mv temp rpoB_counts.txt
```
DESeq
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)

#load data
setwd("~/rna_dohmain/11-perio/01-rpoC-diff-abundance")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./rpoB_counts.txt", header=T, sep="\t", row.names=1)
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
# out of 8942 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 214, 2.4%
# LFC < 0 (down)     : 305, 3.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 5986, 67%
# (mean count < 1)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  519"
# out of 8942 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 273, 3.1%
# LFC < 0 (down)     : 359, 4%
# outliers [1]       : 0, 0%
# low counts [2]     : 5986, 67%
# (mean count < 1)
write.table(resLFC, file="deseq_results_rpoB-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_rpoB-HUUvHI.RData")
```
Valcona Plot
```R
# load("deseq_results_rpoC-HUUvHI.RData")
# add in annotations
homd <- read.table("rpoB.annotations.txt", header=T, sep="\t", quote="")
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
res_ord$GeneInfo <- paste(res_ord$Genus, res_ord$Species,res_ord$gene)

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 1,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-rpoB-HUUvHI.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-rpoB-HUUvHI.pdf")
```
# 2. DNA deseq2
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)

#load data
setwd("~/rna_dohmain/11-perio/01-rpoC-diff-abundance")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1)) # get rid of weird empty column in genecounts
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
# [1] "number of genes with adjusted p value lower than 0.05:  915"
# out of 8437 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 395, 4.7%
# LFC < 0 (down)     : 520, 6.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 6513, 77%
# (mean count < 1)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  859"
# out of 8437 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 612, 7.3%
# LFC < 0 (down)     : 631, 7.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 6038, 72%
# (mean count < 1)
write.table(resLFC, file="deseq_results_rdna-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_dna-HUUvHI.RData")
```
```R
load("deseq_results_dna-HUUvHI.RData")

# add in annotations
homd <- read.table("../../rpoc/taxonomy_bac.txt", header=F, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$V1 %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$V1
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
sigsp <- paste("x", sigloc$V7, sep="_")
sigdf <- as.data.frame(cbind(row.names(sigloc), sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# remove any lists (i.e., genera) with fewer than three hits
# siglist <- Filter(function(x) length(x) >=3, siglist)

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
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

#combine species and gene
res_ord$GeneInfo <- paste0(res_ord$V8, "_", row.names(res_ord))

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 1,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)

pdf("volcano_DNA-HUUvHI.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano_DNA-HUUvHI.red.pdf")

filtered_res <- subset(res_ord, V8 %in% c("Treponema_denticola", "Porphyromonas_gingivalis", "Tannerella_forsythia")) %>%  filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
```
# 2. DNA deseq2 species level
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(phyloseq)
#load data
setwd("~/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("../../rpoc/ps.RData")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
glom <- tax_glom(ps.dat, "V8")
genecounts <- t(otu_table(glom)) # get rid of weird empty column in genecounts
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
subcount <- as.data.frame(subcount)
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
# out of 578 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 39, 6.7%
# LFC < 0 (down)     : 103, 18%
# outliers [1]       : 0, 0%
# low counts [2]     : 202, 35%
# (mean count < 1)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  135"
# out of 578 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 61, 11%
# LFC < 0 (down)     : 111, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 135, 23%
# (mean count < 1)
write.table(resLFC, file="deseq_results_dna_species-HUUvHI.txt", quote=F, sep="\t")
save.image("deseq_results_dna_species-HUUvHI.RData")
```
```R
load("deseq_results_dna_species-HUUvHI.RData")

# add in annotations
homd <- read.table("../../rpoc/taxonomy_bac.txt", header=F, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$V1 %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$V1
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
sigsp <- paste("x", sigloc$V7, sep="_")
sigdf <- as.data.frame(cbind(row.names(sigloc), sigsp))

# split the dataframe as a list by genus
siglist <- split(sigdf$V1, sigdf$sigsp)

# remove any lists (i.e., genera) with fewer than three hits
# siglist <- Filter(function(x) length(x) >=3, siglist)

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
res_ord$Genus[is.na(res_ord$color)] <- NA
# change NA genus to grey
res_ord$color[is.na(res_ord$color)] <- "#808080"

# get key value pairs for plotting
colormap <- setNames(res_ord$color, res_ord$Genus)

#combine species and gene
res_ord$GeneInfo <- paste0(res_ord$V8, "_", row.names(res_ord))

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 1,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)

pdf("volcano_DNA_species-HUUvHI.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano_DNA_species-HUUvHI.red.pdf")

filtered_res <- subset(res_ord, V8 %in% c("Treponema_denticola", "Porphyromonas_gingivalis", "Tannerella_forsythia")) %>%  filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# only p. gingivalis
```






























# 2. Balance of taxa from DNA
```R
library(grid)
library(coda4microbiome)
library(tidyverse)
library(phyloseq)


setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
set.seed(545433543)
ps.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 5% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 10) > (0.5*length(x)), TRUE)
# pull data
dat <- as.data.frame(otu_table(glom))
map <- sample_data(glom)

#get taxonomy of ASV
taxa = as(tax_table(ps.dat), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")

#rename to have V8 level name
dat <- as.data.frame(dat)
dat <- dat %>% rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))
rownames(dat) <- dat$V8
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 3
y <- factor(dat[,dif2]) #geog lcaotion
# z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound
geo_its <- coda_glmnet(x=x,y=y)

sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.HUUvHI.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUUvHI.pdf")
```
Corncob
```R
library(phyloseq)
library(corncob)
library(magrittr)


setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
set.seed(12349)
ps.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 5% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 10) > (0.5*length(x)), TRUE
#choose a model
corncob <- bbdml(formula = ASV1 ~ 1,
phi.formula = ~ 1,
data = glom)

corncob_da <- bbdml(formula = ASV1 ~ hiv_status,
phi.formula = ~ hiv_status,
data = glom)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = glom,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.hiv_status.rpoc.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.rpoc.pdf")
```
Corncob HEU vs HI
```R
library(phyloseq)
library(corncob)
library(magrittr)


setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
set.seed(12349)
ps.sub <- subset_samples(ps.dat, hiv_status == "HEU" | hiv_status == "HI")

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 5% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 10) > (0.5*length(x)), TRUE
#choose a model
corncob <- bbdml(formula = ASV1 ~ 1,
phi.formula = ~ 1,
data = glom)

corncob_da <- bbdml(formula = ASV1 ~ hiv_status,
phi.formula = ~ hiv_status,
data = glom)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = glom,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.hiv_status.HEUvHI.rpoc.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.HEUvHI.rpoc.pdf")
```
Corncob HEU vs HI
```R
library(phyloseq)
library(corncob)
library(magrittr)


setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
set.seed(12349)
ps.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HEU")

# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.sub, "V8")

#diff abundance
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = glom,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.hiv_status.HUUvHEU.rpoc.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.HUUvHEU.rpoc.pdf")
```
ASV Level
```R
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
set.seed(12349)
ps.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")

# collapse data to roughly species level to minimize high sparsity
# glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 5% of samples post merging
# glom <- filter_taxa(glom, function(x) sum(x > 10) > (0.5*length(x)), TRUE
#choose a model
corncob <- bbdml(formula = ASV1 ~ 1,
phi.formula = ~ 1,
data = ps.sub)

corncob_da <- bbdml(formula = ASV1 ~ hiv_status,
phi.formula = ~ hiv_status,
data = ps.sub)

lrtest(mod_null = corncob, mod = corncob_da) #got a p-value of less than 0.05 -> want to use covariate model

#diff abundance
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.sub,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.hiv_status.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.asv.pdf")
```
See distro of P. gingivalis
```R
# Porphyromonas_gingivalis
rel <- microbiome::transform(glom, "compositional")
data <- psmelt(rel)
sub_data <- data[data$V8 == "Porphyromonas_gingivalis",]
pdf("pging_rpoc_dna.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=OTU),stat="identity")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./pging_rpoc_dna.pdf")
pdf("pging_dna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = Abundance)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pging_dna_reads.hist.pdf")
mean(sub_data$Abundance)
```
See distro of T. denticola
```R
# Porphyromonas_gingivalis
rel <- microbiome::transform(glom, "compositional")
data <- psmelt(rel)
sub_data <- data[data$V8 == "Treponema_denticola",]
pdf("tdent_rpoc_dna.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=OTU),stat="identity")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./tdent_rpoc_dna.pdf")
pdf("tdent_dna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = Abundance)) +
  geom_histogram(binwidth = 0.002, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./tdent_dna_reads.hist.pdf")
mean(sub_data$Abundance)
```
