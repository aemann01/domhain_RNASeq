Global DESeq2 Analysis of HIV status groups

First comparing HEU to HUU

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")

library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
setwd("/home/allie/domhain_RNAseq/05-global_diffabund")
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("/home/allie/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)

# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare HEU to HUU
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HUU",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  686804"

# out of 22973393 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 208960, 0.91%
# LFC < 0 (down)     : 477844, 2.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 19487428, 85%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# HEU is positive, HUU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  686804"

# out of 22973393 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 295349, 1.3%
# LFC < 0 (down)     : 593489, 2.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 19487428, 85%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_global-HUUvHEU.txt", quote=F, sep="\t")
save.image("deseq_global.RData")
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

Heatmap of global RNA DE expression 

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")
# install.packages("pheatmap")
library(pheatmap, warn.conflicts = F, quietly = T)





df <- as.data.frame(colData(vld)[,c("ads_nmol_min_mg_of_protein","tooth_type")])
df$ads_nmol_min_mg_of_protein <- log2(as.numeric(df$ads_nmol_min_mg_of_protein))
colnames(df) <- c("ads", "tooth_type")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
library(RColorBrewer)
x <- pheatmap(top25Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_pdvpf.pdf")












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

# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low10 <- head(sortdf$locus_tag, 10)
# negative top 10
top10 <- tail(sortdf$locus_tag, 10)
# concatenate
labgenes <- c(top10, low10)

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