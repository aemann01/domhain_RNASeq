# 1. Run DESeq HI vs HUU
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
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
# [1] "number of genes with adjusted p value lower than 0.05:  11501"

# out of 51006 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 118, 0.23%
# LFC < 0 (down)     : 11383, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 31599, 62%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  11169"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 480, 0.94%
# LFC < 0 (down)     : 12333, 24%
# outliers [1]       : 0, 0%
# low counts [2]     : 28918, 57%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_red-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_red-HIvHUU.RData")
```
Valcona Plot
```R
load("deseq_results_red-HIvHUU.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 
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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHUU.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.red.pdf")

proteases = c("hagA", "prtP", "rgpA", "rgpB", "kgp")
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$gene %in% proteases, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHUU.vince.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.vince.pdf")

#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HIvHUU.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.red_viru.pdf")


# For Vince

res_sub <- res_ord[res_ord$species == "Porphyromonas_gingivalis",]
#sub color
colorsub <- rep("#808080", nrow(res_sub))
names(colorsub) <- rep("nonsignificant", length(colorsub))

# Green: log2FC >= 2 and padj < 0.05
green_idx <- which((res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) & res_sub$padj < pval)

colorsub[green_idx] <- "#0571B0"
names(colorsub)[green_idx] <- "significant"

# Red: gene name is 'rgpB'
red_idx <- which(
  res_sub$gene_name %in% c("rgpB", "rgpA", "kgp", "vimF", "porS",
                           "fimA", "fimC", "fimE", "hagA", "hagC", "mfa1") &
  ( (res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) &
    res_sub$padj < pval )
)
colorsub[red_idx] <- "#CA0020"
names(colorsub)[red_idx] <- "virulence"


overall_plot <- EnhancedVolcano(res_sub,
  lab = ifelse(res_sub$gene %in% "bleep", paste(res_sub$SEQ, res_sub$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colorsub ,
  title = "Porphyromonas gingivalis",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75))
) 
pdf("volcano-HIvHUU.pging.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.pging.pdf")

res_sub <- res_ord[res_ord$species == "Treponema_denticola",]
#sub color
colorsub <- rep("#808080", nrow(res_sub))
names(colorsub) <- rep("nonsignificant", length(colorsub))

# Green: log2FC >= 2 and padj < 0.05
green_idx <- which((res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) & res_sub$padj < pval)

colorsub[green_idx] <- "#0571B0"
names(colorsub)[green_idx] <- "significant"

# Red: gene name is 'rgpB'
red_idx <- which(
  res_sub$gene_name %in% c("msp", "prtP", "oppA", "fhbB", "lrrA",
                           "flaA", "fliL", "fliM", "fliY", "flgB", "flgC", "flgD", "flgE", "cheA", "cheX", "cheY") &
  ( (res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) &
    res_sub$padj < pval )
)
colorsub[red_idx] <- "#CA0020"
names(colorsub)[red_idx] <- "virulence"


overall_plot <- EnhancedVolcano(res_sub,
  lab = ifelse(res_sub$gene %in% "bleep", paste(res_sub$SEQ, res_sub$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colorsub ,
  title = "Treponema denticola",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75))
) 
pdf("volcano-HIvHUU.tdent.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.tdent.pdf")


res_sub <- res_ord[res_ord$species == "Tannerella_forsythia",]
#sub color
colorsub <- rep("#808080", nrow(res_sub))
names(colorsub) <- rep("nonsignificant", length(colorsub))

# Green: log2FC >= 2 and padj < 0.05
green_idx <- which((res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) & res_sub$padj < pval)

colorsub[green_idx] <- "#0571B0"
names(colorsub)[green_idx] <- "significant"

# Red: gene name is 'rgpB'
red_idx <- which(
  res_sub$gene_name %in% c("bspA", "tfsA", "tfsB", "wecC") &
  ( (res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) &
    res_sub$padj < pval )
)
colorsub[red_idx] <- "#CA0020"
names(colorsub)[red_idx] <- "virulence"


overall_plot <- EnhancedVolcano(res_sub,
  lab = ifelse(res_sub$gene %in% "bleep", paste(res_sub$SEQ, res_sub$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colorsub ,
  title = "Tannerella forsythia",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75))
) 
pdf("volcano-HIvHUU.tfor.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.tfor.pdf")


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
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
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
# [1] "number of genes with adjusted p value lower than 0.05:  8566"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 141, 0.28%
# LFC < 0 (down)     : 8425, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 35409, 69%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# HUU is positive, HEU negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  8566"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 776, 1.5%
# LFC < 0 (down)     : 9474, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 35409, 69%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.table(resLFC, file="deseq_results_red-HEUvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_red-HEUvHUU.RData")
```
Valcona Plot
```R
load("deseq_results_red-HEUvHUU.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HEUvHUU.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.red.pdf")
#Create volcano plot

res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

sigloc %>% filter(gene =="prtP")

overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HEUvHUU.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.red_viru.pdf")
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
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
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
# [1] "number of genes with adjusted p value lower than 0.05:  5308"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1316, 2.6%
# LFC < 0 (down)     : 3992, 7.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 31624, 62%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# HEU is positive, HI negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HEU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  5105"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1666, 3.3%
# LFC < 0 (down)     : 4969, 9.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 29616, 58%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.table(resLFC, file="deseq_results_red-HIvHEU.txt", quote=F, sep="\t")
save.image("deseq_results_red-HIvHEU.RData")
```
Valcona Plot
```R
load("deseq_results_red-HIvHEU.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 
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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHEU.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHEU.red.pdf")
#Create volcano plot

res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

sigloc %>% filter(gene =="prtP")

overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HIvHEU.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHEU.red_viru.pdf")
```
# 4. Heatmap
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(microbiome)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU","HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
#pcoa diversity
vld <- varianceStabilizingTransformation(se_star)
pdf("pca.hiv.red.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca.hiv.red.pdf")

#make dendogram
#get virulence factors ids
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
#Get top varying genes
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 50)
 
#make a subset of the log transformed counts for just the top 25 varying genes
topCounts <- assay(vld)[topVarGenes,]

df <- as.data.frame(colData(vld)[,c("hiv_status","age_y")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(hiv_status = c("#8213A0","#FA78FA", "#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HEU","HUU")
#get relative abundance of red complex
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
abund <- t(as.data.frame(otu_table(sub)))
abund_df <- as.data.frame(abund)
abund_with_sum <- abund_df %>%
  rowwise() %>%
  mutate(Sum = sum(c_across(everything()))) %>%
  ungroup()
abund <-as.data.frame(abund_with_sum$Sum)
row.names(abund) <- row.names(t(as.data.frame(otu_table(sub))))

merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "age", "rel_abundance")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$rel_abundance <- as.numeric(merged_df$rel_abundance)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  tag))
rows2$name <- paste0(rows2$species, "_", rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name

x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),fontsize_row= 4, fontsize_col =4)
save_pheatmap_pdf(x, "heatmap.hiv.red.pdf")
system("~/.iterm2/imgcat ./heatmap.hiv.red.pdf")
```
# 5. Upset Plot
```R
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ggplot2)
library(ComplexUpset)

setwd("~/rna_dohmain/11-perio/04-red-diff")
HI <- read.csv("deseq_results_red-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HI),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HI), rownames(ann)))]
HI <- HI[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HI)==rownames(ann)) # should all return true
# if all are true, merge together
HI <- cbind(HI, ann)
#HEU
HEU$hiv_status <- "HEU"
#get annots for HEEU
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HEU),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HEU), rownames(ann)))]
HEU <- HEU[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HEU)==rownames(ann)) # should all return true
# if all are true, merge together
HEU <- cbind(HEU, ann)

#bind together HI and HEU
resdf <- rbind(HI, HEU)

#get only signficant genes
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of -2 since that is the side from their respective deseq analysis
resdf <- resdf %>%
  mutate(mark = if_else(log2FoldChange <= -lfc & padj <= pval, "1", "0"))
resdf$mark <- as.numeric(resdf$mark)
resdf_filtered <- select(resdf, mark, hiv_status, gene_tag, species)
resdf_wide <- resdf_filtered %>%
  tidyr::pivot_wider(names_from = hiv_status, values_from = mark, values_fill = list(mark = 0)) %>%
  arrange(gene_tag)
#make upset plot
resdf_df <- as.data.frame(resdf_wide)
resdf_clean_df <- resdf_df %>%
  filter(!is.na(HI) & !is.na(HEU))
resdf_clean_df$species <- as.factor(resdf_clean_df$species)

pdf("HI_HEU.upset.pdf")
upset(resdf_clean_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially Expressed Genes", sets.x.label="Total Number of Differentially Expressed Genes",
  queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()
system("~/.iterm2/imgcat ./HI_HEU.upset.pdf")
```
# 6. Donut Plot
Reads belonging to each group
```R
library(tidyverse)
library(reshape2)
library(webr)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
#get annots for genes
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(genecounts),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(genecounts), rownames(ann)))]
genecounts <- genecounts[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(genecounts)==rownames(ann)) # should all return true
# if all are true, merge together
genecounts <- cbind(genecounts, ann)
gene_long <- melt(genecounts)
combined_data <- left_join(gene_long, metadata, by = c("variable" = "sample_id"))

#make donut plots
gene_sum <- combined_data %>% 
  group_by(hiv_status, genome, species) %>% 
  summarise(Total_genome = sum(value, na.rm = TRUE))
gene_sum2 <- gene_sum %>% 
  group_by(hiv_status, species) %>% 
  summarise(Total = mean(Total_genome, na.rm = TRUE))
pdf("red_rna.donut.pdf")
PieDonut(gene_sum2, aes("hiv_status", "species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./red_rna.donut.pdf")
```
DNA donut
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(microbiome)
setwd("~/rna_dohmain/11-perio/06-red-complex")
#get relative abundance of dna
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, Abundance, V8)
#make donut plots
gene_sum <- red_dna %>% 
  group_by(hiv_status, V8) %>% 
  summarise(Total = sum(Abundance, na.rm = TRUE))
pdf("red_dna.donut.pdf")
PieDonut(gene_sum, aes("hiv_status", "V8", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./red_dna.donut.pdf")  
```
# 7. Run DESeq HI vs HUU using global diff and subset to red
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
load("../03-global-diff/deseq_results-HIvHUU.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "", fill = TRUE) 
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
resdf <- resdf[resdf$species == "Porphyromonas_gingivalis" | resdf$species == "Treponema_denticola" | resdf$species == "Tannerella_forsythia",]

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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHUU.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.pdf")
#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HIvHUU.viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.viru.pdf")


res_sub <- res_ord[res_ord$species == "Porphyromonas_gingivalis",]
#sub color
colorsub <- rep("#808080", nrow(res_sub))
names(colorsub) <- rep("non-sig", length(colorsub))

# Green: log2FC >= 2 and padj < 0.05
green_idx <- which((res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) & res_sub$padj < pval)

colorsub[green_idx] <- "#00FF00"
names(colorsub)[green_idx] <- "upregulated"

# Red: gene name is 'rgpB'
red_idx <- which(
  res_sub$gene_name %in% c("rgpB", "rgpA", "kgp", "vimF", "porS",
                           "fimA", "fimC", "fimE", "hagA", "hagC", "mfa1") &
  ( (res_sub$log2FoldChange <= -lfc | res_sub$log2FoldChange >= lfc) &
    res_sub$padj < pval )
)
colorsub[red_idx] <- "#FF0000"
names(colorsub)[red_idx] <- "virulence"


overall_plot <- EnhancedVolcano(res_sub,
  lab = ifelse(res_sub$gene %in% "bleep", paste(res_sub$SEQ, res_sub$gene, sep=" "), ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colorsub ,
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
pdf("volcano-HIvHUU.pging.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.pging.pdf")

```
# 8. Run DESeq HEUI vs HUU using global diff and subset to red
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
load("../03-global-diff/deseq_results-HEUvHUU.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "", fill = TRUE) 
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
resdf <- resdf[resdf$species == "Porphyromonas_gingivalis" | resdf$species == "Treponema_denticola" | resdf$species == "Tannerella_forsythia",]

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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HEUvHUU.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.pdf")
#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HEUvHUU.viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.viru.pdf")
```
# 9. Upset Plot using global diff and subset to red
```R
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ggplot2)
library(ComplexUpset)

setwd("~/rna_dohmain/11-perio/04-red-diff")
HI <- read.csv("../03-global-diff/deseq_results-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("../03-global-diff/deseq_results-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HI),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HI), rownames(ann)))]
HI <- HI[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HI)==rownames(ann)) # should all return true
# if all are true, merge together
HI <- cbind(HI, ann)
HI <- HI[HI$species == "Porphyromonas_gingivalis" | HI$species == "Treponema_denticola" | HI$species == "Tannerella_forsythia",]
#HEU
HEU$hiv_status <- "HEU"
#get annots for HEEU
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HEU),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HEU), rownames(ann)))]
HEU <- HEU[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HEU)==rownames(ann)) # should all return true
# if all are true, merge together
HEU <- cbind(HEU, ann)
HEU <- HEU[HEU$species == "Porphyromonas_gingivalis" | HEU$species == "Treponema_denticola" | HEU$species == "Tannerella_forsythia",]

#bind together HI and HEU
resdf <- rbind(HI, HEU)

#get only signficant genes
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of -2 since that is the side from their respective deseq analysis
resdf <- resdf %>%
  mutate(mark = if_else(log2FoldChange <= -lfc & padj <= pval, "1", "0"))
resdf$mark <- as.numeric(resdf$mark)
resdf_filtered <- select(resdf, mark, hiv_status, gene_tag, species)
resdf_wide <- resdf_filtered %>%
  tidyr::pivot_wider(names_from = hiv_status, values_from = mark, values_fill = list(mark = 0)) %>%
  arrange(gene_tag)
#make upset plot
resdf_df <- as.data.frame(resdf_wide)
resdf_clean_df <- resdf_df %>%
  filter(!is.na(HI) & !is.na(HEU))
resdf_clean_df$species <- as.factor(resdf_clean_df$species)

pdf("HI_HEU.upset.pdf")
upset(resdf_clean_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially Expressed Genes", sets.x.label="Total Number of Differentially Expressed Genes",
  queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()
system("~/.iterm2/imgcat ./HI_HEU.upset.pdf")
```
# 10. Run DESeq HI H v D
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI",]
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$tooth_health, levels=c("H", "D"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  858"
# out of 22096 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 15, 0.068%
# LFC < 0 (down)     : 843, 3.8%
# outliers [1]       : 6043, 27%
# low counts [2]     : 7008, 32%
# (mean count < 2)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  858"
# out of 22096 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 90, 0.41%
# LFC < 0 (down)     : 1109, 5%
# outliers [1]       : 6043, 27%
# low counts [2]     : 7008, 32%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_red-HI.HvD.txt", quote=F, sep="\t")
save.image("deseq_results_red-HI.HvD.RData")
```
Valcona Plot
```R
load("deseq_results_red-HI.HvD.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 
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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HI.HvD.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HI.HvD.red.pdf")

#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HI.HvD.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HI.HvD.red_viru.pdf")
```
# 12. Run DESeq HEU HvD
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU",]
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$tooth_health, levels=c("H", "D"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  172"
# out of 15649 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.032%
# LFC < 0 (down)     : 167, 1.1%
# outliers [1]       : 4707, 30%
# low counts [2]     : 5758, 37%
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  168"
# out of 15649 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 34, 0.22%
# LFC < 0 (down)     : 222, 1.4%
# outliers [1]       : 4707, 30%
# low counts [2]     : 5561, 36%
# (mean count < 2)
write.table(resLFC, file="deseq_results_red-HEU.HvD.txt", quote=F, sep="\t")
save.image("deseq_results_red-HEU.HvD.RData")
```
Valcona Plot
```R
load("deseq_results_red-HEU.HvD.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 
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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HEU.HvD.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEU.HvD.red.pdf")

#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HEU.HvD.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEU.HvD.red_viru.pdf")
```
# 13. Run DESeq HUU HvD
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/04-red-diff")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU",]
submap <- submap[submap$tooth_health == "H" | submap$tooth_health == "D",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~tooth_health)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$tooth_health, levels=c("H", "D"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  100"
# out of 3998 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 8, 0.2%
# LFC < 0 (down)     : 92, 2.3%
# outliers [1]       : 1810, 45%
# low counts [2]     : 1627, 41%
# (mean count < 2)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  88"
# out of 3998 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 22, 0.55%
# LFC < 0 (down)     : 108, 2.7%
# outliers [1]       : 1810, 45%
# low counts [2]     : 1473, 37%
# (mean count < 2)
write.table(resLFC, file="deseq_results_red-HUU.HvD.txt", quote=F, sep="\t")
save.image("deseq_results_red-HUU.HvD.RData")
```
Valcona Plot
```R
load("deseq_results_red-HUU.HvD.RData")
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$SEQ <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 
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
res_ord$GeneInfo <- paste(res_ord$SEQ,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$SEQ, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HUU.HvD.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUU.HvD.red.pdf")

#Create volcano plot
res_sub <- res_ord %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
color_viru <- setNames(res_sub$color, res_sub$genus)

res_ord %>% filter(gene =="bspA")


overall_plot <- EnhancedVolcano(res_sub,
  lab = res_sub$GeneInfo,
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = color_viru ,
  title = "",
  subtitle = "",
  caption = "",
  labSize = 1,
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  pointSize = (ifelse(rownames(res_sub) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_sub) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HUU.HvD.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HUU.HvD.red_viru.pdf")
```