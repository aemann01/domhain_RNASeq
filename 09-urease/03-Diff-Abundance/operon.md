# 1. Get operon read counts from palmetto
```sh
cd ~/rna_dohmain/10-urease/03-Diff-Abundance
rsync -a scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/filtered/featurecounts/read_counts.txt ./
egrep -h "ureA|ureB|ureC" ../../homd_map/annotations.merge.txt > ../../homd_map/ure.annotations.txt
head -n 1 ../../homd_map/annotations.merge.txt > ../../temp
cat ../../temp ../../ure.annotations.txt > ../../temp1
mv ../../temp1 ../../ure.annotations.txt
```
# 2. Run DESeq for overall
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)

#load data
setwd("/home/suzanne/rna_dohmain/09-urease/03-diff-abundance")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.operon", replacement = "") 
colnames(genecounts) <- gsub(x = colnames(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
submap <- arrange(submap, sample_id)

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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "H"))

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
# # [1] "number of genes with adjusted p value lower than 0.05:  79"
# out of 6705 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 49, 0.73%
# LFC < 0 (down)     : 30, 0.45%
# outliers [1]       : 0, 0%
# low counts [2]     : 6340, 95%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# # [1] "number of genes with adjusted p value lower than 0.05:  79"
# out of 6705 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 64, 0.95%
# LFC < 0 (down)     : 46, 0.69%
# outliers [1]       : 0, 0%
# low counts [2]     : 6340, 95%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_operon-HvD.txt", quote=F, sep="\t")
save.image("deseq_results_operon-HvD.RData")
```
Valcona Plot
```R
load("deseq_results_operon-HvD.RData")
# add in annotations
homd <- read.table("~/rna_dohmain/homd_map/ure.annotations.txt", header=T, sep="\t", quote="")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene,res_ord$SEQ_ID)

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
	labSize = 2,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-operon-HvD.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-operon-HvD.pdf")
#see which ones have all 3 genes diff expressed
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
ureABC_name <- names(which(table(sig_df$SEQ_ID) >= 3))
names(which(table(sig_df$SEQ_ID) >= 3))
# character(0)
```
Make beta diversity plot
```R
vld <- varianceStabilizingTransformation(se_star)
pdf("pca_HvD_allureABC.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca_HvD_allureABC.pdf")
```
Make dendogram of ureABC activity
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
mat_colors <- list(hiv_status = c("#8213A0", "#FA78FA","#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HEU", "HUU")
mat_colors$tooth_health <- c("#FA918B", "#44CCD0")
names(mat_colors$tooth_health) <- c("D", "H")
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("~/rna_dohmain/homd_map/ure.annotations.txt", header=T, sep="\t", quote="")
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
Make a donut plot
```R
library(reshape2)
library(webr)
resm <- as.matrix(genecounts)
long_res <- melt(resm)
#get metadata
homd <- read.table("~/rna_dohmain/homd_map/ure.annotations.txt", header=T, sep="\t", quote="")
map <- metadata[c("sample_id","hiv_status","tooth_health")]
cords <- read.table("../02-operon-mapping/new.cords2", header=T, sep="\t")

#combine files
homd$Species <- paste0(homd$Genus, "_", homd$Species)
res_seqs <- inner_join(long_res, homd, by = join_by(Var1 == locus_tag))
res_seqs <- inner_join(res_seqs, map, by = join_by(Var2 == sample_id))

#make donut by ureABC
res_seqs$Gene <- gsub(x = res_seqs$gene, pattern = "_[0-9]", replacement = "")
res_seqs$Gene <- gsub(x = res_seqs$Gene, pattern = "[0-9]", replacement = "")
#normalize by operon length

#ureA
res_A <- res_seqs[res_seqs$Gene == "ureA",]
cord_A<- as.data.frame(cords$idA)
cord_A$length <- cords$A2 -cords$A1
cord_A$norm <-  cord_A$length/mean(cord_A$length)
unique(sort(cord_A$length/cord_A$norm))
colnames(cord_A) <- c('id', 'length', 'norm')
res_A <-inner_join(res_A, cord_A, by = join_by(Var1 ==  id))
res_A$norm_Gene<- res_A$value/res_A$norm
#ureB
res_B <- res_seqs[res_seqs$Gene == "ureB", ]
cord_B<- as.data.frame(cords$idB)
colnames(cord_B) <- c('id')
cord_B$length <- cords$B2 -cords$B1
cord_B$norm <-  cord_B$length/mean(cord_A$length) #normalizing by length for A
unique(sort(cord_B$length/cord_B$norm))
res_B <-inner_join(res_B, cord_B, by = join_by(Var1 ==  id))
res_B$norm_Gene<- res_B$value/res_B$norm
#ureC
res_C <- res_seqs[res_seqs$Gene == "ureC",]
cord_C<- as.data.frame(cords$idC)
colnames(cord_C) <- c('id')
cord_C$length <- cords$C2 -cords$C1
cord_C$norm <-  cord_C$length/mean(cord_A$length) #normalizing by length for A
unique(sort(cord_C$length/cord_C$norm))
res_C <-inner_join(res_C, cord_C, by = join_by(Var1 ==  id))
res_C$norm_Gene<- res_C$value/res_C$norm
#combine
res_seqs <- rbind(res_A, res_B, res_C)
#make donut plots
ure_count <- res_seqs %>% 
  group_by(tooth_health, Gene) %>% 
  summarise(Total = sum(norm_Gene, na.rm = TRUE))

pdf("gene_ure_donut.pdf")
PieDonut(ure_count, aes("tooth_health", "Gene", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./gene_ure_donut.pdf")

ure_count <- res_seqs %>% 
  group_by(tooth_health, Species) %>% 
  summarise(Total = sum(norm_Gene, na.rm = TRUE))
pdf("species_ure_donut.pdf", width =30, height =30)
PieDonut(ure_count, aes("tooth_health", "Species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./species_ure_donut.pdf")

ure_count <- res_seqs %>% 
  group_by(hiv_status, Species) %>% 
  summarise(Total = sum(norm_Gene, na.rm = TRUE))
pdf("hiv_species_ure_donut.pdf")
PieDonut(ure_count, aes("hiv_status", "Species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./hiv_species_ure_donut.pdf")

ure_count <- res_seqs %>% 
  group_by(hiv_status, Gene) %>% 
  summarise(Total = sum(norm_Gene, na.rm = TRUE))
pdf("hiv_ure_donut.pdf")
PieDonut(ure_count, aes("hiv_status", "Gene", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./hiv_ure_donut.pdf")

ure_count <- res_seqs %>% 
  group_by(Var2, Species) %>% 
  summarise(Total = sum(norm_Gene, na.rm = TRUE))
pdf("sample_ure_donut.pdf", width =50, height =50)
PieDonut(ure_count, aes("Var2", "Species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./sample_ure_donut.pdf")
#find prportion of transcripts that are A. neaslundii
data_sorted <- ure_count[order(ure_count$Total,
                                decreasing = TRUE), ]
 
# select top 1 values from each sample
data_mod <- Reduce(rbind,                                
                    by(data_sorted,
                       data_sorted["Var2"],
                       head,
                       n = 1))
ure_count$a_neas <- ifelse(ure_count$Species == "Actinomyces_naeslundii",  "yes", "no")
ure_count <- ure_count %>% group_by(Var2, a_neas) %>% summarise(reads = sum(Total, na.rm = TRUE))
total_reads <- ure_count %>% group_by(Var2) %>% summarise(sum = sum(reads, na.rm = TRUE))
a_count <- ure_count[ure_count$a_neas == 'yes',]
reads <-inner_join(total_reads, a_count, by = join_by(Var2 ==  Var2))
reads$pro <- reads$reads/reads$sum
```
# 3. HUU DESeq
```R
# filter metadata so that we only compare H to D
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
submap <- submap[submap$hiv_status == "HUU",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "H"))

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
# # [1] "number of genes with adjusted p value lower than 0.05:  27"
# out of 214 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 27, 13%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 88, 41%
# low counts [2]     : 53, 25%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# # [1] "number of genes with adjusted p value lower than 0.05:  25"
# out of 214 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 39, 18%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 88, 41%
# low counts [2]     : 42, 20%
write.table(resLFC, file="deseq_results_operon-HvD-HUU.txt", quote=F, sep="\t")
save.image("deseq_results_operon-HvD-HUU.RData")
```
HUU Valcona Plot
```R
# filter by locus tag 
load("deseq_results_ure-HvD-HUU.RData")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene,res_ord$SEQ_ID)

#Create volcano plot
HUU_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 2,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HvD-operon-HUU.pdf", width=15, height=10)
HUU_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HvD-operon-HUU.pdf")
#see which ones have all 3 genes diff expressed
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
ureABC_name <- names(which(table(sig_df$SEQ_ID) >= 3))
names(which(table(sig_df$SEQ_ID) >= 3))
# [1] "SEQF5153.1" "SEQF5162.1"
```
# 4. HEU DESeq
```R
# filter metadata so that we only compare H to D
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
submap <- submap[submap$hiv_status == "HEU",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "H"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  27"
# out of 243 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 14, 5.8%
# LFC < 0 (down)     : 13, 5.3%
# outliers [1]       : 83, 34%
# low counts [2]     : 57, 23%

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  23"
# out of 243 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 18, 7.4%
# LFC < 0 (down)     : 11, 4.5%
# outliers [1]       : 83, 34%
# low counts [2]     : 0, 0%
write.table(resLFC, file="deseq_results_operon-HvD-HEU.txt", quote=F, sep="\t")
save.image("deseq_results_operon-HvD-HEU.RData")
```
HEU Valcona Plot
```R
# filter by locus tag 
load("deseq_results_operon-HvD-HEU.RData")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene, res_ord$SEQ_ID)

#Create volcano plot
HEU_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 3,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HvD-operon-HEU.pdf", width=15, height=10)
HEU_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HvD-operon-HEU.pdf")
#see which ones have all 3 genes diff expressed
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
ureABC_name <- names(which(table(sig_df$SEQ_ID) >= 3))
names(which(table(sig_df$SEQ_ID) >= 3))
# character(0)
```
# 6. HI DESeq
```R
# filter metadata so that we only compare H to D
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
submap <- submap[submap$hiv_status == "HI",]
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
star_results$tooth_health <- factor(star_results$tooth_health, levels=c("D", "H"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  0"
# out of 187 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 53, 28%
# low counts [2]     : 0, 0%

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# # [1] "number of genes with adjusted p value lower than 0.05:  2"
# # out of 187 with nonzero total read count
# LFC > 0 (up)       : 12, 6.4%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 53, 28%
# low counts [2]     : 114, 61%
write.table(resLFC, file="deseq_results_operon-HvD-HI.txt", quote=F, sep="\t")
save.image("deseq_results_operon-HvD-HI.RData")
```
Valcona Plot
```R
# filter by locus tag 
load("deseq_results_operon-HvD-HI.RData")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene,res_ord$SEQ_ID)

#Create volcano plot
HI_plot <- EnhancedVolcano(res_ord,
	lab = res_ord$GeneInfo,
	x = 'log2FoldChange',
	y = 'padj',
	FCcutoff = lfc,
	pCutoff = pval,
	colCustom = colormap ,
	title = "",
	subtitle = "",
	caption = "",
	labSize = 3,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HvD-operon-HI.pdf", width=15, height=10)
HI_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HvD-operon-HI.pdf")
#see which ones have all 3 genes diff expressed
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
ureABC_name <- names(which(table(sig_df$SEQ_ID) >= 3))
names(which(table(sig_df$SEQ_ID) >= 3))
# character(0)
```
# 7. Combine plots
```R
library(cowplot)
library("gridExtra")
library(ggpubr)
#combine plots
overall_plot <- overall_plot + theme(legend.position = c(0.8, 0.8))
HUU_plot <- HUU_plot + theme(legend.position = c(0.8, 0.8))
HEU_plot <- HEU_plot + theme(legend.position = c(0.8, 0.8))
HI_plot  <- HI_plot + theme(legend.position = c(0.8, 0.8))

gt <- arrangeGrob(overall_plot, 
             HUU_plot, HEU_plot, HI_plot
             layout_matrix = rbind(c(1,1,1), c(2,3,4)))
pdf("combine.pdf")
as_ggplot(gt)# Add labels
dev.off()
system("~/.iterm2/imgcat ./combine.pdf")
```