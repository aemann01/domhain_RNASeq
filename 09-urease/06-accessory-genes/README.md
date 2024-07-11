# 1. Get seq ids for ureABC
```sh
cd ~/rna_dohmain/09-urease/06-accessory-genes
grep "ureD" ../../homd_map/ALL_genomes.gff| grep Urease | sed 's/|.*//' > ureD.seqs
grep "ureE" ../../homd_map/ALL_genomes.gff| grep Urease| sed 's/|.*//' > ureE.seqs
grep "ureF" ../../homd_map/ALL_genomes.gff| grep Urease| sed 's/|.*//' > ureF.seqs
grep "ureG" ../../homd_map/ALL_genomes.gff| grep Urease | sed 's/|.*//' > ureG.seqs
#these are less common
grep "ureH" ../../homd_map/ALL_genomes.gff| grep Urease| sed 's/|.*//' > ureH.seqs
grep "ureR" ../../homd_map/ALL_genomes.gff| grep Urease| sed 's/|.*//' > ureR.seqs
wc -l *seqs
```
# 2. See which sequences have ureDEFG
```sh 
python3 common_elements.py #ureDEFG.seqs
```
# 3. Get frame with upregulated seqids
```sh
cat ureDEFG.seqs | while read line; do grep $line ../../homd_map/annotations.merge.txt | grep ureD; done > ureD.seqid
cat ureDEFG.seqs | while read line; do grep $line ../../homd_map/annotations.merge.txt | grep ureE; done > ureE.seqid
cat ureDEFG.seqs | while read line; do grep $line ../../homd_map/annotations.merge.txt | grep ureF; done > ureF.seqid
cat ureDEFG.seqs | while read line; do grep $line ../../homd_map/annotations.merge.txt | grep ureG; done > ureG.seqid
# get all seqids together
cat ureD.seqid ureE.seqid ureF.seqid ureG.seqid | awk '{print $1}' >seqids
parallel -a seqids -j 7 -k "grep -wm 1 '{}' ../../homd_map/read_counts.txt"> ureDEFG_counts.txt
head -n 1 ../../homd_map/read_counts.txt > ./temp
cat ./temp ./ureDEFG_counts.txt > ./temp1
mv ./temp1 ./ureDEFG_counts.txt
#make annotation file
cat ureD.seqid ureE.seqid ureF.seqid ureG.seqid > ./ureDEFG.annotations.txt
head -n 1 ../../homd_map/annotations.merge.txt > ./temp
cat ./temp ./ureDEFG.annotations.txt > ./temp1
mv ./temp1 ./ureDEFG.annotations.txt
```
# 4. Run DESeq for ureDEFG accesory genes
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)

#load data
setwd("/home/suzanne/rna_dohmain/09-urease/06-accessory-genes")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./ureDEFG_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$tooth_health == "H" | metadata$tooth_health == "D",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  118"
# out of 7938 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 82, 1%
# LFC < 0 (down)     : 36, 0.45%
# outliers [1]       : 0, 0%
# low counts [2]     : 7521, 95%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# # [1] "number of genes with adjusted p value lower than 0.05:  118"
# out of 7938 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 72, 0.91%
# LFC < 0 (down)     : 68, 0.86%
# outliers [1]       : 0, 0%
# low counts [2]     : 7521, 95%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?result
write.table(resLFC, file="deseq_results_ure_acces-HvD.txt", quote=F, sep="\t")
save.image("deseq_results_ure_acces-HvD.RData")
```
Valcona Plot
```R
load("deseq_results_ure_acces-HvD.RData")
# add in annotations
homd <- read.table("./ureDEFG.annotations.txt", header=T, sep="\t", quote="")
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
	labSize = 1,
	shape = 19,
	legendPosition = 'right',
	boxedLabels = TRUE,
	drawConnectors = TRUE,
	pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
	colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
)
pdf("volcano-HvD.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HvD.pdf")

#see which ones have all 4 genes diff expressed
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
ureDEFG_name <- names(which(table(sig_df$SEQ_ID) >= 4))
names(which(table(sig_df$SEQ_ID) >= 3))
```