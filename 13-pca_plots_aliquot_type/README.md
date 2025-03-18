## Ordination of RNA and DNA from different teeth

### 1. Load environment

```bash
cd /home/allie/domhain_RNAseq/08-pca_plots_aliquot_type
conda activate 2024-HIV_RNASeq
```

### 2. Load R libraries

```R
# remotes::install_github("gauravsk/ranacapa")
# remotes::install_github("Russel88/MicEco")
library(phyloseq, verbose=F)
library(ggplot2, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(ranacapa, verbose=F)
library(MicEco, verbose=F)
```

### 3. Load and clean up rpoC data

```R
# metadata & clean up
map <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
map$aliquot_type <- sub("-", "", map$aliquot_type)
row.names(map) <- map$sample_id
# sequence table
seqtab <- read.table("~/domhain_RNAseq/homd_rpoc_suzanne-7.2.24/sequence_table.merged.txt", header=T, sep="\t", row.names=1)
tax <- read.table("~/domhain_RNAseq/homd_rpoc_suzanne-7.2.24/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
```

Check that sample names match across sequencing and metadata file

```R
notinmeta <- setdiff(row.names(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), row.names(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw
# should both come back as character(0)  
```

Create initial phyloseq object

```R
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
```

Remove low abundance ASVs

```R
system("mkdir img")
# compute prevalence dataframe
prevdf <- apply(X=otu_table(ps.dat), MARGIN=ifelse(taxa_are_rows(ps.dat), yes=1, no=2), FUN=function(x){sum(x>0)})
# add taxa and total read counts to dataframe
prevdf <- data.frame(Prevalence=prevdf, TotalAbundance=taxa_sums(ps.dat), tax_table(ps.dat))
# which phyla are comprised as mostly low prevalence ASVs?
lowprev <- ggplot(prevdf, aes(TotalAbundance, Prevalence, nsamples(ps.dat), color="V4")) + geom_hline(yintercept=0.05, alpha=0.5, linetype=2) + geom_point(size=2, alpha=0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~V4) + theme(legend.position="none")
lowprev
pdf("img/totalabund_vs_prevalence.pdf")
lowprev
dev.off()
# ASVs must be found in minimum number of samples 2, minimum reads 50
ps.dat <- ps_prune(ps.dat, min.samples = 2, min.reads = 50)
ps.dat
# write out filtered sequence and taxonomy table
write.table(as.data.frame(otu_table(ps.dat)), "sequence_table.filt.txt", sep="\t", row.names=T, col.names=T, quote=F)
write.table(as.data.frame(tax_table(ps.dat)), "taxonomy_bac.filt.txt", sep="\t", row.names=T, col.names=T, quote=F)
```

### 4. Beta diversity rpoC

```R
# set up some color palettes for major groupings
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
# clr transformation
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")  
```

Distance based redundancy analysis

```R
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~hiv_status)
# capscale plot by HIV status group
pdf("img/bdiv_cap.hiv_status.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="hiv_status") + 
    theme_minimal() + 
    scale_color_manual(values=hivCols)
dev.off()
system("/home/allie/.iterm2/imgcat img/bdiv_cap.hiv_status.pdf")
# tooth health category
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~tooth_health)
pdf("img/bdiv_cap.tooth_health.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="tooth_health") + 
    theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat img/bdiv_cap.tooth_health.pdf")
# aliquot type
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~aliquot_type)
pdf("img/bdiv_cap.aliquot_type.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="aliquot_type") + 
    theme_minimal() +
    scale_color_manual(values=healthCols)
dev.off()
system("/home/allie/.iterm2/imgcat img/bdiv_cap.aliquot_type.pdf")
```

### 5. DESeq2 Analysis of all genes comparing health categories

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")

library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# Load data 
metadata <- read.table("~/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("~/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
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
colnames(genecounts)
rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true

# create deseq object -- first comparing H vs D
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
# [1] "number of genes with adjusted p value lower than 0.05:  595"

# out of 4208 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 537, 13%
# LFC < 0 (down)     : 58, 1.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 2396, 57%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# health is positive, dentin cavity negative
resLFC <- lfcShrink(se_star, coef="tooth_health_H_vs_D", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  600"

# out of 4208 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 596, 14%
# LFC < 0 (down)     : 119, 2.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 2497, 59%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-HvD.txt", quote=F, sep="\t")
# save.image()
```

### 3. PCA plot comparing H to D all ADS activity


```R
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-HvD.pdf")
plotPCA(vld, intgroup=c("tooth_health")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-HvD.pdf")
# aliquot type?
pdf("pca_pdvpf_ADS-aliquot_type.pdf")
plotPCA(vld, intgroup=c("aliquot_type")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-aliquot_type.pdf")
# hiv group?
pdf("pca_pdvpf_ADS-hiv_status.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-hiv_status.pdf")
# no realy clustering patterns
```
