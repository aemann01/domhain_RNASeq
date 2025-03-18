# Global differential abundance 

Download and install mamba

```bash
# wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
# bash Miniforge3-Linux-x86_64.sh
```

Install R and other packges with mamba (running on Hillary)

```bash
# conda create --name 2024-HIV_RNASeq python=3.8
conda activate 2024-HIV_RNASeq
# mamba install -c conda-forge r-base
# mamba install -c bioconda subread
# mamba install -c conda-forge r-magrittr
# mamba install -c conda-forge r-corncob
```

Install required R libraries

```R
install.packages("pheatmap")
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "gplots", "RColorBrewer", "NMF", "Glimma"))
```

Load required R libaries

```R
library(pheatmap)
library(corncob)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(limma)
library(Glimma) # only works on pickles
```

Load data 

```R
setwd("/home/allie/domhain_RNAseq/05-global_diffabund")
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("/home/allie/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")  
# reorder columns by metadata 
metadata <- metadata[order(colnames(genecounts)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(genecounts)==metadata$sample_id) # should return all true
```

Convert gene counts to DGEList object for DE analysis

```R
dge <- DGEList(genecounts)
# check library sizes across samples
dge$samples
# set groups you want to compare
group <- as.factor(metadata$aliquot_type) # starting with aliquot type
# reset levels in group factor
levels(dge$samples$group) <- c("HCF", "HCE", "HCD", "ECE", "ECD", "DCD")
# add to DEG object
dge$samples$group <- group
# how many samples per group?
table(dge$samples$group)
```

Filter lowly expressed genes

```R
keep <- filterByExpr(dge, legacy = FALSE)
# how many genes removed due to low expression?
table(keep) 
# filter from dge object
dge <- dge[keep, keep.lib.sizes = FALSE]
# estimate dispersons and check that genes follow trend line
dge <- estimateDisp(dge)
# visualize dispersion estimates
png("BCV.png")
plotBCV(dge)
dev.off()
```

Get annotation file from HOMD

```bash
# wget https://www.homd.org/ftp/genomes/PROKKA/current/tsv/ALL_genomes.tsv
awk -F"\t" '{print $1 "\t" $4 "\t" $7}' ALL_genomes.tsv > annotations.txt
# remove spaces
sed -i 's/ /_/g' annotations.txt
# if no gene identifier, add NA
awk -F"\t" '{for(i=1; i<=3; i++) if ($i=="") $i="NA"; print $1 "\t" $2 "\t" $3}' annotations.txt > temp
mv temp annotations.txt
# add species to genes 
wget https://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
awk '{print $1}' annotations.txt | sed 's/_.*//' | while read line; do grep -w -m1 $line SEQID_info.txt | awk -F'\t' '{print $3 "_" $4}'; done > species


```

Adding annotations to our data

```R
# download associated data ()
homd.tsv <- read.table("~/domhain_RNAseq/03-star_map/homd_db/annotations.txt", header=T)
# filter by locus tag 
ann <- homd.tsv[homd.tsv$locus_tag %in% rownames(dge$counts),]
# check that locus tags match between the two dataframes
table(ann$locus_tag==rownames(dge$counts)) # should all return true
# add genes to dge object
dge$genes <- ann
```

Library normalization using TMM (Trimmed Mean of M-Values)

```R
dge <- calcNormFactors(dge, method="TMM")
# normalization factors will be displayed in table as 'norm.factors'
dge$samples
```

Create design matrix

```R
# create design matrix
design <- model.matrix(~ 0 + group)
# design
# clean up column names
colnames(design) <- levels(group)
design
```

Voom transform the data (adjust library sizes using norm.factors) -- this should be a smooth distribution around the line

```R
png("voom.png")
par(mfrow=c(1,1))
v <- voom(dge,design,plot = TRUE)
dev.off()
```

Test for differential expression

```R
# Fit the linear model
fit <- lmFit(v)
# names(fit)
```

Specify which groups you want to compare 

```R
cont.matrix <- makeContrasts(HCFvDCD=HCF - DCD, 
	levels=design)
# apply contrasts matrix to fit object to get stats and estimated parameters of our comparison
fit.cont <- contrasts.fit(fit, cont.matrix)
# perform Bayes shrinkage on the variances and get p values for our DE genes
fit.cont <- eBayes(fit.cont)
# summarize DE genes in contrast 
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
```

Test volcano plot

```R
png("volcano.png")
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$gene, main="HCFvDCD")
dev.off()
# save image so we can run glimma on pickles
save.image("global.RData")
```

Use glimma to create interactive volcano plot (on pickles) -- NOTE: this is really slow and maybe not so useful for the full dataset but keeping here for when we subset

```R
# scp allie@hillary.clemson.edu:/home/allie/domhain_RNAseq/05-global_diffabund/global.RData .
library(pheatmap)
library(corncob)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(limma)
library(Glimma)
load("global.RData")
# run interactive plot (will output html in ./volcano folder)
group2 <- group
levels(group2) <- c("HCF", "HCE", "HCD", "ECE", "ECD", "DCD")
# set browser URL on remote
# utils::browseURL("http://google.com", browser="/usr/bin/firefox")
glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", counts=v$E, groups=group2, status=summa.fit[,1],
         anno=fit.cont$product, side.main="gene", folder="volcano", launch=FALSE)
```

Get top DE genes

```R
topTable(fit.cont,coef=1,sort.by="p", number=50)
```










Ok, now that we have the process down, try running a comparison of H-CF vs D-CD
-------------------------------------------------------------------------------

```R
# subset metadata to only include wanted categories
sub.map <- metadata[metadata$aliquot_type == "H-CF" | metadata$aliquot_type == "D-CD",]
# subset genecounts file
sub.count <- genecounts[,colnames(genecounts) %in% row.names(sub.map)]
# reorder by rownames
sub.map <- sub.map[order(colnames(sub.count)),]
# check
table(colnames(sub.count)==sub.map$sample_id) # should return all true

# create dge object for edgeR
dge <- DGEList(sub.count, group=as.factor(sub.map$aliquot_type))
# reset levels in group factor
levels(dge$samples$group) <- c("H-CF", "D-CD")
# check library sizes across samples
dge$samples
# how many samples per group?
table(dge$samples$group)
keep <- filterByExpr(dge)
# how many genes removed due to low expression?
table(keep) 
# filter from dge object
dge <- dge[keep, ]

dge <- estimateDisp(dge, robust = TRUE)
# visualize dispersion estimates
png("BCV.H-CFvD-CD.png")
plotBCV(dge)
dev.off()
# create color scheme so you can see samples better in MDS plot
col.cell <- c("green", "red")[dge$samples$group]
# check
data.frame(dge$samples$group, col.cell)
# generate plot
png("MDS.H-CFvD-CD.png")
plotMDS(dge, col=col.cell)
dev.off()
```

DE calculation

```R
# fit model with robust estimation
fit <- glmQLFit(dge, design, robust = TRUE)
# differential expression 
qlf <- glmQLFTest(fit, coef=2)
# top differntially expressed genes
topTags(qlf)
```

















