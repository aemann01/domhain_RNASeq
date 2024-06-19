# HIV differential abundance 

Activate environment (running on hillary)

```bash
conda activate 2024-HIV_RNASeq
```

Load required R libaries

```R
library(pheatmap)
library(corncob)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(limma)
```

Load data 

```R
setwd("/home/allie/domhain_RNAseq/06-hiv_diffabund")
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
group <- as.factor(metadata$hiv_status) 
# reset levels in group factor
levels(dge$samples$group) <- c("HI", "HEU", "HUU")
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
awk '{print $1}' annotations.txt | sed 's/_.*//' > species
sed -i 's/locus/speciesID/' species
# merge
paste annotations.txt species > temp
mv temp annotations.txt
```

Get species names from seqid file

```bash
mamba install ipython
mamba install pandas
```

```python
# ipython
import pandas as pd 
ann = pd.read_csv("annotations.txt", delimiter="\t")
seqid = pd.read_csv("SEQID_info.txt", delimiter="\t")
# what column names?
for col in seqid.columns:
    print(col)
# merge together on species id column
merge = pd.merge(ann, seqid, left_on="speciesID", right_on="SEQ_ID", how="inner", indicator=True)
# sanity check -- did we have matches in both dataframes?
merge[merge._merge == "left_only"] # should return empty dataframe
merge[merge._merge == "right_only"] # should return empty dataframe
# clean up to only keep wanted columns
merge = merge.drop(["Sequence_Source", "Combined_Size", "Contigs", "speciesID", "Strain", "HMT_ID", "product", "_merge"], axis=1)
# strip all whitespace
merge = merge.apply(lambda x: x.str.strip() if x.dtype=="object" else x)
# write to file
merge.to_csv("annotations.merge.txt", sep="\t", index=False, na_rep="none")
```

Another check to make sure there are no spaces to screw up R

```bash
sed -i 's/ /_/g' annotations.merge.txt
```

Adding annotations to our data

```R
homd <- read.table("~/domhain_RNAseq/03-star_map/homd_db/annotations.merge.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$locus_tag %in% rownames(dge$counts),]
# check that locus tags match between the two dataframes
table(rownames(dge$counts)==ann$locus_tag) # should all return true
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
cont.matrix <- makeContrasts(HIvHUU=HI - HUU, HIvHEU=HI - HEU, HEUvHUU=HEU - HUU,
	levels=design)
# apply contrasts matrix to fit object to get stats and estimated parameters of our comparison
fit.cont <- contrasts.fit(fit, cont.matrix)
# perform Bayes shrinkage on the variances and get p values for our DE genes
fit.cont <- eBayes(fit.cont)
# summarize DE genes in contrast 
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
#        HIvHUU HIvHEU HEUvHUU
# Down     3866      0       0
# NotSig 672615 677062  677059
# Up        581      0       3
```

Test volcano plot

```R
png("volcano-HIvHUU.png")
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$gene, main="HI v HUU")
dev.off()
# save image so we can run glimma on pickles
save.image("hivDE.RData")
```

Get all DE genes (over estimating so that we can filter the table later)

```R
allDE <- topTable(fit.cont,coef=1, adjust.method="BH", sort.by="P", number=5000)
write.csv(allDE, "limma_DEresults.txt")
```

Test heatmap

```R
coolmap(fit.cont)
```








```R
install.packages("ggpubr")
library(ggpubr)

options(ggrepel.max.overlaps = Inf)
p <- ggmaplot(fit.cont, main = expression("HI" %->% "HUU"),
   fdr = 0.05, fc = 4, size = 0.4,
   palette = c("#1465AC", "#B31B21", "darkgray"),
   genenames = as.vector(row.names(res)),
   legend = "top", top = 0,
   font.label = c("bold", 11),
   font.legend = "bold",
   font.main = "bold",
   label.rectangle = T,
   label.select = c("SEQF2748_01243", "SEQF2748_01244", "SEQF2748_01245", "SEQF2016_01287", "SEQF2016_01288", "SEQF2016_01289", "SEQF1998_00217", "SEQF1998_00218", "SEQF1998_00219", "SEQF1964_00410", "SEQF1964_00411", "SEQF1964_00412", "SEQF1964_00414", "SEQF2625_00202", "SEQF2625_00203", "SEQF2625_00204", "SEQF1766_00037", "SEQF1766_00038", "SEQF1766_00039"),
   ggtheme = ggplot2::theme_minimal())
p




















