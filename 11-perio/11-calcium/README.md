# 1. Run DESeq low vs normal calcium
```sh
awk '{print $1}' map.txt | while read -r line; do
    result=$(grep -w "$line" map_domhain_long_2.txt)
    if [ -n "$result" ]; then
        echo "$result" | awk -F "\t" '{print $69, $70, $71, $72}'
    else
        echo "BLEAK"
    fi
done
```
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
metadata$Ca_category <- ifelse(metadata$total_Ca_mg <= 400, "low", ifelse(metadata$total_Ca_mg >= 600, "high", "normal"))
metadata <- na.omit(metadata)

# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$Ca_category == "normal" | metadata$Ca_category == "low",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~Ca_category)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$Ca_category <- factor(star_results$Ca_category, levels=c("low", "normal"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  5386"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 6, 0.012%
# LFC < 0 (down)     : 5380, 11%
# outliers [1]       : 1111, 2.2%
# low counts [2]     : 40159, 79%
# (mean count < 3)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="Ca_category_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  5013"
# out of 51006 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 7, 0.014%
# LFC < 0 (down)     : 6974, 14%
# outliers [1]       : 1111, 2.2%
# low counts [2]     : 38270, 75%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_Ca-NormalvLow.txt", quote=F, sep="\t")
save.image("deseq_results_Ca-NormalvLow.RData")
```
Valcona Plot
```R
load("deseq_results_Ca-NormalvLow.RData")
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
pdf("volcano-LowvNormal.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-LowvNormal.red.pdf")
```
# 2. Using global abundance
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
metadata$Ca_category <- ifelse(metadata$total_Ca_mg <= 400, "low", ifelse(metadata$total_Ca_mg >= 600, "high", "normal"))
metadata <- na.omit(metadata)
table(metadata$Ca_category)
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$Ca_category == "normal" | metadata$Ca_category == "low",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~Ca_category)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$Ca_category <- factor(star_results$Ca_category, levels=c("low", "normal"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  111418"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 118, 0.0018%
# LFC < 0 (down)     : 111300, 1.7%
# outliers [1]       : 68939, 1%
# low counts [2]     : 6210710, 94%
# (mean count < 11)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="Ca_category_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  210857"# out of 51006 with nonzero total read count
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 123978, 1.9%
# LFC < 0 (down)     : 147149, 2.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 5564940, 84%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_global-LowvNormal.txt", quote=F, sep="\t")
save.image("deseq_results_global-LowvNormal.RData")
```
Valcona Plot
```R
load("deseq_results_global-LowvNormal.RData")
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
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$species, res_ord$gene, sep=" "), ""),
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
pdf("volcano-LowvHigh.global.pdf", width=30, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-LowvHigh.global.pdf")
```
# 3. Look at beta-diversity from full analysis
```R
library(tidyr, verbose=F)
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(coda4microbiome, verbose=F)
library(microbiome, verbose=F)
library(ecole, verbose=F)
library(ggpubr)

#load data
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
load("~/long_oral/master_phyloseq.RData")
hiv_stat <- c("HI", "HEU", "HUU")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
meta$Ca_category <- ifelse(meta$total_Ca_mg <= 440, "low", ifelse(meta$total_Ca_mg >= 600, "high", "normal"))
# compare for unique individuals
meta$sample_visit <- paste0(meta$study_id,"V",meta$visit_num)
meta_sub <- meta[!is.na(meta$total_Ca_mg), ]
meta_unique <- meta_sub[!duplicated(meta_sub$sample_visit), ]
table(meta_unique$Ca_category) %>% (table(meta_unique$Ca_category) / sum(table(meta_unique$Ca_category))) * 100
cor.test(meta_unique$total_Ca_mg, meta_unique$Oral_Hygiene_Score, method = "spearman")

pdf("CavHIV.long.pdf")
ggplot(meta_unique, aes(x=factor(hiv_status, levels=hiv_stat),y=total_Ca_mg))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./CavHIV.long.pdf")

ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
ps.sub <- subset_samples(ps.dat, Ca_category=="low" | Ca_category=="normal" | Ca_category=="high")

hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

ps.dat.clr <- microbiome::transform(ps.sub, transform="clr", target="OTU")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~Ca_category)
pdf("./bdiv_cap.long_Ca.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="Ca_category") + 
    theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.long_Ca.pdf")

metadata <- as(sample_data(ps.dat.clr), "data.frame")
clr.dist <- dist(otu_table(ps.dat.clr), method="euclidean")
adonis2(clr.dist ~ hiv_status * aliquot_type * sex * Ca_category, data=metadata)
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$Ca_category, method="euclidean", padj= "fdr") # signficant 

#            pairs SumOfSqs  F.Model          R2  pval  p.adj
# 1  low vs normal 9544.660 1.409394 0.001212477 0.001 0.0015
# 2    low vs high 9330.741 1.376650 0.001276595 0.001 0.0015
# 3 normal vs high 8812.104 1.404631 0.011382322 0.002 0.0020
```
# 4. coda4microbiome long data set Ca
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
meta$Ca_category <- ifelse(meta$total_Ca_mg <= 440, "low", ifelse(meta$total_Ca_mg >= 600, "high", "normal"))
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
ps.sub <- subset_samples(ps.dat, Ca_category=="low" | Ca_category=="normal" | Ca_category=="high")

# overall t
temp <- tax_glom(ps.sub, taxrank=rank_names(ps.sub)[7])
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
temp <- subset_samples(temp, !is.na(sample_data(temp)$total_Ca_mg))
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.overall_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.overall_Ca.long.pdf")

# HI oral hygiene score
temp <- tax_glom(ps.sub, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HI")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HI_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HI_Ca.long.pdf")

# HEU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HEU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HEU_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HEU_Ca.long.pdf")

# HUU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HUU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HUU_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUU_Ca.long.pdf")
```
# 5. Differences between V1 and V2 for Calcium
```R
library(tidyr, verbose=F)
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(coda4microbiome, verbose=F)
library(microbiome, verbose=F)
library(ecole, verbose=F)
library(ggpubr)
library(ggdist)
set.seed(12345)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
load("~/long_oral/master_phyloseq.RData")
hiv_stat <- c("HI", "HEU", "HUU")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
hiv_stat <- c("HI", "HEU", "HUU")
meta$hiv_status <- factor(meta$hiv_status, levels = hiv_stat)
meta$Ca_category <- ifelse(meta$total_Ca_mg <= 440, "low", ifelse(meta$total_Ca_mg >= 600, "high", "normal"))
# compare for unique individuals
meta$sample_visit <- paste0(meta$study_id,"V",meta$visit_num)
meta_sub <- meta[!is.na(meta$total_Ca_mg), ]
meta_unique <- meta_sub[!duplicated(meta_sub$sample_visit), ]
meta_unique$visit_num <- as.factor(meta_unique$visit_num)
meta_unique <- meta_unique[meta_unique$total_Ca_mg < 2000, ]

pdf("total_Ca_mg.visit.pdf")
ggplot(meta_unique, aes(x = visit_num, y = as.numeric(total_Ca_mg))) + 
    geom_point(aes(color = Ca_category), position =  position_jitter(height = 0.05)) +
    scale_color_manual(values = c("high" = "#006164", "normal" = "#EDA247", "low" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns = TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    ylab("Simplified Oral Hygiene Score") +
    xlab("HIV status") +
    theme(legend.title = element_blank()) +
    facet_wrap(~hiv_status)
dev.off()
system("~/.iterm2/imgcat ./total_Ca_mg.visit.pdf")

pdf("Oral_Hygiene_Score.visit.pdf")
ggplot(meta_unique, aes(x=visit_num,y=Oral_Hygiene_Score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~hiv_status)
dev.off()
system("~/.iterm2/imgcat ./Oral_Hygiene_Score.visit.pdf")

diff.v1v2 <- meta_unique %>%
  select(study_id, hiv_status, visit_num, total_Ca_mg) %>%
  filter(!is.na(total_Ca_mg)) %>%
  pivot_wider(names_from = visit_num, values_from = total_Ca_mg, names_prefix = "visit_") %>%
  mutate(delta_Ca = visit_2 - visit_1)

diff.v1v2 <- diff.v1v2[!is.na(diff.v1v2$delta_Ca), ]
diff.v1v2$hiv_status <- factor(diff.v1v2$hiv_status, levels = c("HI", "HEU", "HUU"))
# what is the standard deviation
diff.v1v2 %>%
  dplyr::group_by(hiv_status) %>%
  dplyr::summarise(
    count = dplyr::n(),
    median = median(delta_Ca, na.rm = TRUE),
    SD = sd(delta_Ca, na.rm = TRUE)
  )


pdf("./calcium.bplot.v1v2.pdf")
options(repr.plot.width = 7, repr.plot.height =5)
ggplot(diff.v1v2, aes(x=hiv_status, y=delta_Ca, fill=hiv_status)) + 
    geom_pwc(label = "{p.adj.format}", hide.ns =FALSE, method = "wilcox_test", p.adjust.method = "none") +
    theme_bw() +
    stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA, normalize="all", scale=0.5) + 
    geom_boxplot(width = 0.12, alpha = 0.5, outlier.color=NA) +
    stat_dots(side = "left", justification = 1.1, binwidth = 0.25) +
    scale_fill_manual(values=c("HI"="#8213A0", "HEU"="#FA78FA", "HUU"="#40A0FA")) +
    scale_colour_manual(values=c("HI"="#8213A0", "HEU"="#FA78FA", "HUU"="#40A0FA")) +
    scale_y_continuous(limits=c(10,45)) +
    coord_flip()
dev.off()
system("~/.iterm2/imgcat ./calcium.bplot.v1v2.pdf")
```
# 2. Using global abundance and treating as continuous variable
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-calcium")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
metadata$Ca_category <- ifelse(metadata$total_Ca_mg <= 440, "low", ifelse(metadata$total_Ca_mg >= 600, "high", "normal"))
metadata <- na.omit(metadata)
table(metadata$Ca_category)
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$Ca_category == "normal" | metadata$Ca_category == "low",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~total_Ca_mg)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$Ca_category <- factor(star_results$Ca_category, levels=c("low", "normal"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  111418"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 118, 0.0018%
# LFC < 0 (down)     : 111300, 1.7%
# outliers [1]       : 68939, 1%
# low counts [2]     : 6210710, 94%
# (mean count < 11)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="Ca_category_normal_vs_low", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  210857"# out of 51006 with nonzero total read count
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 123978, 1.9%
# LFC < 0 (down)     : 147149, 2.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 5564940, 84%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_global-calcium.txt", quote=F, sep="\t")
save.image("deseq_results_global-calcium.RData")
```
Valcona Plot
```R
load("deseq_results_global-LowvNormal.RData")
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
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
res_ord$GeneInfo <- paste(res_ord$species,res_ord$gene)
res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$tag %in% labgenes, paste(res_ord$species, res_ord$gene, sep=" "), ""),
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
pdf("volcano-HIvHUU.global.pdf", width=30, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.global.pdf")
```
