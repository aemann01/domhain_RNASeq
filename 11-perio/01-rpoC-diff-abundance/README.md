# 1. DNA deseq2
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
# [1] "number of genes with adjusted p value lower than 0.05:  108"# adjusted p-value < 0.05
# out of 418 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 33, 7.9%
# LFC < 0 (down)     : 75, 18%
# outliers [1]       : 0, 0%
# low counts [2]     : 98, 23%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  105"
# out of 418 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 52, 12%
# LFC < 0 (down)     : 78, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 65, 16%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
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
filtered_res
```
# 3. Balance of taxa from DNA
```R
library(grid)
library(coda4microbiome)
library(tidyverse)
library(phyloseq)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
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
dat <- as.data.frame(t(as.data.frame(dat)))
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
y <- factor(dat$hiv_status) 
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

ps.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HEU")

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
dat <- as.data.frame(t(as.data.frame(dat)))
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
y <- factor(dat$hiv_status) 
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

pdf("./bal.HUUvHEU.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUUvHEU.pdf")
```
Corncob
```R
library(phyloseq)
library(corncob)
library(magrittr)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)
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
pdf("./diffab.hiv_status.dna.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.dna.pdf")
```
Corncob HEU vs HI
```R
library(phyloseq)
library(corncob)
library(magrittr)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)
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
pdf("./diffab.hiv_status.HEUvHI.dna.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.hiv_status.HEUvHI.dna.pdf")
```
Corncob HEU vs HI
```R
library(phyloseq)
library(corncob)
library(magrittr)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
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
# 4. rpoC bubble plots conditional mean (shout out to Allie)
Note: mean abundance is not a good way of looking at microbiome data because of a high prevalence of zeros -- therefore the mean abundance by itself can be misleading as it 1) dilutes true biological signals by averaging with potentially higher numbers of zeros in the data, 2) it doesn't distinguish between true absence versus technical zeros (undetected but present) and 3) it underrepresents low-abundance but consistently present taxa

Instead, we are using a log10 transformed conditional mean --> average abundance of microbial taxa calculated only from samples where the taxon is actually present (i.e. ignoring zeros)
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)

# change sp. HMT to oral taxon to match nomenclature differences
filter_sp <- c("Treponema_denticola", "Porphyromonas_gingivalis", "Tannerella_forsythia")
# do these actually exist in my phyloseq object?
species_df <- data.frame(QuerySpecies = filter_sp, stringsAsFactors = FALSE)
# Fuzzy match with taxonomy table
tax_levels <- paste0("V", 2:9)
# Find out if these species exist anywhere in the taxonomy and pull asvids
tax_df <- as.data.frame(tax_table(ps.dat)) %>%
  rownames_to_column("ASVID")

# Find matching ASVs for each query species
result_list <- lapply(species_df$QuerySpecies, function(x) {
  # Find rows where the species appears in any taxonomic level
  matches <- tax_df %>%
    filter(if_any(all_of(tax_levels), ~ . == x))
  
  if(nrow(matches) > 0) {
    data.frame(QuerySpecies = x, 
              ASVID = matches$ASVID,
              TaxLevel = apply(matches[, tax_levels], 1, function(row) {
                names(which(row == x))[1]
              }))
  } else {
    NULL
  }
})

# Combine results
matches_df <- bind_rows(result_list)
# Get all unique matching ASV IDs
matching_asvs <- unique(matches_df$ASVID)
# filter by matching ASVIDs
ps.dat.filt <- prune_taxa(matching_asvs, ps.dat)
ps.dat.filt
# first need to get a standardized species column (since the species show up in different levels)
matches_df <- matches_df %>%
  mutate(StandardizedSpecies = QuerySpecies)

# First calculate total number of samples in each comparison group
sample_counts <- sample_data(ps.dat.filt) %>% 
  as_tibble() %>% 
  count(hiv_status, tooth_health, name = "total_samples")
# Conditional mean abundance
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status", "tooth_health")) %>%
  # Calculate metrics
  group_by(Species, hiv_status, tooth_health) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
# create bubble plot
df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)

# I want to add in missing species that do not show up in our rpoC data to make the figures comparable
df <- df %>%
  # Ensure all species are included (even missing ones)
  complete(
    Species = filter_sp,
    hiv_status = unique(df$hiv_status),
    tooth_health = unique(df$tooth_health),
    fill = list(conditional_mean = 0, prevalence = 0)
  ) %>%
  # Reapply factor levels to Species
  mutate(Species = factor(Species, levels = filter_sp))

# order by species level
df <- df[order(df$Species),]
# set levels of x axis
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))
# df <- df %>%
  # filter(tooth_health != "E")
# and order of grid
# df$tooth_health <- factor(df$tooth_health, levels = c("H", "D"))
hivCols <- c("#40A0FA", "#FA78FA", "#8213A0")

pdf("HIVstat_rpoC_bubble_plot.pdf", width = 10)
ggplot(df,
  aes(
    x = hiv_status,
    y = Species,
    size = log10(conditional_mean + 1),
    color = hiv_status
  )
) +
  geom_point(alpha = ifelse(df$conditional_mean == 0, 0.1, 0.7)) +
  scale_size(range = c(0.5, 10), name = "Log10 Conditional Mean") +  
  scale_y_discrete(limits = rev) +  # Reverse y-axis order (top-to-bottom)
  scale_color_manual(values = hivCols) +
  labs(
    x = "Source",
    y = "Species",
    color = "Tooth hiv_status"
  ) +
  # facet_grid(. ~ tooth_health, switch = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Improve x-axis label readability
  )
dev.off()
system("~/.iterm2/imgcat HIVstat_rpoC_bubble_plot.pdf")

# make bar chart
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status")) %>%
  # Calculate metrics
  group_by(Species, hiv_status) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
# create bubble plot
df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)


gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
df$Species <- factor(df$Species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia"))
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("HIVstat_rpoC_bar_plot.pdf")
ggplot(df, aes(x = hiv_status, y = log10(conditional_mean + 1), fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = gencols) +
  theme_minimal() +
  # coord_flip() +  # This flips the x and y axes
  labs(x = "HIV Status", y = "log10 Count") +  # Custom axis labels
  theme(axis.text.x = element_text(size =10))  # Rotate x-axis labels
dev.off()
system("~/.iterm2/imgcat HIVstat_rpoC_bar_plot.pdf")

# do it by visit number
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status")) %>%
  # Calculate metrics
  group_by(Species, hiv_status, visit_num) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
# create bubble plot
df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)


gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
df$Species <- factor(df$Species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia"))
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("HIVstat_rpoC_bar_plot.visit_num.pdf")
ggplot(df, aes(x = visit_num, y = log10(conditional_mean + 1), fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = gencols) +
  theme_minimal() +
  # coord_flip() +  # This flips the x and y axes
  labs(x = "HIV Status", y = "log10 Count") +  # Custom axis labels
  theme(axis.text.x = element_text(size =10))  # Rotate x-axis labels
dev.off()
system("~/.iterm2/imgcat HIVstat_rpoC_bar_plot.visit_num.pdf")
```
# 5. non conditional mean
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)

setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)

# change sp. HMT to oral taxon to match nomenclature differences
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
rel <- microbiome::transform(glom, "log10")
ps.dat.filt <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
# do it by visit number
abundmean <- psmelt(ps.dat.filt) # %>%
#   # Join with standardized species names
#   left_join(
#     matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
#     by = "OTU"
#   ) %>%
#   # Filter for target species
#   filter(Species %in% filter_sp) %>% 
#   # Join with sample counts
#   left_join(sample_counts, by = c("hiv_status")) %>%
#   # Calculate metrics
#   group_by(Species, hiv_status, visit_num) %>%
#   summarize(
#     conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
#     prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
#     .groups = "drop"
#   ) %>%
#   # Replace NaN (from 0/0) with 0
#   mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
# create bubble plot
df <- abundmean
# reorder species column to match RNA seq figures
df$V8 <- factor(df$V8, levels = filter_sp)


gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
df$V8 <- factor(df$V8, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia"))
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("HIVstat_rpoC_bar_plot.log10.hiv_status.pdf")
ggplot(df, aes(x = hiv_status, y = Abundance, fill = V8)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = gencols) +
  theme_minimal() +
  # coord_flip() +  # This flips the x and y axes
  labs(x = "HIV Status", y = "log10 Count") +  # Custom axis labels
  theme(axis.text.x = element_text(size =10))  # Rotate x-axis labels
dev.off()
system("~/.iterm2/imgcat HIVstat_rpoC_bar_plot.log10.hiv_status.pdf")
```
# 6. DNA vs RNA comparison
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio//01-rpoC-diff-abundance")
#get relative abundance of dna
load("../../rpoc/ps.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, Abundance, V8)
red_dna$nucl <- "dna"
red_dna <- red_dna %>%
  rename(sample = Sample)
red_dna <- red_dna %>%
  rename(species = V8)
red_dna <- red_dna %>%
  rename(value = Abundance)
red_dna[red_dna$sample == "DM00008V1PQ16-2", ]
#get relative abundance of rna
rna_counts <- read.csv("../03-global-diff/species_reads.txt", sep='\t')
rna_counts$sample <- gsub(x = rna_counts$sample, pattern = "\\.red", replacement = "") 
rna_counts$sample <- gsub(x = rna_counts$sample, pattern = "\\.", replacement = "-") 
red_rna <- select(rna_counts, sample, Tannerella_forsythia, Porphyromonas_gingivalis, Treponema_denticola)
red_rna$nucl <- "rna"
red_rna <- melt(red_rna)
red_rna <- red_rna %>%
  rename(species = variable)
map$sample <- row.names(map)
meta <- as.data.frame(as.matrix(map)) %>% dplyr::select(sample, hiv_status)
red_rna <- left_join(meta, red_rna, by = "sample")
red_rna$value <- red_rna$value / 100
#combine
combined_df <- rbind(red_rna, red_dna)
# plot
combined_df$hiv_status <- factor(combined_df$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("red.RNAvDNA.bar.pdf")
ggplot() + geom_bar(data = combined_df, aes(x = sample, y = value, fill = nucl), position = "dodge", stat = "identity")+
	facet_grid(~ hiv_status, switch = "x", scales = "free_x")+
	theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.RNAvDNA.bar.pdf")
#corr
#plot
comb <- left_join(red_dna, red_rna, by = c("sample"="sample","species"="species", "hiv_status"="hiv_status"))
comb$hiv_status <- factor(comb$hiv_status, levels = c("HUU", "HEU", "HI"))
pdf("red.RNAvDNA.corr.pdf")
ggscatter(comb, x = "value.x", y = "value.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 facet_wrap(~species, scales = "free")+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.RNAvDNA.corr.pdf")
cor.test(comb$value.x, comb$value.y, method = "spearman") # not sig

```
# 7. Porportion of samples with DNA
```R 
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio//01-rpoC-diff-abundance")
#get relative abundance of dna
load("../../rpoc/ps.RData")
# proportion of different types of strep in pd vs pf
rel <- microbiome::transform(ps.dat, "compositional")
strep <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
# plot
# red complex by sample
data$Sample <- with(data, reorder(Sample, oral_hygiene_score))

pdf("samples.red.bar.pdf")
ggplot() + geom_bar(data = data, aes(x = Sample, y = Abundance, fill = V8), position = "stack", stat = "identity")+
	facet_grid(~ visit_num, switch = "x", scales = "free_x")+
	scale_color_manual(values=hivCols)+
	scale_fill_manual(values=hivCols)+
	theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./samples.red.bar.pdf")

# Calculate total number of samples per HIV status
total_samples_per_hiv <- filtered_data %>%
  group_by(hiv_status) %>%
  summarize(total_samples = n_distinct(Sample))

# Calculate count of samples per species and HIV status
species_counts_by_hiv <- filtered_data %>%
  group_by(V8, hiv_status) %>%
  summarize(sample_count = n_distinct(Sample), .groups = 'drop')

# Join total counts with species counts and calculate percentages
percentages <- species_counts_by_hiv %>%
  left_join(total_samples_per_hiv, by = "hiv_status") %>%
  mutate(percentage = (sample_count / total_samples) * 100)

# Print the results
print(percentages)
percentages$hiv_status <- factor(percentages$hiv_status, levels = c("HUU", "HEU", "HI"))
hivCols <- c("#40A0FA", "#FA78FA", "#8213A0")

pdf("samples.red.prop.pdf")
ggplot() + geom_bar(data = percentages, aes(x = hiv_status, y = percentage, fill = hiv_status), position = "dodge", stat = "identity")+
	facet_grid(~ V8, switch = "x", scales = "free_x")+
	scale_color_manual(values=hivCols)+
	scale_fill_manual(values=hivCols)+
	theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./samples.red.prop.pdf")
#find signficance
percentages$none <- percentages$total_samples - percentages$sample_count
sub_precent <- as.data.frame(select(percentages, V8, hiv_status, sample_count, none))
sub_precent <- filter(sub_precent, V8 == "Porphyromonas_gingivalis" & hiv_status != "HEU") %>% select(hiv_status, sample_count, none)
row.names(sub_precent) <- sub_precent$hiv_status
sub_precent <- select(sub_precent, sample_count, none)
chisq.test(sub_precent)


ps.rare <- rarefy_even_depth(ps.dat, rngseed=1, sample.size=0.99*min(sample_sums(ps.dat)), replace=F)

redcomplex <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(redcomplex, taxrank=rank_names(redcomplex)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
data$Abundance > 0
data_more <- data[data$Abundance > 20,]
hiv_group <- data_more %>%
  group_by(V8, hiv_status) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
hiv_precent <- hiv_group %>%
  group_by(V8, hiv_status) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)
pdf("sub.hiv.prop.pdf")
ggplot() + geom_bar(data = hiv_precent, aes(x = V8, y = percentage, fill = hiv_status), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.hiv.prop.pdf")

oral_group <- data_more %>%
  group_by(V8, oral_hygiene_score_remark) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
oral_precent <- oral_group %>%
  group_by(V8, oral_hygiene_score_remark) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)

oralCols <-c("#AA0A3B", "#F0F032", "#24B45A")
oral_precent$hiv_status <- factor(oral_precent$oral_hygiene_score_remark, levels=c("Good", "Fair", "Poor"))

pdf("sub.oral.prop.pdf")
ggplot() + geom_bar(data = oral_precent, aes(x = V8, y = percentage, fill = oral_hygiene_score_remark), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.oral.prop.pdf")

# find proportaion of samples
precent_samples <- data %>%
  group_by(V8) %>% 
  summarise(
    total_samples = n(),
    samples_above_20 = sum(Abundance >= 20, na.rm = TRUE),  
    percentage_above_20 = (samples_above_20 / total_samples) * 100  
  )

pdf("sub.red_all.prop.pdf")
ggplot() + geom_bar(data = precent_samples, aes(x = V8, y = percentage_above_20), position = "stack", stat = "identity")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.red_all.prop.pdf")

# ASV level
ps.rare <- rarefy_even_depth(ps.dat, rngseed=1, sample.size=0.99*min(sample_sums(ps.dat)), replace=F)

redcomplex <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
data <- psmelt(redcomplex) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
data$Abundance > 0
data$ASV <- paste0(data$V8,"_", data$OTU)
data_more <- data[data$Abundance > 10,]
hiv_group <- data_more %>%
  group_by(OTU, V8, hiv_status) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
hiv_precent <- hiv_group %>%
  group_by(OTU, V8, hiv_status) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)
pdf("sub.hiv.asv.prop.pdf")
ggplot() + geom_bar(data = hiv_precent, aes(x = OTU, y = percentage, fill = hiv_status), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  facet_wrap(~V8, scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.hiv.asv.prop.pdf")

oral_group <- data_more %>%
  group_by(OTU, V8, oral_hygiene_score_remark) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
oral_precent <- oral_group %>%
  group_by(OTU, V8, oral_hygiene_score_remark) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)

oralCols <-c("#8213A0", "#FA78FA", "#40A0FA")

pdf("sub.oral.asv.prop.pdf")
ggplot() + geom_bar(data = oral_precent, aes(x = OTU, y = percentage, fill = oral_hygiene_score_remark), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  facet_wrap(~V8, scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.oral.asv.prop.pdf")

# find proportaion of samples
precent_samples <- data %>%
  group_by(OTU, V8) %>% 
  summarise(
    total_samples = n(),
    samples_above_20 = sum(Abundance >= 20, na.rm = TRUE),  
    percentage_above_20 = (samples_above_20 / total_samples) * 100  
  )

pdf("sub.red_all.asv.prop.pdf")
ggplot() + geom_bar(data = precent_samples, aes(x = OTU, y = percentage_above_20), position = "stack", stat = "identity")+
  facet_wrap(~V8, scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./sub.red_all.asv.prop.pdf")
```
# 8. Porportion of all 1960 samples with DNA
```R 
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio//01-rpoC-diff-abundance")
#get relative abundance of dna
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))

rel <- microbiome::transform(ps.dat, "compositional")
strep <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
# plot
filtered_data <- data %>% filter(Abundance > 20)

# Calculate total number of samples per HIV status
total_samples_per_hiv <- filtered_data %>%
  group_by(hiv_status) %>%
  summarize(total_samples = n_distinct(Sample))

# Calculate count of samples per species and HIV status
species_counts_by_hiv <- filtered_data %>%
  group_by(V8, hiv_status) %>%
  summarize(sample_count = n_distinct(Sample), .groups = 'drop')

# Join total counts with species counts and calculate percentages
percentages <- species_counts_by_hiv %>%
  left_join(total_samples_per_hiv, by = "hiv_status") %>%
  mutate(percentage = (sample_count / total_samples) * 100)

# Print the results
print(percentages)
percentages$hiv_status <- factor(percentages$hiv_status, levels = c("HI", "HEU", "HUU"))
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

pdf("all_samples.red.prop.pdf")
ggplot() + geom_bar(data = percentages, aes(x = hiv_status, y = percentage, fill = hiv_status), position = "dodge", stat = "identity")+
  facet_grid(~ V8, switch = "x", scales = "free_x")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.red.prop.pdf")

strep <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$hiv_status <- factor(data$hiv_status, levels = c("HI", "HEU", "HUU"))
data$Abundance > 0
data_more <- data[data$Abundance > 0,]
# plot
pdf("all_samples.red.distro.pdf")
ggplot(data_more, aes(x = Abundance, fill = hiv_status)) +
  geom_density(position ="dodge", alpha = 0.4) +  # Adjust binwidth as needed
  facet_grid(~ V8, scales ="free") +
  scale_fill_manual(values = hivCols) + 
  theme_minimal() +
  labs(x = "Abundance", y = "Count", fill = "HIV Status")
dev.off()
system("~/.iterm2/imgcat ./all_samples.red.distro.pdf")

# Group by V8 and hiv_status, then summarize
# Compute summary statistics
summary_stats <- data %>%
  group_by(V8, hiv_status) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    Median = median(Abundance, na.rm = TRUE),
    # Mode = mode_function(Abundance),
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )
print(summary_stats)
summary_stats <- data_more %>%
  group_by(V8, hiv_status) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    Median = median(Abundance, na.rm = TRUE),
    # Mode = mode_function(Abundance),
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

#find signficance
percentages$none <- percentages$total_samples - percentages$sample_count
sub_precent <- as.data.frame(select(percentages, V8, hiv_status, sample_count, none))
sub_precent <- filter(sub_precent, V8 == "Porphyromonas_gingivalis" & hiv_status != "HEU") %>% select(hiv_status, sample_count, none)
row.names(sub_precent) <- sub_precent$hiv_status
sub_precent <- select(sub_precent, sample_count, none)
chisq.test(sub_precent)


#create barchart that shows species proportion for each hiv category
hiv_percentage_per_V8 <- data_more %>%
  group_by(V8, hiv_status) %>% 
  summarise(total_abundance_per_status = mean(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(V8) %>%
  mutate(percentage = (total_abundance_per_status / sum(total_abundance_per_status)) * 100) %>%
  ungroup()

pdf("all_samples.hiv.prop.pdf")
ggplot() + geom_bar(data = hiv_percentage_per_V8, aes(x = V8, y = percentage, fill = hiv_status), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.hiv.prop.pdf")
#create barchart that shows species proportion for each oral hygenince category
ps.rare <- rarefy_even_depth(ps.dat, rngseed=1, sample.size=0.99*min(sample_sums(ps.dat)), replace=F)

redcomplex <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(redcomplex, taxrank=rank_names(redcomplex)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
data$Abundance > 0
data_more <- data[data$Abundance > 20,]
hiv_group <- data_more %>%
  group_by(V8, hiv_status) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
hiv_precent <- hiv_group %>%
  group_by(V8, hiv_status) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)
pdf("all_samples.hiv.prop.pdf")
ggplot() + geom_bar(data = hiv_precent, aes(x = V8, y = percentage, fill = hiv_status), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.hiv.prop.pdf")

oral_group <- data_more %>%
  group_by(V8, Oral_Hygiene_Score_Remark) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
oral_precent <- oral_group %>%
  group_by(V8, Oral_Hygiene_Score_Remark) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)

oralCols <-c("#24B45A", "#F0F032", "#AA0A3B")
oral_precent$Oral_Hygiene_Score_Remark <- factor(oral_precent$Oral_Hygiene_Score_Remark, levels=c("Good", "Fair", "Poor"))
pdf("all_samples.oral.prop.pdf")
ggplot() + geom_bar(data = oral_precent, aes(x = V8, y = percentage, fill = Oral_Hygiene_Score_Remark), position = "stack", stat = "identity")+
  scale_color_manual(values=oralCols)+
  scale_fill_manual(values=oralCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.oral.prop.pdf")


# find proportaion of samples
precent_samples <- data %>%
  group_by(V8) %>% 
  summarise(
    total_samples = n(),
    samples_above_20 = sum(Abundance >= 20, na.rm = TRUE),  
    percentage_above_20 = (samples_above_20 / total_samples) * 100  
  )

pdf("all_samples.red_all.prop.pdf")
ggplot() + geom_bar(data = precent_samples, aes(x = V8, y = percentage_above_20), position = "stack", stat = "identity")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.red_all.prop.pdf")


# ASV level
ps.rare <- rarefy_even_depth(ps.dat, rngseed=1, sample.size=0.99*min(sample_sums(ps.dat)), replace=F)

redcomplex <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
data <- psmelt(redcomplex) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
data$Abundance > 0
data$ASV <- paste0(data$V8,"_", data$OTU)
data_more <- data[data$Abundance > 10,]
hiv_group <- data_more %>%
  group_by(OTU, V8, hiv_status) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
hiv_precent <- hiv_group %>%
  group_by(OTU, V8, hiv_status) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)
hiv_precent$hiv_status <- factor(hiv_precent$hiv_status, levels= c("HI", "HEU", "HUU"))
pdf("all_samples.hiv.asv.prop.pdf")
ggplot() + geom_bar(data = hiv_precent, aes(x = OTU, y = percentage, fill = hiv_status), position = "stack", stat = "identity")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  facet_wrap(~V8, scales = "free_x")+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.hiv.asv.prop.pdf")

oral_group <- data_more %>%
  group_by(OTU, V8, Oral_Hygiene_Score_Remark) %>% 
  summarise(mean_Abudnance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
oral_precent <- oral_group %>%
  group_by(OTU, V8, Oral_Hygiene_Score_Remark) %>%
  summarise(total_abundance = sum(mean_Abudnance)) %>%
  mutate(percentage = (total_abundance / sum(total_abundance)) * 100)

pdf("all_samples.oral.asv.prop.pdf")
ggplot() + geom_bar(data = oral_precent, aes(x = OTU, y = percentage, fill = Oral_Hygiene_Score_Remark), position = "stack", stat = "identity")+
  scale_color_manual(values=oralCols)+
  scale_fill_manual(values=oralCols)+
  facet_wrap(~V8, scales = "free_x")+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.oral.asv.prop.pdf")

# find proportaion of samples
precent_samples <- data %>%
  group_by(OTU, V8) %>% 
  summarise(
    total_samples = n(),
    samples_above_20 = sum(Abundance >= 20, na.rm = TRUE),  
    percentage_above_20 = (samples_above_20 / total_samples) * 100  
  )

pdf("all_samples.red_all.asv.prop.pdf")
ggplot() + geom_bar(data = precent_samples, aes(x = OTU, y = percentage_above_20), position = "stack", stat = "identity")+
  facet_wrap(~V8, scales = "free_x")+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.red_all.asv.prop.pdf")
```
# 9. Differences in CLR Abundance
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio/01-rpoC-diff-abundance")
#get relative abundance of dna
load("../../rpoc/ps.RData")
# tax_table(ps.dat) <- tax_table(ps.dat)[, c("V8")]
glom <- tax_glom(ps.dat, taxrank="V8")
rel <- microbiome::transform(glom, "clr")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU) 
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
pdf("porphyromonas_gingivalis.boxplot.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./porphyromonas_gingivalis.boxplot.pdf")

actino <- subset_taxa(rel, V8=="Treponema_denticola")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU)
pdf("treponema_denticola.boxplot.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./treponema_denticola.boxplot.pdf")

actino <- subset_taxa(rel, V8=="Tannerella_forsythia")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU)
pdf("tannerella_forsythia.boxplot.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./tannerella_forsythia.boxplot.pdf")
load("../../rpoc/ps.RData")

# do it for all samples
load("~/long_oral/master_phyloseq.RData")
glom <- tax_glom(ps.dat, taxrank="V8")
rel <- microbiome::transform(glom, "clr")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU) 
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
pdf("porphyromonas_gingivalis.boxplot.long.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./porphyromonas_gingivalis.boxplot.long.pdf")

actino <- subset_taxa(rel, V8=="Treponema_denticola")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU)
pdf("treponema_denticola.boxplot.long.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./treponema_denticola.boxplot.long.pdf")

actino <- subset_taxa(rel, V8=="Tannerella_forsythia")
data_bp <- psmelt(actino) #getting object into a dataframe
data_pb_2 <- select(data_bp, Sample, Abundance, hiv_status, OTU)
pdf("tannerella_forsythia.boxplot.long.pdf")
ggplot(data_pb_2, aes(x=factor(hiv_status, levels=hiv_stat),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+ #color dots by sample
  labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./tannerella_forsythia.boxplot.long.pdf")
```
# 10. ASVs across the 3 time period
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)
library(ggpubr)
setwd("/home/suzanne/rna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)

# change sp. HMT to oral taxon to match nomenclature differences
# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
rel <- microbiome::transform(ps.dat, "log10")
ps.dat.filt <- subset_taxa(rel, V8=="Porphyromonas_gingivalis")
# do it by visit number
abundmean <- psmelt(ps.dat.filt)

pdf("DNA.log10.visit_num.pdf")
ggplot(abundmean, aes(x = visit_num, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(values = gencols) +
  theme_minimal() +
  # coord_flip() +  # This flips the x and y axes
  labs(x = "HIV Status", y = "log10 Count") +  # Custom axis labels
  theme(axis.text.x = element_text(size =10))  # Rotate x-axis labels
dev.off()
system("~/.iterm2/imgcat DNA.log10.visit_num.pdf")

# glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
rel <- microbiome::transform(ps.dat, "clr")
ps.dat.filt <- subset_taxa(rel, V8=="Porphyromonas_gingivalis")
# do it by visit number
abundmean <- psmelt(ps.dat.filt)

pdf("pging_ASV.boxplot.pdf")
ggplot(abundmean, aes(x=as.factor(visit_num),y=Abundance))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  facet_wrap(~OTU) +
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./pging_ASV.boxplot.pdf")
```