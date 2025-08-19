# 1. Random forest for predicting oral hygiene score for RNAseq data
```R
# install.packages("kernelshap")
# install.packages("ranger")
# install.packages("vip")
# devtools::install_github("ModelOriented/shapviz")
# install.packages("reshape2")
library(ranger) #random forest package
library(kernelshap) #shapley
library(vip)
library(shapviz)
library(reshape2)
library(ggplot2)
setwd("/home/suzanne/rna_dohmain/11-perio/13-random-forest")
set.seed(546446)

metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
# add pseudocount to avoid errors with size factor estimation
# subcount <- subcount + 1
# reorder columns by metadata 
submap <- submap[order(colnames(subcount)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(subcount)==submap$sample_id) # should return all true
# normalize read counts
count_tab_norm <- sweep(subcount, 2, colSums(subcount), '/')*100
#rename to species, gene, and SEQF
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")

# run random forest
fit <- ranger(var ~ ., data = asv_tab_var,
	scale.permutation.importance = TRUE, importance = 'permutation')
```
# 2. Using long dataset
```R
library(ranger) #random forest package
library(kernelshap) #shapley
library(vip)
library(shapviz)
library(reshape2)
library(ggplot2)
library(phyloseq, verbose=F)
library(tidyverse, verbose=F)
library(caTools, verbose=F)
library(ranger, verbose=F)
library(fastshap, verbose=F)
library(kernelshap, verbose=F)
library(iml, verbose=F)
library(vip, verbose=F)
library(shapviz, verbose=F)
library(reshape2, verbose=F)
library(janitor)
setwd("/home/suzanne/rna_dohmain/11-perio/13-random-forest")
set.seed(546446)
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
# glom and rename to species
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
temp
# save copy to reduce time on previous command
glom <- temp
rel <- microbiome::transform(glom, "compositional")
# pull data
dat <- t(as.data.frame(otu_table(rel)))
map <- as.data.frame(as.matrix(sample_data(rel))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(rel), "matrix")
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
# add the oral hygienge score
meta_core <- map[, c("rowname", "Oral_Hygiene_Score")]
rownames(meta_core) <- meta_core$rowname
meta_core$rowname <- NULL  # Not needed anymore
table(rownames(dat)==rownames(meta_core)) # should all return true
asv_tab_var <- cbind(dat, meta_core)
asv_tab_var$Oral_Hygiene_Score <- as.numeric(asv_tab_var$Oral_Hygiene_Score)
#make random forrest model
colnames(asv_tab_var) <- gsub("\\.", "", colnames(asv_tab_var))
asv_tab_var <- janitor::clean_names(asv_tab_var)
fit <- ranger(oral_hygiene_score ~ ., data = asv_tab_var,
	scale.permutation.importance = TRUE, importance = 'permutation')

#variable importance
pdf("./rf.OHS.importance.pdf")
vip(fit, title = "Variable Importance")
dev.off()
system("~/.iterm2/imgcat ./rf.OHS.importance.pdf")
ranger::importance(fit)

#goodness of fit
p.ra <- predict(fit, data=asv_tab_var)
str(p.ra)
par(mfrow=c(1,2))
pdf("goodfit.OHS.pdf")
plot(asv_tab_var$Oral_Hygiene_Score ~ p.ra$predictions, asp=1, pch=20, xlab="fitted", ylab="actual", xlim=c(0,3.3),
          ylim=c(0,3.7),main="OHS, Random Forest")
grid(); abline(0,1)
dev.off()
system("~/.iterm2/imgcat ./goodfit.OHS.pdf")

#out of bag prediction
summary(fit$predictions)
summary(fit$predictions - p.rf.oob)  # difference
```
# 3. Using long dataset with the remark
```R
library(ranger) #random forest package
library(kernelshap) #shapley
library(vip)
library(shapviz)
library(reshape2)
library(ggplot2)
library(phyloseq, verbose=F)
library(tidyverse, verbose=F)
library(caTools, verbose=F)
library(ranger, verbose=F)
library(fastshap, verbose=F)
library(kernelshap, verbose=F)
library(iml, verbose=F)
library(vip, verbose=F)
library(shapviz, verbose=F)
library(reshape2, verbose=F)
setwd("/home/suzanne/rna_dohmain/11-perio/13-random-forest")
set.seed(546446)
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
# glom and rename to species
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
temp
# save copy to reduce time on previous command
glom <- temp
rel <- microbiome::transform(glom, "compositional")
# pull data
dat <- t(as.data.frame(otu_table(rel)))
map <- as.data.frame(as.matrix(sample_data(rel))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(rel), "matrix")
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
# add the oral hygienge score
meta_core <- map[, c("rowname", "Oral_Hygiene_Score_Remark")]
rownames(meta_core) <- meta_core$rowname
meta_core$rowname <- NULL  # Not needed anymore
table(rownames(dat)==rownames(meta_core)) # should all return true
asv_tab_var <- cbind(dat, meta_core)
asv_tab_var$Oral_Hygiene_Score_Remark <- as.factor(asv_tab_var$Oral_Hygiene_Score_Remark)
#make random forrest model
colnames(asv_tab_var) <- gsub("\\.", "", colnames(asv_tab_var))
asv_tab_var <- janitor::clean_names(asv_tab_var)
fit <- ranger(oral_hygiene_score_remark ~ ., data = asv_tab_var,
	scale.permutation.importance = TRUE, importance = 'permutation')
fit
fit$confusion.matrix
