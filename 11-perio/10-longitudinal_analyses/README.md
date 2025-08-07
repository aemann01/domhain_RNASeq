# 1. How does red complex RNA change over V1 to V2
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
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

# list samples of interest
twice <- c("DM00101V2PQ16","DM00101V3PQ16","DM00305V2PQ16","DM00305V3PQ16","DM00338V2PQ16","DM00338V3PQ16","DM00008V1PQ16-2","DM00008V2PQ16","DM00023V1PQ16-1","DM00023V2PQ16","DM00035V1PQ16","DM00035V2PQ16","DM00074V1PQ55","DM00074V2PQ55","DM00142V1PQ16","DM00142V2PQ16","DM00255V1PQ26","DM00255V2PQ26","DM00266V1PQ16","DM00266V2PQ16","DM00464V1PQ65","DM00464V2PQ65","DM00563V1PQ55","DM00563V2PQ55")

homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")

# create log2 expression of red complex genes 
# filter metadata
map.twice <- submap[rownames(submap) %in% twice,]
count.twice <- subcount[,colnames(subcount) %in% twice]
# remove any genes with a sum count of zero
count.twice <- count.twice[rowSums(count.twice) != 0,]
# filter annotations by those found in our genecounts file
ann <- homd[homd$tag %in% rownames(count.twice),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(count.twice), rownames(ann)))]
count.twice <- count.twice[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.twice)==rownames(ann)) # should all return true
# if all are true, merge together
count.twice <- cbind(count.twice, ann)
# get list of genera to pull from rpoC data later
count.twice.gen <- unique(count.twice$species)

# collapse by species and sum across rows
merge.count <- count.twice %>% group_by(species) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# get total count across all rows at the species level 
group_sums <- group_sums %>% mutate(total=rowSums(across(where(is.numeric))))
# get total counts to use as denominator 
tot <- sum(group_sums$total)
thresh <- 0.01 * tot # less than 1% of total

# collapse low abundant groups
collapsed <- group_sums %>%
  mutate(species = ifelse(total < thresh, "Other", species)) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
# remove total column
collapsed <- subset(collapsed, select = -total)

# Create stacked histograms for each sample colored by species
# set color palette
gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("RNA_log10_Red_200days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat RNA_log10_Red_200days.pdf")
```
# 2. How does red complex RNA change over V1 to V3
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
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

# list samples of interest
triple <- c("DM00044V1PQ16","DM00044V2PQ16","DM00044V3PQ16","DM00103V1PQ16","DM00103V2PQ16","DM00103V3PQ16","DM00254V1PQ55","DM00254V2PQ55","DM00254V3PQ55","DM00519V1PQ16","DM00519V2PQ16","DM00519V3PQ16")

homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")

# create log2 expression of red complex genes 
# filter metadata
map.triple <- submap[rownames(submap) %in% triple,]
count.triple <- subcount[,colnames(subcount) %in% triple]
# remove any genes with a sum count of zero
count.triple <- count.triple[rowSums(count.triple) != 0,]
# filter annotations by those found in our genecounts file
ann <- homd[homd$tag %in% rownames(count.triple),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(count.triple), rownames(ann)))]
count.triple <- count.triple[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.triple)==rownames(ann)) # should all return true
# if all are true, merge together
count.triple <- cbind(count.triple, ann)
# get list of genera to pull from rpoC data later
count.triple.gen <- unique(count.triple$species)

# collapse by species and sum across rows
merge.count <- count.triple %>% group_by(species) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# get total count across all rows at the species level 
group_sums <- group_sums %>% mutate(total=rowSums(across(where(is.numeric))))
# get total counts to use as denominator 
tot <- sum(group_sums$total)
thresh <- 0.01 * tot # less than 1% of total

# collapse low abundant groups
collapsed <- group_sums %>%
  mutate(species = ifelse(total < thresh, "Other", species)) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
# remove total column
collapsed <- subset(collapsed, select = -total)

# Create stacked histograms for each sample colored by species
# set color palette
gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("RNA_log10_Red_400days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat RNA_log10_Red_400days.pdf")
```
# 3. Load in rpoC ASV count data for V1 to V2
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("../../rpoc/ps.RData")
ps.dat.filt  <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
twice <- c("DM00101V2PQ16","DM00101V3PQ16","DM00305V2PQ16","DM00305V3PQ16","DM00338V2PQ16","DM00338V3PQ16","DM00008V1PQ16-2","DM00008V2PQ16","DM00023V1PQ16-1","DM00023V2PQ16","DM00035V1PQ16","DM00035V2PQ16","DM00074V1PQ55","DM00074V2PQ55","DM00142V1PQ16","DM00142V2PQ16","DM00255V1PQ26","DM00255V2PQ26","DM00266V1PQ16","DM00266V2PQ16","DM00464V1PQ65","DM00464V2PQ65","DM00563V1PQ55","DM00563V2PQ55")

# filter samples to only include our twice samples
ps.dat.twice <- prune_samples(twice, ps.dat.filt)
ps.dat.twice
# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.twice, taxrank = "V8")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))
# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V8)
asvtab <- as.data.frame(asvtab)
# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by genus (V25) and sum across rows
merge.asvtab <- asvtab %>% group_by(V25) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# rename columns
merge.asvtab <- merge.asvtab %>% rename(species = V25)

# groups to collapse by
keepgroup <- c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia")
# Collapse and sum rows
collapsed <- merge.asvtab %>%
  mutate(species = ifelse(species %in% keepgroup, species, "Other")) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))

# set color palette
gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)

# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("rpoC_log10_Red_200days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat rpoC_log10_Red_200days.pdf")
```
# 4. Load in rpoC ASV count data for V1 to V3
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("../../rpoc/ps.RData")
ps.dat.filt  <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
triple <- c("DM00044V1PQ16","DM00044V2PQ16","DM00044V3PQ16","DM00103V1PQ16","DM00103V2PQ16","DM00103V3PQ16","DM00254V1PQ55","DM00254V2PQ55","DM00254V3PQ55","DM00519V1PQ16","DM00519V2PQ16","DM00519V3PQ16")

# filter samples to only include our triple samples
ps.dat.triple <- prune_samples(triple, ps.dat.filt)
ps.dat.triple
# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.triple, taxrank = "V8")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))
# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V8)
asvtab <- as.data.frame(asvtab)
# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by genus (V13) and sum across rows
merge.asvtab <- asvtab %>% group_by(V13) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# rename columns
merge.asvtab <- merge.asvtab %>% rename(species = V13)

# groups to collapse by
keepgroup <- c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia")
# Collapse and sum rows
collapsed <- merge.asvtab %>%
  mutate(species = ifelse(species %in% keepgroup, species, "Other")) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))

# set color palette
gencols <- c(Porphyromonas_gingivalis = "#340043", 
			Treponema_denticola = "#FBE51F", 
			Tannerella_forsythia = "#1E7F7A")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)

# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("rpoC_log10_Red_400days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat rpoC_log10_Red_400days.pdf")
```
# 5. Compare oral hygiene score by visit
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))

# V1 oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HI")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# visit one first
glom <- subset_samples(glom, visit_num == "1")
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
y <- as.numeric(datmerge$Oral_Hygiene_Score)
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

pdf("./bal.V1_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.V1_oral.long.pdf")

#V2
# save copy to reduce time on previous command
glom <- temp
# visit one first
glom <- subset_samples(glom, visit_num == "2")
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
y <- as.numeric(datmerge$Oral_Hygiene_Score)
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

pdf("./bal.V2_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.V2_oral.long.pdf")

#V2
# save copy to reduce time on previous command
glom <- temp
# visit one first
glom <- subset_samples(glom, visit_num == "3")
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
y <- as.numeric(datmerge$Oral_Hygiene_Score)
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

pdf("./bal.V3_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.V3_oral.long.pdf")


# for HI
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
glom <- subset_samples(glom, hiv_status == "HI")
glom <- filter_taxa(glom, function(x) sum(x > 500) > (0.01*length(x)), TRUE)
# pull data
dat <- as.data.frame(otu_table(glom))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# first need to get comparable variable to merge our data on
# head(dist_diff)
map$Var1 <- paste(map$study_id, map$FDI_code, sep=".")
# head(map)
# head(dist_diff)
# merge dist_diff with map
# map <- merge(as.data.frame(map), dist_diff, by="Var1")
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
# x <- select_if(x, is.numeric)
# dim(x)
x <- as.matrix(x) # ASV abundance table
x_time <- as.numeric(datmerge$visit_num) # time point
subject_id <- datmerge$study_id # subject id
y <- datmerge$sex # response value (here is oral hygience score)
y <- as.numeric(as.character(y))
ini_time <- 1
end_time <- 3
set.seed(852)

#balance of taxa
bal <- coda_glmnet_longitudinal(x, y, x_time, subject_id, ini_time, end_time, nfolds=10, showPlots = FALSE)
```
# 6. How does oral hygiene score change over V1 to V3
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
map <- sample_data(ps.dat)
map$studyID_FDI <- paste(map$study_id, ".", map$FDI_code, sep="")
map.v1v3 <- map[map$visit_num == 1 | map$visit_num == 3,]
# count the number of times that a study ID + FDI code appears in the data
counts <- data.frame(table(map.v1v3$studyID_FDI))
map.v1v3 <- map.v1v3[map.v1v3$studyID_FDI %in% counts$Var1[counts$Freq>1],]
# reorder
map.v1v3 <- map.v1v3[order(map.v1v3$studyID_FDI),]
# number of samples
dim(map.v1v3)
oral_score <- as.data.frame(matrix(map.v1v3$Oral_Hygiene_Score, ncol=2, byrow=TRUE))
mean(oral_score$V2 - oral_score$V1)
# average difference between date of visit one and visit two
# Step 1: Filter for V1 and V3 visits
map.v1 <- as.data.frame(as.matrix((map.v1v3 [map.v1v3 $visit_num == 1,])))
map.v1 <- map.v1 %>% rename_with(~ ifelse(. == "study_id", ., paste0(. , "_v1")))
map.v3 <- as.data.frame(as.matrix(map.v1v3 [map.v1v3 $visit_num == 3,]))
map.v3 <- map.v3 %>% rename_with(~ ifelse(. == "study_id", ., paste0(. , "_v3")))

# Step 2: Merge the V1 and V3 dataframes on 'study_id' to compare the same subjects
oral_v1v3 <- left_join(map.v1, map.v3, by = join_by(study_id ==  study_id))
oral_v1v3 <- oral_v1v3[!(oral_v1v3$Oral_Hygiene_Score_v1 %in% NA),]
oral_v1v3 <- oral_v1v3[!(oral_v1v3$Oral_Hygiene_Score_v3 %in% NA),]
oral_v1v3$diff_OHS <- as.numeric(oral_v1v3$Oral_Hygiene_Score_v3) - as.numeric(oral_v1v3$Oral_Hygiene_Score_v1)

hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

oral_v1v3$hiv_status_v1 <- factor(oral_v1v3$hiv_status_v1, levels = hiv_stat)
pdf("OHS_diff,v1v3.long.pdf")
ggplot(oral_v1v3, aes(x=factor(hiv_status_v1, levels=hiv_stat),y=diff_OHS))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./OHS_diff,v1v3.long.pdf")
```
# 7. DESEQ for V1
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
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
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "1",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  134122"
# out of 713739 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 64271, 9%
# LFC < 0 (down)     : 69851, 9.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 207567, 29%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  130752"
# out of 713739 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 79770, 11%
# LFC < 0 (down)     : 99675, 14%
# outliers [1]       : 0, 0%
# low counts [2]     : 179886, 25%
# (mean count < 2)
write.table(resLFC, file="deseq_results_V1-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_V1-HIvHUU.RData")
```
# 8. DESEQ for V2
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
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
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "2",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  115976"
# out of 691904 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 54749, 7.9%
# LFC < 0 (down)     : 61227, 8.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 241457, 35%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  115976"
# out of 691904 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 78420, 11%
# LFC < 0 (down)     : 86468, 12%
# outliers [1]       : 0, 0%
# low counts [2]     : 241457, 35%
# (mean count < 3)
write.table(resLFC, file="deseq_results_V2-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_V2-HIvHUU.RData")
```
# 9. DESEQ for V3
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
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
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HI",]
submap <- submap[submap$visit_num == "3",]
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
# [1] "number of genes with adjusted p value lower than 0.05:  92224"
# out of 554488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 45147, 8.1%
# LFC < 0 (down)     : 47077, 8.5%
# outliers [1]       : 172074, 31%
# low counts [2]     : 96749, 17%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  86021"
# out of 554488 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 49828, 9%
# LFC < 0 (down)     : 58603, 11%
# outliers [1]       : 172074, 31%
# low counts [2]     : 31241, 5.6%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_V3-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_V3-HIvHUU.RData")
```
# 10. Make tree heatmaps looking at the differences in expression across all 3 visits
## 3.1 P. ginigvalis
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
load("./deseq_results_V1-HIvHUU.RData")

tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
# virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
# virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# now make it for virulence factors
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)



pdf("pging_heat.viru_genes.pdf", width = 8, height = 8)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))
dev.off()
system("~/.iterm2/imgcat ./pging_heat.viru_genes.pdf")
```
### 3.1.1 Split by gene
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "rpoC")
for (virulence2 in virulence_genes_pg) {
load("./deseq_results_V1-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# for p. gingivalis
# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_heat_", virulence2, ".pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_heat_", virulence2 , ".pdf", sep=""))
}
```
### 3.1.2 Split by gene combined scale
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "rpoC")
for (virulence2 in virulence_genes_pg) {
load("./deseq_results_V1-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# for p. gingivalis
# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_wrap(~type) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6,  labeller = labeller(type = facet_labels)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_heat_combined_", virulence2, ".pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_heat_combined_", virulence2 , ".pdf", sep=""))
}
```
## 3.2 T. denticola
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
load("./deseq_results_V1-HIvHUU.RData")

tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_denticola/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
virulence_genes_td <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Treponema_denticola")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Treponema_denticola")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# now make it for virulence factors
virulence_genes_td <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Treponema_denticola")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Treponema_denticola")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors
virulence_genes_td <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tips,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf("tdent_heat.viru_genes.pdf", width = 12, height = 8)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.83, .9))
dev.off()
system("~/.iterm2/imgcat ./tdent_heat.viru_genes.pdf")
```
### 3.2.1 Split by gene
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_td <- c("prtP", "oppA")
for (virulence2 in virulence_genes_td) {
load("./deseq_results_V1-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_denticola/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./tdent_heat_", virulence2, ".pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.83, .9))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./tdent_heat_", virulence2 , ".pdf", sep=""))
}
```
## 3.3 T. forsythia
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
load("./deseq_results_V1-HIvHUU.RData")

tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_forsythia/align/tforsythia.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
virulence_genes_tf <- c("bspA", "clpB", "clpP", "clpX", "lipA", "fimA", "lpxA", "lpxB", "lpxD", "lpxK", "tfsA","tfsB", "fsa", "hptA" ,"hpdB", "hpdC" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Tannerella_forsythia")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Tannerella_forsythia")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_tf & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Tannerella_forsythia")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Tannerella_forsythia")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_tf & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_tf & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .008)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1.5,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf("tforsythia_heat.viru_genes.pdf", width = 10, height = 8)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.9, .6))
dev.off()
system("~/.iterm2/imgcat ./tforsythia_heat.viru_genes.pdf")
```
### 3.1.2 Split by gene combined scale
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_pg <- c("bspA", "rpoC")
for (virulence2 in virulence_genes_pg) {
load("./deseq_results_V1-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_forsythia/align/tforsythia.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# for p. gingivalis
# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .075)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_wrap(~type) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  facet_wrap(~type, ncol =6,  labeller = labeller(type = facet_labels)) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./tfor_heat_combined_", virulence2, ".pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./tfor_heat_combined_", virulence2 , ".pdf", sep=""))
}
```

# 11. DESEQ for V1 vs V3
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses")
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
submap <- metadata[metadata$visit_num == "1" | metadata$visit_num == "3",]
# submap <- submap[submap$visit_num == "1",]
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
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~visit_num)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$visit_num <- factor(star_results$visit_num, levels=c("1", "3"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  191003"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 43339, 0.66%
# LFC < 0 (down)     : 147664, 2.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 5551410, 84%
# (mean count < 1)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="visit_num_3_vs_1", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  191003"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 78646, 1.2%
# LFC < 0 (down)     : 169992, 2.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 5551410, 84%
# (mean count < 1)
write.table(resLFC, file="deseq_results-V1vV3.txt", quote=F, sep="\t")
save.image("deseq_results-V1vV3.RData")
```
Valcano
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("./deseq_results-V1vV3.RData")
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
pdf("volcano-V1vV3.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-V1vV3.pdf")
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
pdf("volcano-V1vV3.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-V1vV3.red_viru.pdf")
```
## 11. P. gingvalis
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
load("./deseq_results-V1vV3.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)
tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .035)
tip_labels <- rev(get_taxa_name(tree_plot))
# make log10 for RNA
ann <- homd[homd$tag %in% rownames(subcount),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(subcount), rownames(ann)))]
subcount <- subcount[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(subcount)==rownames(ann)) # should all return true
# if all are true, merge together
subcount <- cbind(subcount, ann)
# get list of genera to pull from rpoC data later

# collapse by species and sum across rows
merge.count <- subcount %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# long format
collapsed_long <- group_sums %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
pgings_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis",]$genome))

collapsed_long2 <- collapsed_long[collapsed_long$genome %in% pgings_seqs,]
rna_abund <- left_join(collapsed_long2, submap, by = c("sample" = "sample_id"))

summary_log10_counts <- rna_abund %>%
  group_by(visit_num, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))

#loop
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "rpoC")
for (virulence2 in virulence_genes_pg) {

# add in annotations
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_200days <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_200days <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .055)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
summary_log10_counts$metric <- "log10_count"
combined_metrics<- bind_rows(
  sig_average_baseMean_200days %>% mutate(type = "baseMean"),
  sig_average_log2_200days %>% mutate(type = "log2FoldChange"),
  summary_log10_counts %>% mutate(type = "log10(count)")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean" = "baseMean",
  "log2FoldChange" = "log2FoldChange",
  "log10(count)" = "log10(count)"
)
combined_metrics_complete[combined_metrics_complete$genome  == "SEQF3212.1",]
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_gradient2(
    low = "#5AB4AC", mid = "white", high = "#D8B365", 
    midpoint = 0,
    guide = guide_colorbar(title = "Average log fold change", title.position = "top")
  ) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "log10_count"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, 
                       guide = guide_colorbar(title = "Average log10(count)", title.position = "top")) +
  facet_wrap(~ interaction(metric, ifelse(metric == "log10_count", as.character(visit_num), as.character(type))),
             scales = "free", ncol = 4,
             labeller = labeller(
               `type` = facet_labels,
               `sampling_date` = label_value
             )) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_V1_V3_combined_", virulence2, ".pdf"), width = 10, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .5))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_V1_V3_combined_", virulence2 , ".pdf", sep=""))
}
```
# 12. DESEQ for Samples we have over 200 days
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
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

# list samples of interest
twice <- c("DM00101V2PQ16","DM00101V3PQ16","DM00305V2PQ16","DM00305V3PQ16","DM00338V2PQ16","DM00338V3PQ16","DM00008V1PQ16-2","DM00008V2PQ16","DM00023V1PQ16-1","DM00023V2PQ16","DM00035V1PQ16","DM00035V2PQ16","DM00074V1PQ55","DM00074V2PQ55","DM00142V1PQ16","DM00142V2PQ16","DM00255V1PQ26","DM00255V2PQ26","DM00266V1PQ16","DM00266V2PQ16","DM00464V1PQ65","DM00464V2PQ65","DM00563V1PQ55","DM00563V2PQ55")

# create log2 expression of red complex genes 
# filter metadata
map.twice <- submap[rownames(submap) %in% twice,]
count.twice <- subcount[,colnames(subcount) %in% twice]
# remove any genes with a sum count of zero
count.twice <- count.twice[rowSums(count.twice) != 0,]
table(colnames(count.twice)==map.twice$sample_id) # should return all true
map.twice <- map.twice[order(map.twice$study_id, map.twice$visit_num), ]

# create a new column to label sampling order
map.twice$sampling_order <- ave(map.twice$visit_num, map.twice$study_id, 
                                FUN = function(x) {
                                  order(x)
                                })
table(colnames(count.twice)==map.twice$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = count.twice, colData = map.twice, design = ~sampling_order)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$sampling_order <- factor(star_results$sampling_order, levels=c("1", "2"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  89179"
# out of 695471 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 38636, 5.6%
# LFC < 0 (down)     : 50543, 7.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 242705, 35%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="sampling_order_2_vs_1", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  84212"
# out of 695471 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 56410, 8.1%
# LFC < 0 (down)     : 70676, 10%
# outliers [1]       : 0, 0%
# low counts [2]     : 202175, 29%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results-200days.txt", quote=F, sep="\t")
save.image("deseq_results-200days.RData")
```
Valcano
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("./deseq_results-200days.RData")
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
pdf("volcano-200days.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-200days.pdf")
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
pdf("volcano-200days.red_viru.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-200days.red_viru.pdf")
```
### 12.1.1 Split by gene combined scale P. ging
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
load("./deseq_results-200days.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)
tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .035)
tip_labels <- rev(get_taxa_name(tree_plot))
# make log10 for RNA
ann <- homd[homd$tag %in% rownames(count.twice),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(count.twice), rownames(ann)))]
count.twice <- count.twice[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.twice)==rownames(ann)) # should all return true
# if all are true, merge together
count.twice <- cbind(count.twice, ann)
# get list of genera to pull from rpoC data later
count.twice.gen <- unique(count.twice$species)
# check that locus tags match between the two dataframes
table(rownames(count.twice)==rownames(ann)) # should all return true
# if all are true, merge together
count.twice <- cbind(count.twice, ann)

# collapse by species and sum across rows
merge.count <- count.twice %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# long format
collapsed_long <- group_sums %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
pgings_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis",]$genome))

collapsed_long2 <- collapsed_long[collapsed_long$genome %in% pgings_seqs,]
rna_abund <- left_join(collapsed_long2, map.twice, by = c("sample" = "sample_id"))

summary_log10_counts <- rna_abund %>%
  group_by(sampling_order, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))

# pdf("test.pdf") 
# ggplot(summary_log10_counts, aes(x = as.factor(sampling_order), y = factor(genome, levels = tip_labels), fill = mean_log10)) +
#   geom_tile() +
#   scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
#                        guide = guide_colorbar(title = "Average log10(count)", title.position = "top"))
# dev.off()
# system("~/.iterm2/imgcat ./test.pdf")

# pdf("test.pdf") 
# ggplot(rna_abund, aes(x =  factor(genome, levels = tip_labels), y = log10_count, fill = as.factor(sampling_order))) +
#   geom_bar(stat="identity", position="stack") +
#   theme_minimal() +
#   facet_wrap(~study_id)
#   # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   # scale_fill_manual(values = gencols) +
#   # geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
# dev.off()
# system("~/.iterm2/imgcat ./test.pdf")


# pdf("test.pdf") 
# ggplot(summary_log10_counts, aes(x = as.factor(sampling_order), y = median_log10, fill = genome)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = some_palette_vector) +  # Optional: define colors for genomes
#   labs(x = "Sampling Order", y = "Average log10(count)", fill = "Genome") +
#   theme_minimal()
# dev.off()
# system("~/.iterm2/imgcat ./test.pdf")

#  # doing by relative abudnance
# count.twice$SEQF <- sub("^(SEQF\\d+)-.*", "\\1", rownames(count.twice))

# seqf_counts <- count.twice %>%
#   group_by(SEQF) %>%
#   summarise(across(where(is.numeric), sum))
# total_counts <- colSums(genecounts)
# percent_seqf <- seqf_counts
# percent_seqf[,-1] <- sweep(percent_seqf[,-1], 2, total_counts, FUN = "/") * 100
# homd_reduced <- unique(homd[, c("species", "genome", "seq")])

# percent_seqf2 <- left_join(percent_seqf, homd_reduced, by = c("SEQF" = "seq"))

# pgings_seqs <- percent_seqf2[percent_seqf2$species  == "Porphyromonas_gingivalis",]
# # submap <- submap[submap$sample_id %in% sample_list, ]
# pging_count <- melt(pgings_seqs)
# rna_abund <- left_join(pging_count, map.twice, by = c("variable" = "sample_id"))

# pdf("test.pdf") 
# ggplot(rna_abund, aes(x = as.factor(sampling_order), y = factor(genome, levels = tip_labels), fill = value)) +
#   geom_tile() +
#   scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
#                        guide = guide_colorbar(title = "Average base mean", title.position = "top"))+
#   facet_wrap(~study_id)
# dev.off()
# system("~/.iterm2/imgcat ./test.pdf")

#loop
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "rpoC")
for (virulence2 in virulence_genes_pg) {

# add in annotations
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_200days <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_200days <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .055)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
summary_log10_counts$metric <- "log10_count"
combined_metrics<- bind_rows(
  sig_average_baseMean_200days %>% mutate(type = "baseMean"),
  sig_average_log2_200days %>% mutate(type = "log2FoldChange"),
  summary_log10_counts %>% mutate(type = "log10(count)")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean" = "baseMean",
  "log2FoldChange" = "log2FoldChange",
  "log10(count)" = "log10(count)"
)
combined_metrics_complete[combined_metrics_complete$genome  == "SEQF3212.1",]
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_gradient2(
    low = "#5AB4AC", mid = "white", high = "#D8B365", 
    midpoint = 0,
    guide = guide_colorbar(title = "Average log fold change", title.position = "top")
  ) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "log10_count"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, 
                       guide = guide_colorbar(title = "Average log10(count)", title.position = "top")) +
  facet_wrap(~ interaction(metric, ifelse(metric == "log10_count", as.character(sampling_order), as.character(type))),
             scales = "free", ncol = 4,
             labeller = labeller(
               `type` = facet_labels,
               `sampling_date` = label_value
             )) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_200days_combined_", virulence2, ".pdf"), width = 10, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .5))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_200days_combined_", virulence2 , ".pdf", sep=""))
}
```
### 12.1.1 Split by gene combined scale T. dent
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
load("./deseq_results-200days.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_denticola/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)
tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .085)
tip_labels <- rev(get_taxa_name(tree_plot))
# make log10 for RNA
ann <- homd[homd$tag %in% rownames(count.twice),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(count.twice), rownames(ann)))]
count.twice <- count.twice[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.twice)==rownames(ann)) # should all return true
# if all are true, merge together
count.twice <- cbind(count.twice, ann)
# get list of genera to pull from rpoC data later
count.twice.gen <- unique(count.twice$species)
# check that locus tags match between the two dataframes

# collapse by species and sum across rows
merge.count <- count.twice %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# long format
collapsed_long <- group_sums %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long1 <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
pgings_seqs <- unique(sort(homd[homd$species  == "Treponema_denticola",]$genome))

collapsed_long2 <- collapsed_long1[collapsed_long1$genome %in% pgings_seqs,]
rna_abund <- left_join(collapsed_long2, map.twice, by = c("sample" = "sample_id"))

summary_log10_counts <- rna_abund %>%
  group_by(sampling_order, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))


#loop
virulence_genes_pg <- c("prtP", "prcB", "prcA", "oppA", "rpoC")
for (virulence2 in virulence_genes_pg) {

# add in annotations
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "") 

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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_200days <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_200days <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")

#base mean
summary_log10_counts$metric <- "log10_count"
combined_metrics<- bind_rows(
  sig_average_baseMean_200days %>% mutate(type = "baseMean"),
  sig_average_log2_200days %>% mutate(type = "log2FoldChange"),
  summary_log10_counts %>% mutate(type = "log10(count)")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean" = "baseMean",
  "log2FoldChange" = "log2FoldChange",
  "log10(count)" = "log10(count)"
)
combined_metrics_complete[combined_metrics_complete$genome  == "SEQF3212.1",]
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_gradient2(
    low = "#5AB4AC", mid = "white", high = "#D8B365", 
    midpoint = 0,
    guide = guide_colorbar(title = "Average log fold change", title.position = "top")
  ) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "log10_count"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, 
                       guide = guide_colorbar(title = "Average log10(count)", title.position = "top")) +
  facet_wrap(~ interaction(metric, ifelse(metric == "log10_count", as.character(sampling_order), as.character(type))),
             scales = "free", ncol = 4,
             labeller = labeller(
               `type` = facet_labels,
               `sampling_date` = label_value
             )) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./tdent_200days_combined_", virulence2, ".pdf"), width = 10, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.45, .5))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./tdent_200days_combined_", virulence2 , ".pdf", sep=""))
}
```
# 13. DESEQ for Samples we have over 400 days
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
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

# list samples of interest
triple <- c("DM00044V1PQ16","DM00044V2PQ16","DM00044V3PQ16","DM00103V1PQ16","DM00103V2PQ16","DM00103V3PQ16","DM00254V1PQ55","DM00254V2PQ55","DM00254V3PQ55","DM00519V1PQ16","DM00519V2PQ16","DM00519V3PQ16")

# create log2 expression of red complex genes 
# filter metadata
map.triple <- submap[rownames(submap) %in% triple,]
count.triple <- subcount[,colnames(subcount) %in% triple]
# remove any genes with a sum count of zero
count.triple <- count.triple[rowSums(count.triple) != 0,]
table(colnames(count.triple)==map.triple$sample_id) # should return all true
map.triple <- map.triple[order(map.triple$study_id, map.triple$visit_num), ]

# create a new column to label sampling order
map.triple$sampling_order <- ave(map.triple$visit_num, map.triple$study_id, 
                                FUN = function(x) {
                                  order(x)
                                })
table(colnames(count.triple)==map.triple$sample_id) # should return all true

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = count.triple, colData = map.triple, design = ~sampling_order)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$sampling_order <- factor(star_results$sampling_order, levels=c("1", "2"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  89179"
# out of 695471 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 38636, 5.6%
# LFC < 0 (down)     : 50543, 7.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 242705, 35%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="sampling_order_2_vs_1", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  84212"
# out of 695471 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 56410, 8.1%
# LFC < 0 (down)     : 70676, 10%
# outliers [1]       : 0, 0%
# low counts [2]     : 202175, 29%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results-400days.txt", quote=F, sep="\t")
save.image("deseq_results-400days.RData")
```
Valcano
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("~/rna_dohmain/11-perio/10-longitudinal_analyses")
load("./deseq_results-200days.RData")
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
pdf("volcano-200days.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-200days.pdf")
```

