# Get the data frames
all
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
setwd("/home/suzanne/rna_dohmain/09-urease/03-diff-abundance")
# add in annotations
resLFC <- read.table("deseq_results_operon-HvD.txt",)

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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

dfDall <- subset(sortdf, log2FoldChange < 0)
dfHall <- subset(sortdf, log2FoldChange > 0)
```
HUU
```R
resLFC <- read.table("deseq_results_operon-HvD-HUU.txt")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

dfDHUU <- subset(sortdf, log2FoldChange < 0)
dfHHUU <- subset(sortdf, log2FoldChange > 0)
```
HEU
```R
resLFC <- read.table("deseq_results_operon-HvD-HEU.txt")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

dfDHEU <- subset(sortdf, log2FoldChange < 0)
dfHHEU <- subset(sortdf, log2FoldChange > 0)
```
HI
```R
resLFC <- read.table("deseq_results_operon-HvD-HI.txt")
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
res_ord$GeneInfo <- paste(res_ord$Species,res_ord$gene)
# finally get a list of the top 10 genes (based on highest p value) to label on the volcano plot
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$locus_tag, 5)
# negative top 10
top <- tail(sortdf$locus_tag, 5)
# concatenate
labgenes <- c(top, low)

dfDHI <- subset(sortdf, log2FoldChange < 0)
dfHHI <- subset(sortdf, log2FoldChange > 0)
```
# 2. Make bubble plot
```R
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)

# rename the input dataframes so that you don't have to run all the code above
dfHall -> all
dfHHI -> HI
dfHHEU -> HEU
dfHHUU -> HUU 

# subset for clarity
all <- all %>% select('baseMean', 'Genus', 'Species', "SEQ_ID")
HI <- HI %>% select('baseMean', 'Genus', 'Species', "SEQ_ID")
HEU <- HEU %>% select('baseMean', 'Genus', 'Species', "SEQ_ID")
HUU <- HUU %>% select('baseMean', 'Genus', 'Species', "SEQ_ID")

# group and get average base mean
group_all <- all %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HI <- HI %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_all <- group_all %>% mutate(Source = "All")
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged <- rbind(group_all, group_HI, group_HEU, group_HUU)
# add health indicator
merged <- merged %>% mutate(Condition = "H")

# do the same thing with the disease results
dfDall -> all
dfDHI -> HI
dfDHEU -> HEU
dfDHUU -> HUU 

# subset for clarity
all <- all %>% select('baseMean', 'Genus', 'Species', 'SEQ_ID')
HI <- HI %>% select('baseMean', 'Genus', 'Species', 'SEQ_ID')
HEU <- HEU %>% select('baseMean', 'Genus', 'Species', 'SEQ_ID')
HUU <- HUU %>% select('baseMean', 'Genus', 'Species', 'SEQ_ID')

# group and get average base mean
group_all <- all %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HI <- HI %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HEU <- HEU %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))
group_HUU <- HUU %>% group_by(Genus, Species, SEQ_ID) %>% summarise(average_baseMean = mean(baseMean))

# add source column so you can track them back
group_all <- group_all %>% mutate(Source = "All")
group_HI <- group_HI %>% mutate(Source = "HI")
group_HEU <- group_HEU %>% mutate(Source = "HEU")
group_HUU <- group_HUU %>% mutate(Source = "HUU")

# merge all into one
merged2 <- rbind(group_all, group_HI, group_HEU, group_HUU)
# add health indicator
merged2<- merged2 %>% mutate(Condition = "D")
# merge them both together!
df <- rbind(merged, merged2)

# now that the data is formatted properly, can make a bubble plot
# Count number of species per genus
species_counts <- df %>%
  group_by(Genus) %>%
  summarise(num_species = n_distinct(Species)) %>%
  arrange(desc(num_species))  # Order by number of species descending

# Reorder Genus based on the number of species in each group
df$Genus <- factor(df$Genus, levels = species_counts$Genus)
df$Species <- as.factor(df$Species)
# Order dataframe based on Genus factor levels
df <- df[order(df$Genus), ]
# Create a nested label combining Genus and Species
df$Species_nested <- paste(df$Genus, df$Species, df$SEQ_ID, sep = "_")
# get order of nested species
ordered_species <- rev(unique(df$Species_nested))
df$Species_nested <- as.factor(df$Species_nested)
df$Species_nested <- factor(df$Species_nested, levels = ordered_species)
# set levels of the x axis as well
df$Source <- factor(df$Source, levels = c("All", "HUU", "HEU", "HI"))
# and order of grid
df$Condition <- factor(df$Condition, levels = c("H", "D"))
viridis(4)

order <- rev(sort(unique(as.character(df$Species_nested))))
df$Species_nested <- ordered(df$Species_nested, levels=order)

pdf("HvD_bubble_plot.pdf", width = 10)
  ggplot(df, aes(x = Source, y = Species_nested, size = average_baseMean, color = Genus)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  labs(
    x = "Source",
    y = "Species",
    size = "Base Mean"
  ) +
  facet_grid(. ~ Condition, switch = "y") +
  scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")) +
  scale_y_discrete("Species_nested")+
   theme_minimal()
dev.off()
system("~/.iterm2/imgcat HvD_bubble_plot.pdf")
```
