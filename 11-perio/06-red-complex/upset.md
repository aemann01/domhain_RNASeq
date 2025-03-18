# 1. Upset Plot
```R
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ggplot2)
library(ComplexUpset)


setwd("~/rna_dohmain/11-perio/06-red-complex")
HI <- read.csv("deseq_results_red-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HI),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HI), rownames(ann)))]
HI <- HI[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HI)==rownames(ann)) # should all return true
# if all are true, merge together
HI <- cbind(HI, ann)
#HEU
HEU$hiv_status <- "HEU"
#get annots for HEEU
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HEU),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HEU), rownames(ann)))]
HEU <- HEU[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HEU)==rownames(ann)) # should all return true
# if all are true, merge together
HEU <- cbind(HEU, ann)

#bind together HI and HEU
resdf <- rbind(HI, HEU)

#get only signficant genes
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of -2 since that is the side from their respective deseq analysis
resdf <- resdf %>%
  mutate(mark = if_else(log2FoldChange <= -lfc & padj <= pval, "1", "0"))
resdf$mark <- as.numeric(resdf$mark)
resdf_filtered <- select(resdf, mark, hiv_status, gene_tag, species)
resdf_wide <- resdf_filtered %>%
  pivot_wider(names_from = hiv_status, values_from = mark, values_fill = list(mark = 0)) %>%
  arrange(gene_tag)
#make upset plot
resdf_df <- as.data.frame(resdf_wide)
resdf_clean_df <- resdf_df %>%
  filter(!is.na(HI) & !is.na(HEU))
resdf_clean_df$species <- as.factor(resdf_clean_df$species)

pdf("HI_HEU.upset.pdf")
upset(resdf_clean_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially Expressed Genes", sets.x.label="Total Number of Differentially Expressed Genes",
	queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()
system("~/.iterm2/imgcat ./HI_HEU.upset.pdf")
```
# 2. Donut Plot
Reads belonging to each group
```R
library(tidyverse)
library(reshape2)
library(webr)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
#get annots for genes
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(genecounts),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(genecounts), rownames(ann)))]
genecounts <- HEU[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(genecounts)==rownames(ann)) # should all return true
# if all are true, merge together
genecounts <- cbind(genecounts, ann)
gene_long <- melt(genecounts)
combined_data <- merge(gene_long, metadata, by.x = "variable", by.y = "sample_id", all.x = TRUE)

#make donut plots
gene_sum <- combined_data %>% 
  group_by(hiv_status, species) %>% 
  summarise(Total = sum(value, na.rm = TRUE))

pdf("red_rna.donut.pdf")
PieDonut(gene_sum, aes("hiv_status", "species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./red_rna.donut.pdf")
```
DNA donut
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio/06-red-complex")
#get relative abundance of dna
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, Abundance, V8)
#make donut plots
gene_sum <- red_dna %>% 
  group_by(hiv_status, V8) %>% 
  summarise(Total = sum(Abundance, na.rm = TRUE))
pdf("red_dna.donut.pdf")
PieDonut(gene_sum, aes("hiv_status", "V8", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./red_dna.donut.pdf")  
```
