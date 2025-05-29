# 1. Look at HI for upregulated across genome
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
load("../06-red-complex/deseq_results_red-HIvHUU.RData")
# add in annotations
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
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

sigloc$GeneInfo <- paste(sigloc$species,sigloc$gene)
sigloc$gene_name <- gsub(x = sigloc$gene, pattern = "_.", replacement = "") 
sig_low <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange <= -lfc) 

count_sig <- sig_low %>%
  group_by(contig, species) %>%
  summarise(Total_sig = sum(!is.na(tag)), .groups = "drop")
count_genome <- homd %>%
  group_by(contig, species) %>%
  summarise(Total = sum(!is.na(tag)), .groups = "drop")

# combine
combined_data <- count_genome %>%
  left_join(count_sig, by = c("contig", "species"))
combined_data$Total_sig <- ifelse(is.na(combined_data$Total_sig), 0, combined_data$Total_sig)
combined_data$prop <-  combined_data$Total_sig / combined_data$Total
combined_data
sorted_combined_data <- combined_data %>%
  arrange(desc(prop))
print(sorted_combined_data, n = 52)

specCols<-  c("Porphyromonas_gingivalis" = "#340043", "Treponema_denticola" = "#FBE51F", "Tannerella_forsythia" = "#1E7F7A")
pdf("HIvHUU.density_coverage.pdf")
ggplot(sorted_combined_data, aes(x = prop, fill =species)) +
  geom_density( alpha = 0.5) +
  labs(title = "Density Plot of Upregulated Genes",
       x = "Proportion of Genes",
       y = "Density") +
  theme_minimal()+
  scale_fill_manual(values = specCols)
dev.off()
system("~/.iterm2/imgcat ./HIvHUU.density_coverage.pdf")

```
# 2. Look at HEU for upregulated across genome
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
load("../06-red-complex/deseq_results_red-HEUvHUU.RData")
# add in annotations
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
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

sigloc$GeneInfo <- paste(sigloc$species,sigloc$gene)
sigloc$gene_name <- gsub(x = sigloc$gene, pattern = "_.", replacement = "") 
sig_low <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange <= -lfc) 

count_sig <- sig_low %>%
  group_by(contig, species) %>%
  summarise(Total_sig = sum(!is.na(tag)), .groups = "drop")
count_genome <- homd %>%
  group_by(contig, species) %>%
  summarise(Total = sum(!is.na(tag)), .groups = "drop")

# combine
combined_data <- count_genome %>%
  left_join(count_sig, by = c("contig", "species"))
combined_data$Total_sig <- ifelse(is.na(combined_data$Total_sig), 0, combined_data$Total_sig)
combined_data$prop <-  combined_data$Total_sig / combined_data$Total  
sorted_combined_data <- combined_data %>%
  arrange(prop)
print(sorted_combined_data, n = 52)

specCols<-  c("Porphyromonas_gingivalis" = "#340043", "Treponema_denticola" = "#FBE51F", "Tannerella_forsythia" = "#1E7F7A")
pdf("HEUvHUU.density_coverage.pdf")
ggplot(sorted_combined_data, aes(x = prop, fill =species)) +
  geom_density( alpha = 0.5) +
  labs(title = "Density Plot of Upregulated Genes",
       x = "Proportion of Genes",
       y = "Density") +
  theme_minimal()+
  scale_fill_manual(values = specCols)
dev.off()
system("~/.iterm2/imgcat ./HEUvHUU.density_coverage.pdf")
```
# 3. How many genes have reads mapped to them
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
load("../06-red-complex/deseq_results_red-HIvHUU.RData")
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
# add annotations to gene counts
genecounts$gene_name <- row.names(genecounts)
annot_gene <- left_join(genecounts,homd, by = c("gene_name"= "tag"))
# melt dataframe
long_counts <- melt(annot_gene)
long_counts <- left_join(long_counts,metadata, by = c("variable"= "sample_id"))

# precent mapped by sample
percent_genes_mapped <- long_counts %>%
  group_by(contig, variable, species, hiv_status) %>%
  summarise(
    total_genes = n(),
    genes_with_reads = sum(value > 0),  
    percent_mapped = (genes_with_reads / total_genes) * 100 
  )
print(percent_genes_mapped$percent_mapped)
sort_percent_genes_mapped <- percent_genes_mapped %>%
  arrange(desc(percent_mapped))
print(sort_percent_genes_mapped, n =200)

# highest for each species in each sample
highest_per_species_sample <- sort_percent_genes_mapped %>%
  group_by(species, variable) %>% 
  slice_max(percent_mapped, n = 1)
sort_highest_per_species_sample <- highest_per_species_sample %>%
  arrange(desc(percent_mapped))
print(sort_highest_per_species_sample, n =10)

# now find by hiv status mean
percent_genes_mean <- long_counts %>%
  group_by(variable, species, hiv_status) %>%
  summarise(
    total_genes = n(),
    genes_with_reads = sum(value > 0),  
    percent_mapped = (genes_with_reads / total_genes) * 100 
  )
sort_percent_genes_mean <- percent_genes_mean %>%
  arrange(desc(percent_mapped))
print(sort_percent_genes_mean, n =200)

mean_hiv_status <- percent_genes_mean %>%
  group_by(species, hiv_status) %>%
  summarise(mean_mapped = mean(percent_mapped))
mean_hiv_status
sort_percent_genes_mean <- percent_genes_mean %>%
  arrange(desc(percent_mapped))
print(sort_percent_genes_mean, n =93)


percent_genes_mean <- percent_genes_mapped %>%
  group_by(contig, species, hiv_status) %>%
  summarise(
    total_geness = mean(percent_mapped) 
  ) 
print(percent_genes_mean)
print(percent_genes_mean, n=156)
```
# 4. Look at precent of upregulated genes vs genome that has gene mapped
```R
load("../06-red-complex/deseq_results_red-HIvHUU.RData")
# add in annotations
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
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

sigloc$GeneInfo <- paste(sigloc$species,sigloc$gene)
sigloc$gene_name <- gsub(x = sigloc$gene, pattern = "_.", replacement = "") 
sig_low <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange <= -lfc) 
# get sig genes
sig_genes <- row.names(sig_low)

genecounts$gene_name <- row.names(genecounts)
annot_gene <- left_join(genecounts,homd, by = c("gene_name"= "tag"))
# melt dataframe
long_counts <- melt(annot_gene)
long_counts <- left_join(long_counts,metadata, by = c("variable"= "sample_id"))
#total number of genes per genome that has one value
genes_with_one_count <- long_counts %>%
  group_by(variable, contig, hiv_status, species) %>%
  summarise(genes_with_count = sum(value > 0))
#total number of sig genes per genome that has one value
significant_genes_in_sample <- long_counts %>%
  filter(gene_name %in% sig_genes) %>%
  group_by(variable, contig, hiv_status) %>%
  summarise(significant_genes_count = sum(value > 0))
significant_genes_in_sample[is.na(significant_genes_in_sample)] <- 0

prop_sig <- left_join(genes_with_one_count, significant_genes_in_sample, by = c("variable", "contig", "hiv_status"))
prop_sig$ratio <- prop_sig$significant_genes_count / prop_sig$genes_with_count
prop_sig %>% arrange(desc(ratio))
prop_sig[is.na(prop_sig)] <- 0

prop_sig_species <- prop_sig %>%
  group_by(variable, species, hiv_status) %>%
  summarise(ratio_mean = mean(ratio))

specCols<- c("Porphyromonas_gingivalis" = "#340043", "Treponema_denticola" = "#FBE51F", "Tannerella_forsythia" = "#1E7F7A")
mat_colors <- c("HI" = "#8213A0", "HEU" = "#FA78FA", "HUU" = "#40A0FA")
prop_sig_species$hiv_status <- factor(prop_sig_species$hiv_status, levels=c("HI", "HEU","HUU"))
pdf("species.density_coverage.pdf", width = 7)
ggplot(prop_sig_species, aes(x = ratio_mean, fill =hiv_status)) +
  geom_density( alpha = 0.5, position = "stack") +
  labs(title = "Density Plot of Upregulated Genes",
       x = "Proportion of Upregulated Genes",
       y = "Density") +
  theme_minimal()+
  scale_fill_manual(values = mat_colors)+
  facet_wrap(~species, , ncol = 1, scales = "free")+
  theme(strip.position = "top")
dev.off()
system("~/.iterm2/imgcat ./species.density_coverage.pdf")
prop_sig %>%
  group_by( species, hiv_status) %>%
  summarise(ratio_mean = mean(ratio))
```
# 5. Look at which genes are found in most sanples
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(UpSetR)
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
load("../06-red-complex/deseq_results_red-HIvHUU.RData")
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
# add annotations to gene counts
genecounts$gene_name <- row.names(genecounts)
annot_gene <- left_join(genecounts,homd, by = c("gene_name"= "tag"))
# melt dataframe
long_counts <- melt(annot_gene)
long_counts <- left_join(long_counts,metadata, by = c("variable"= "sample_id"))

# get the number of samples for each hiv category 
total_samples_per_hiv <- long_counts %>%
  group_by(hiv_status) %>%
  summarize(total_samples = n_distinct(variable))

# count non-zero occurrences of each gene_name per HIV category
gene_counts_per_sample <- long_counts %>%
  group_by(hiv_status, contig, gene_name, species, gene) %>%
  summarize(count_non_zero = sum(value > 0)) %>%
  ungroup()

# calculate the treshold for each sample
threshold_50 <- gene_counts_per_sample %>%
  inner_join(total_samples_per_hiv, by = "hiv_status") %>%
  mutate(threshold = total_samples * 0.50)  # 50%

# #do it for over all 50$
# gene_counts_per_sample <- long_counts %>%
#   group_by(contig, gene_name, species) %>%
#   summarize(count_non_zero = sum(value > 0)) %>%
#   ungroup()
# gene_counts_per_sample$threshold <- 93*.5
# genes_in_50_percent <- gene_counts_per_sample %>%
#   filter(count_non_zero >= threshold) %>%
#   select(gene_name, contig, species)
# gene_50 <- unique(genes_in_50_percent$gene_name)

# filter by threshold
genes_in_50_percent <- threshold_50 %>%
  filter(count_non_zero >= threshold) %>%
  select(hiv_status, gene_name, contig, species, gene)
gene_50 <- unique(genes_in_50_percent$gene_name)
#total number of core genes per genome that has one value
core_gene <- genes_in_50_percent %>%
  group_by(contig, hiv_status, species, gene) %>%
  summarise(core_genes_count = n_distinct(gene_name))

count_genome <- homd %>%
  group_by(contig, species) %>%
  summarise(total_genes = sum(!is.na(tag)), .groups = "drop")
# find proprotion
prop_core <- left_join(core_gene, count_genome, by = c("contig", "species"))
prop_core$proportion_core <- prop_core$core_genes_count / prop_core$total_genes
#how many core genes per species and hiv status
prop_core %>%
  group_by(species, hiv_status) %>%
  summarise(num_genes = sum(core_genes_count))
specCols<-  c("Porphyromonas_gingivalis" = "#340043", "Treponema_denticola" = "#FBE51F", "Tannerella_forsythia" = "#1E7F7A")
pdf("species_core.density_coverage.pdf")
ggplot(prop_core, aes(x = proportion_core, fill =species)) +
  geom_density( alpha = 0.5) +
  labs(title = "Density Plot of Core Genes",
       x = "Proportion of Core Genes",
       y = "Density") +
  theme_minimal()+
  facet_wrap(~hiv_status)+
  scale_fill_manual(values = specCols)
dev.off()
system("~/.iterm2/imgcat ./species_core.density_coverage.pdf")

# make upset plot to show shared core
genes_in_50_percent$gene_tag <- genes_in_50_percent$gene_name
binary_core <- genes_in_50_percent %>%
  group_by(gene_name, gene, species) %>%
  pivot_wider(names_from = hiv_status, 
              values_from = gene_name,
              values_fn = list(gene_name = length), 
              values_fill = list(gene_name = 0)) %>%
  select(gene_tag, species, HUU, HEU, HI) %>% as.data.frame()

pdf("core_genes.upset.pdf")
upset(binary_core, order.by="freq", sets = c("HI", "HEU", "HUU"), mainbar.y.label="Number of Shared Core Genes", sets.x.label="Total Number of Core Genes",
  queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()
system("~/.iterm2/imgcat ./core_genes.upset.pdf")
binary_core %>%
  filter(HUU == 1, HEU == 1, HI == 1)
# count how many unique genes
gene_count_per_hiv <- genes_in_50_percent %>%
  group_by(hiv_status) %>%
  summarize(genes_in_50_percent = n_distinct(gene_name))
print(gene_count_per_hiv)

```
# 4. Look at normalized distro of transcript for genome by sample
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(reshape2)
library(ggpubr)
library(ggside)
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../06-red-complex/red_counts.txt", header=T, sep="\t", row.names=1)
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

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~hiv_status)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HEU","HUU"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
#pcoa diversity
vld <- varianceStabilizingTransformation(se_star)
vld_df <- assay(vld)
vld_long <- melt(vld_df)
# add annotations
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
tags_not_in_vld_long <- setdiff(homd$tag, vld_long$Var1)
tags_not_in_vld_long
merged_data <- merge(vld_long, homd, by.x = "Var1", by.y = "tag")
merged_df2 <- merge(merged_data, metadata, by.x = "Var2", by.y = "sample_id")
red_mean <- merged_df2 %>% 
  group_by(Var2, species, Var1, hiv_status) %>% 
  summarise(Total = mean(value, na.rm = TRUE))

#get just one species
sub_por <- red_mean %>% filter(species == "Porphyromonas_gingivalis")

pdf("por.distro.pdf", width = 50, height = 50)
ggplot(sub_por, aes(x = Total, fill = hiv_status)) + 
  geom_histogram(binwidth = 0.5, position = "dodge") + 
  facet_wrap(species ~ Var2, scales = "free_y") +  # Separate by contig
  theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Histograms for each Contig", 
       x = "Value", 
       y = "Frequency") +
  theme(legend.position = "top")
dev.off()
system("~/.iterm2/imgcat ./por.distro.pdf")

sub_trep <- red_mean %>% filter(species == "Treponema_denticola")

pdf("Treponema_denticola.distro.pdf", width = 50, height = 50)
ggplot(sub_trep, aes(x = Total, fill = hiv_status)) + 
  geom_histogram(binwidth = 0.5, position = "dodge") + 
  facet_wrap(species ~ Var2, scales = "free_y") +  # Separate by contig
  theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Histograms for each Contig", 
       x = "Value", 
       y = "Frequency") +
  theme(legend.position = "top")
dev.off()
system("~/.iterm2/imgcat ./Treponema_denticola.distro.pdf")

sub_tan <- red_mean %>% filter(species == "Tannerella_forsythia")

pdf("Tannerella_forsythia.distro.pdf", width = 50, height = 50)
ggplot(sub_tan, aes(x = Total, fill = hiv_status)) + 
  geom_histogram(binwidth = 0.5, position = "dodge") + 
  facet_wrap(species ~ Var2, scales = "free_y") +  # Separate by contig
  theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Histograms for each Contig", 
       x = "Value", 
       y = "Frequency") +
  theme(legend.position = "top")
dev.off()
system("~/.iterm2/imgcat ./Tannerella_forsythia.distro.pdf")

pdf("test.pdf", width = 100, height = 100)
ggplot(red_mean, aes(x = Total, fill = hiv_status)) + 
  geom_density(position = "dodge", alpha = 0.6) +  # Density plot with dodge position and transparency
  facet_wrap(species ~ Var2, scales = "free") +  # Free x-axis scale for each plot
  theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold", hjust = 0.5),  # Make species label bold
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  labs(title = "Density Plots for Each Contig", 
       x = "Value", 
       y = "Density") +
  theme(legend.position = "top")
dev.off()
system("~/.iterm2/imgcat ./test.pdf")
```
# 5. How many reads per species in each sample
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(reshape2)
library(ggpubr)
library(ggside)
setwd("/home/suzanne/rna_dohmain/11-perio/12-genome-coverage")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../06-red-complex/red_counts.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# add annotations
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
# add metadata
tags_not_in_gene_count <- setdiff(homd$tag, row.names(genecounts))
tags_not_in_gene_count
merged_data <- melt(merge(genecounts, homd, by.x = "row.names", by.y = "tag"))
merged_df2 <- merge(merged_data, metadata, by.x = "variable", by.y = "sample_id")

#sum by sample and species
red_mean <- merged_df2 %>% 
  group_by(variable, species, hiv_status) %>% 
  summarise(Total = sum(value, na.rm = TRUE))
print(red_mean, n =279)

red_mean2 <- group_by(red_mean, hiv_status, species) %>% 
  summarise(Total = mean(Total, na.rm = TRUE))
red_mean2

filtered_data <- red_mean %>%
  filter(species == "Porphyromonas_gingivalis" & Total > 5000 & hiv_status == "HI")
count_samples <- nrow(filtered_data)
print(count_samples)
