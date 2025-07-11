# 1. Get proportion of transcripts that belong to red complex
```sh
cp annotations.merge.txt ../09-urease/09-global-distro
paste annotations.merge.txt <(awk '{print $4, $5}' annotations.merge.txt | sed 's/ /_/') | sed 's/ /\t/g' > temp
mv temp annotations.merge.txt
```
```py
import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/homd_map/read_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.merge.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus_Species


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('species_reads.txt', sep="\t", index=True, header = True)
```
# 1. DNA vs RNA comparison
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio/08-DNAvRNA")
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
rna_counts <- read.csv("~/rna_dohmain/09-urease/09-global-distro/species_rpoC.txt", sep='\t')
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
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor(aes(color = hiv_status))+
 facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.RNAvDNA.corr.pdf")