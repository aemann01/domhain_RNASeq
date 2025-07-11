# Look at output from DNA rpoC
```R
library(phyloseq)
library(ggplot2)
setwd("~/rna_dohmain/11-perio/07-orange-complex")
#get relative abundance of p. gingivalis
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Campylobacter_gracilis" | V8=="Campylobacter _rectus" | V8=="Campylobacter_showae" | V8=="Eubacterium_nodatum" | V8=="Fusobacterium_nucleatum" | V8=="Fusobacterium_periodonticum" | V8=="Prevotella_micros" | V8=="Prevotella_intermedia" | V8=="Prevotella_nigrescens"| V8=="Streptococcus_constellatus")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("orange.dna.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./orange.dna.pdf")
test <- data %>% group_by(Sample) %>% summarise(average_baseMean = across(c(Abundance)))
glom <- tax_glom(actino, taxrank=rank_names(actino)[6])
data <- psmelt(glom) # create dataframe from phyloseq object
mean(data$Abundance)
```
# Look at output from RNA rpoC
```R
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
meta <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
counts <- read.csv("../../09-urease/09-global-distro/species_rpoC.txt", header=T,sep = "\t")
data <- melt(counts)
data <- left_join(data, meta, by = join_by(sample ==  sample_id))
data$genus <- gsub(x = data$variable, pattern = "_.*", replacement = "")
#actinomyces naeslundii plot
sub_data <- data[data$variable == "Campylobacter_gracilis" | data$variable=="Campylobacter _rectus" | data$variable=="Campylobacter_showae" | data$variable=="Eubacterium_nodatum" | data$variable=="Fusobacterium_nucleatum" | data$variable=="Fusobacterium_periodonticum" | data$variable=="Prevotella_micros" | data$variable=="Prevotella_intermedia" | data$variable=="Prevotella_nigrescens"| data$variable=="Streptococcus_constellatus",]
pdf("orange.rna.rpoc.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./orange.rna.rpoc.pdf")
