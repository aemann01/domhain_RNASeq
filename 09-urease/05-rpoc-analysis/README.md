# 1. Load libraries
```R
library(philr, warn.conflicts = F, quietly = T)
library(RColorBrewer, warn.conflicts = F, quietly = T)
library(UpSetR, warn.conflicts = F, quietly = T)
library(ggfortify, warn.conflicts = F, quietly = T)
library(randomForest, warn.conflicts = F, quietly = T)
# library(rfUtilities, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(gridExtra, warn.conflicts = F, quietly = T)
library(microbiome, warn.conflicts = F, quietly = T)
library(phylofactor, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(pairwiseAdonis, warn.conflicts = F, quietly = T)
library(ape, warn.conflicts = F, quietly = T)
library(metagMisc, warn.conflicts = F, quietly = T)
library(ranacapa, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)
library(ggdendro, warn.conflicts = F, quietly = T)
#set working directory
setwd("~/rna_dohmain/09-urease/05-rpoc-analysis")
```
# 2. Make phyloseq object
```R
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
tree <- read.tree("RAxML_bestTree.ref.tre")
tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
```
# 3. Filter out low frequency ASVs
```R
# filter out samples with fewer than 1000 reads (based on ASV rareness, this shouldn't be an issue)
ps.dat <- prune_samples(sample_sums(ps.dat) > 1000, ps.dat)
ps.dat
```
# 4. Proportion of different A. naeslundii
```R
#convert to relative abundance
rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V8=="Actinomyces_naeslundii")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("actino.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./actino.pdf")
```
# 5. Proportion of different H. parainfluenzae
```R
#convert to relative abundance
rel <- microbiome::transform(ps.dat, "compositional")
haem <- subset_taxa(rel, V8=="Haemophilus_parainfluenzae")
glom <- tax_glom(haem, taxrank=rank_names(haem)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("h_para.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./h_para.pdf")

```