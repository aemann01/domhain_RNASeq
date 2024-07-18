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
test <- data %>% group_by(Sample) %>% summarise(average_baseMean = across(c(Abundance)))
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
# 6. urease compptent bacteria compared to RNA
Get the ASVs with bacteia that haveu urease
```sh
awk '{print $1}' ../01-operon-identify/ureABC.operon | sed '1d' >seqs 
parallel -a seqs -j 7 -k  "grep '{}' -w ../../homd_map/ure.annotations.txt" | awk '{print $4, $5}' | sed 's/ /_/' | sort | uniq > ure_taxa
parallel -a ure_taxa -j 7 -k  "grep '{}' ../../rpoc/taxonomy.txt" | sed 's/\t.*//' > ure_ASVs
```
Get ureABC locus tags
```sh
cat <(awk '{print $2}' ../01-operon-identify/ureABC.operon | sed '1d') <(awk '{print $3}' ../01-operon-identify/ureABC.operon | sed '1d') <(awk '{print $4}' ../01-operon-identify/ureABC.operon | sed '1d') > ureABC_seqs
```
Make plot
```R
library(tidyverse)
library(reshape2)
setwd("~/rna_dohmain/09-urease/05-rpoc-analysis")
seqtab_rna <- read.table("~/rna_dohmain/homd_map/read_counts.txt", header=T, row.names=1)
seqtab_rpoc <- read.table("~/rna_dohmain/rpoc/sequence_table.merged.txt", header=T, row.names=1)

seqtab_rpoc$sample <- row.names(seqtab_rpoc)
rpoc_df <- melt(as.data.frame(seqtab_rpoc))

seqtab_rna$SEQ <- row.names(seqtab_rna)
rna_df <- melt(seqtab_rna)

ure_taxa <- as.character(read.table("ure_ASVs", header = FALSE)$V1)
ure_seqs <- as.character(read.table("ureABC_seqs", header = FALSE)$V1)

# Annotating urease presence in rpoc_df
rpoc_df$ure_comp <- ifelse(rpoc_df$variable %in% ure_taxa, "yes", "no")
rpoc_pres <- rpoc_df %>% group_by(sample, ure_comp) %>% summarise(across(c(value), sum))

# Annotating urease presence in rna
rna_df$ure_comp <- ifelse(rna_df$SEQ %in% rna_df, "yes", "no")
rna_pres <- rna_df %>% group_by(sample, ure_comp) %>% summarise(across(c(value), sum))
save.image("test.RData")
```
# 7. See if A. naeslundii is picked uo by primers
```sh
wget https://github.com/egonozer/in_silico_pcr/raw/master/in_silico_PCR.pl
grep -w 1655 ~/rna_dohmain/rpoc/database/rpoc_ref.fa -A 1 > a_naeslundii.fa

perl in_silico_PCR.pl -s ~/rna_dohmain/rpoc/database/rpoc_ref.fa -a MAYGARAARMGNATGYTNCARGA -b GMCATYTGRTCNCCRTCRAA > results.txt 2> amplicons.fasta

perl in_silico_PCR.pl -s a_naeslundii.fa -a MAYGARAARMGNATGYTNCARGA -b GMCATYTGRTCNCCRTCRAA > results.txt 2> amplicons.fasta

#make a tree
grep Actinomyces -w ../../rpoc/taxonomy_bac.txt | awk '{print $1}' | while read line; 
seqtk subseq ../../rpoc/rep_set.fa <(grep Actinomyces -w ../../rpoc/taxonomy_bac.txt | awk '{print $1}') > actinomyces.fa
#trime out primers from amplicon
cutadapt -g MAYGARAARMGNATGYTNCARGA -o front_trim.fasta amplicons.fasta
cutadapt -l 478 -o trimmed_amplicon.fasta amplicons.fasta

cat actinomyces.fa trimmed_amplicon.fasta > tmp
mv tmp actinomyces.fa
mafft --thread 7 actinomyces.fa > actinomyces.align.fa

raxmlHPC-PTHREADS-SSE3 -T 7 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n actino.tre -s actinomyces.align.fa

grep Actinomyces -w ../../rpoc/taxonomy_bac.txt | awk '{print $1, $8}' | sed 's/ /\t/' > actino_annots.txt
awk '{print $1,$2}' results.txt | sed '1d' | sed 's/ /\t/' >> actino_annots.txt
```
