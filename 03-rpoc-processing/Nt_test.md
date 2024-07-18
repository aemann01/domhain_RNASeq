# 1. Assign taxonomy using nt database and kraken
```sh
#!/bin/bash

#SBATCH --job-name rpoc_nt_taxonomy
#SBATCH --nodes 1
#SBATCH --tasks-per-node 3
#SBATCH --cpus-per-task 1
#SBATCH --mem 950gb
#SBATCH --time 48:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

cd /scratch/scrull/hiv_rnaseq/rpoc
module add kraken2/2.1.2
kraken2 --db ../denovo-taxonomy/nt_kraken \
  --threads 7 \
  --use-names \
  --output rep_set.kraken.out rep_set.fa \
  --unclassified-out rep_set.unclassified.out --confidence 0.01
#   Loading database information... done.
# 8437 sequences (4.03 Mbp) processed in 3.698s (136.9 Kseq/m, 65.46 Mbp/m).
#   8418 sequences classified (99.77%)
#   19 sequences unclassified (0.23%)
```
Get full taxonomy
```sh
scp  scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/rpoc/rep_set.kraken.out ./
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ../database/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq ../rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ../fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
See how common A. naeslundii is for nt dabase using kraken2
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
setwd("~/rna_dohmain/rpoc/nt")
seqtab <- t(read.table("..//sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("./taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
#check prportion of A. naeslundii
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
# 2. Make local blastdb from homd
```sh
cat taxids | sed 's/$/|/' | while read line; do grep -m 1 ^$line ../rankedlineage_clean || echo $line "no lineage" ; done > lineages
parallel -a taxids -j 7 -k "grep -wm 1 ^'{}' ../rankedlineage_clean"> lineages
awk -F "|" '{print $1, $2}' lineages | sed 's/ /_/g' > lineas

paste rpoc_seqs lineas | sed 's/\t/|blast:taxid|/' | sed 's/\t/ /' | sed 's/^/>/' > blast_heads
paste blast_heads seqs | sed 's/\t/\n/' > rpoc_blast.fa

makeblastdb -in ../rpoc_blast.fa -dbtype nucl -out homd_blastdb
cd ~/rna_dohmain/rpoc/database/homd_blastdb
blastn -query ../../rep_set.fa -db ~/rna_dohmain/rpoc/database/homd_blastdb/homd_blastdb -out blastn_output.alignment -num_threads 7 -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
grep ASV blastn_output.alignment | awk '{print $1}' | sort | uniq | while read line; do grep -wm 1 $line blastn_output.alignment; done > best_match
tr -cd '\11\12\15\40-\176' < best_match > temp
grep -v suzanne temp > best_match
grep Actinomyces_naeslundii best_match -c
```
Get taxonomy
```sh
sed 's/.*taxid|//' best_match | sed 's/_.*//' | sed 's/ //'> taxids
sed "s/\t//g" ../rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $1}' best_match > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq ../../rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ./fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
See how common A. naeslundii is for homd dabase using blast
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
setwd("~/rna_dohmain/rpoc/database/homd_blastdb")
seqtab <- t(read.table("../../sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("../../../homd_map/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
#check prportion of A. naeslundii
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