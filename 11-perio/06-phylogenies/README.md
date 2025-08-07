# 1. Porphyromonas gingivalis tree
```sh
cd ~/bin
#setup iqtree3
wget https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux-intel.tar.gz 
tar -xzvf iqtree-3.0.1-Linux-intel.tar.gz
mv ~/bin/iqtree-3.0.1-Linux-intel/bin/iqtree ./
#setup 3seq
wget https://mol.ax/rxiv/3seq_build_170612.zip
unzip 3seq_build_170612.zip
mv 3seq\ build\ 170612 3seq_build_170612
cd 3seq_build_170612
make
mv 3seq ../
./3seq -g myPvalueTable500 500

# strart making phylogenies
cd ~/rna_dohmain/11-perio/06-phylogenies
mkdir p_gingivalis && cd p_gingivalis
grep Porphyromonas ../../02-pgap/SEQID_info.csv | grep gingivalis | grep -v cangingivalis | awk -F "," '{print $1}' | sed 's/\"//g' | awk '{print $1, "_modified_annot_cds_from_genomic.fna"}' | sed 's/ //' | while read line; do cp ~/rna_dohmain/11-perio/02-pgap/annotations/$line ./; done
# clean up headers
sed -i 's/ .*//' *_modified_annot_cds_from_genomic.fna
sed -i 's/lcl|//' *_modified_annot_cds_from_genomic.fna
cat *_modified_annot_cds_from_genomic.fna > all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 60 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
# ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
Rscript ../../pairwise.R
# 3seq -f c2502.align.fa -ptable ~/bin/myPvalueTable500 -id c2502.align.rec
# 3seq -f c1153.align.fa -id c1153.align.rec
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
cat recomb | while read line; do mv $line $line.rd; done
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 > 0.05 )' > no_rec
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 < 0.05 )' > recom

# ls c*align.fa | sed 's/.align.fa//' | while read line; do run_gubbins.py $line.align.fa --threads 130 --model GTRCAT --prefix $line.gubbins > $line.rec; done
# run_gubbins.py c3118.align.fa --threads 60 --model GTRCAT --prefix c3118.gubbins

# find . -type f -name "*.embl" -empty

# combine single copy core
# combine alignments
python3 ../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree2 -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core -safe

# make annotation file
cat rpoC.ids | sed 's/\..*//' | sed 's/lcl|//' | sort | uniq | while read line; do grep -m 1 $line ../../red_annots.txt | awk -F "\t" '{print $1, $6}' | sed 's/ /\t/' ; done | sed '1s/^/genome\tspecies\n/'> combined_core.annots.txt
```
# 2. Expression analysis for P. ginigvalis
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
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
load("../04-red-diff/deseq_results_red-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
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
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")

sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

sig_average_long <- sig_average %>%
  pivot_longer(cols = c(avg_baseMean, avg_log2FoldChange), 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# make graph
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

tr2 <- phylo4d(tree.root, res_ord_sub$species[match(tree.root$tip.label, res_ord_sub$genome)])
tr2 <- phylo4d(tree.root)

tree_plot <- ggtree(tr2)+ 
  geom_tiplab(align = TRUE, size =2.8,offset = .0033) +
  # geom_tippoint(aes(color=dt, x = .0085), alpha=1) +
  scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043")) +
  theme(legend.position = "none")+
  xlim(0, .015)

pdf("pging_heat.all_genes.pdf", width =10)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(3, 3))
dev.off()
system("~/.iterm2/imgcat ./pging_heat.all_genes.pdf")

# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_sub_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub$genome

tree.sub <- keep.tip(tree.root, sig_average_baseMean_sub$genome)
tr3 <- phylo4d(tree.sub, res_ord$species[match(tree.sub$tip.label, res_ord$genome)])
root_num <- getRoot(tree.sub)

tree_plot <- ggtree(tr3, ladderizes = FALSE) + 
  geom_tiplab(align = TRUE, size =2,offset = .0001) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .015)
tip_labels <- rev(get_taxa_name(tree_plot))

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_log2_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, guide = guide_colorbar(title = "Average log fold change for virulence factors", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, guide = guide_colorbar(title = "Average base mean for virulence factors", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)


pdf("pging_heat.viru_genes.pdf", width = 8, height = 8)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.75, .45))
dev.off()
system("~/.iterm2/imgcat ./pging_heat.viru_genes.pdf")
```
## 2.1 How common are most virulent strains
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
library(reshape2)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
#get relative abundance of rna
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
# select columns for a certain genome 
gene.sub <- genecounts[grep("SEQF3203", row.names(genecounts), invert = FALSE, ignore.case = TRUE),]
gene.df <- melt(gene.sub)
gene.annot <- left_join(gene.df, metadata, by = c("variable" = "sample_id"))
gene.annot $hiv_status <- factor(gene.annot$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("SEQF3203.RNA.sample.pdf")
ggplot() + geom_bar(data = gene.annot, aes(x = variable, y = value), position = "dodge", stat = "identity")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./SEQF3203.RNA.sample.pdf")

# select columns for a certain genome 
gene.sub <- genecounts[grep("SEQF3220", row.names(genecounts), invert = FALSE, ignore.case = TRUE),]
gene.df <- melt(gene.sub)
gene.annot <- left_join(gene.df, metadata, by = c("variable" = "sample_id"))
gene.annot $hiv_status <- factor(gene.annot$hiv_status, levels = c("HUU", "HEU", "HI"))

pdf("SEQF3220.RNA.sample.pdf")
ggplot() + geom_bar(data = gene.annot, aes(x = variable, y = value), position = "dodge", stat = "identity")+
  facet_grid(~ hiv_status +visit_num, switch = "x", scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./SEQF3220.RNA.sample.pdf")
```
## 2.2 Genome Cluster for P. ginigivalis
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
library(vegan)
library(ggdendro)
library(plotly)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
load("../04-red-diff/deseq_results_red-HIvHUU.RData")
# read in data
dat=t(read.table("./p_gingivalis/gene_cluster.tab", row.names = 1, header = TRUE))
d=vegdist(dat,method="jaccard")
hc=hclust(d,method="average")
hc
dhc <- as.dendrogram(hc)
data <- dendro_data(dhc, type = "rectangle")

# get heatmpa info# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\..*", replacement = "") 

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
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")

# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_sub_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")

ddata <- dendro_data(hc)
tip_labels <- ddata$labels$label
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_log2_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, guide = guide_colorbar(title = "Average log fold change for virulence factors", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, guide = guide_colorbar(title = "Average base mean for virulence factors", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  # scale_y_discrete(labels = ) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

tree_plot <- ggplot(segment(data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(data), aes(x = x, y = y, label = label), 
          hjust = 0, size = 3)  +
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))+
  theme_dendro()

pdf("gene_cluster.pging.pdf", width = 15)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(5, .75))
dev.off()
system("~/.iterm2/imgcat ./gene_cluster.pging.pdf")
```
# 3. Treponema denticola tree
```sh
cd ~/rna_dohmain/11-perio/06-phylogenies
mkdir t_denticola && cd t_denticola
grep Treponema ../../02-pgap/SEQID_info.csv | grep denticola | grep -v cangingivalis | awk -F "," '{print $1}' | sed 's/\"//g' | awk '{print $1, "_modified_annot_cds_from_genomic.fna"}' | sed 's/ //' | while read line; do cp ~/rna_dohmain/11-perio/02-pgap/annotations/$line ./; done
# clean up headers
sed -i 's/ .*//' *_modified_annot_cds_from_genomic.fna
sed -i 's/lcl|//' *_modified_annot_cds_from_genomic.fna
cat *_modified_annot_cds_from_genomic.fna > all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 60 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | while read line; do mafft --thread -1 $line > $line.align.fa; done
# combine alignments
sed -i 's/_.*//g' *align.fa
# test for recombination
Rscript ../../pairwise.R
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done # test for recombination
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
wc -l recomb
ls *align.fa | wc -l
cat recomb | while read line; do mv $line $line.rd; done
# combine non recombinant genes
python3 ../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tdent.core_genome -safe
iqtree2 -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tdent.core -safe
```
# 4. Expression analysis for T. denticola
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
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
load("../04-red-diff/deseq_results_red-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_denticola/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))

# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
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
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")

sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

sig_average_long <- sig_average %>%
  pivot_longer(cols = c(avg_baseMean, avg_log2FoldChange), 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# make graph
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

tr2 <- phylo4d(tree.root, res_ord_sub$species[match(tree.root$tip.label, res_ord_sub$genome)])
tr2 <- phylo4d(tree.root)

tree_plot <- ggtree(tr2)+ 
  geom_tiplab(align = TRUE, size =2.8,offset = .0033) +
  # geom_tippoint(aes(color=dt, x = .0085), alpha=1) +
  scale_color_manual(values = c("Treponema_denticola" = "#340043")) +
  theme(legend.position = "none")+
  xlim(0, .06)

pdf("tdent_heat.all_genes.pdf", width =10)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(3, 3))
dev.off()
system("~/.iterm2/imgcat ./tdent_heat.all_genes.pdf")

# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
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

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2_sub<- filter(sig_average_log2, species == "Treponema_denticola")
sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Treponema_denticola")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_sub_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub$genome

tree.sub <- keep.tip(tree.root, sig_average_baseMean_sub$genome)
tr3 <- phylo4d(tree.sub, res_ord$species[match(tree.sub$tip.label, res_ord$genome)])
root_num <- getRoot(tree.sub)

tr3 <- phylo4d(tree.sub, res_ord$species[match(tree.sub$tip.label, res_ord$genome)])
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
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_log2_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, guide = guide_colorbar(title = "Average log fold change for virulence factors", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, guide = guide_colorbar(title = "Average base mean for virulence factors", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)


pdf("tdent_heat.viru_genes.pdf", width = 12, height = 6)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(3, .75))
dev.off()
system("~/.iterm2/imgcat ./tdent_heat.viru_genes.pdf")
```
## 4.1 Genome Clustering
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
library(vegan)
library(ggdendro)
library(plotly)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
load("../04-red-diff/deseq_results_red-HIvHUU.RData")
# read in data
dat=t(read.table("./t_denticola/gene_cluster.tab", row.names = 1, header = TRUE))
d=vegdist(dat,method="jaccard")
hc=hclust(d,method="average")
hc
dhc <- as.dendrogram(hc)
data <- dendro_data(dhc, type = "rectangle")

# get heatmpa info# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\..*", replacement = "") 

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
res_ord_sub <- filter(res_ord, species == "Treponema_denticola")

# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
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

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2_sub<- filter(sig_average_log2, species == "Treponema_denticola")
sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Treponema_denticola")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_sub_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")

ddata <- dendro_data(hc)
tip_labels <- ddata$labels$label

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_log2_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, guide = guide_colorbar(title = "Average log fold change for virulence factors", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, guide = guide_colorbar(title = "Average base mean for virulence factors", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  # scale_y_discrete(labels = ) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

tree_plot <- ggplot(segment(data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(data), aes(x = x, y = y, label = label), 
          hjust = 0, size = 3)  +
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))+
  theme_dendro()

pdf("gene_cluster.tdent.pdf", width = 15)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(5, .75))
dev.off()
system("~/.iterm2/imgcat ./gene_cluster.tdent.pdf")
```
# 5. Tannerella forsythia tree
```sh
cd ~/rna_dohmain/11-perio/06-phylogenies
mkdir t_forsythia && cd t_forsythia
grep Tannerella ../../02-pgap/SEQID_info.csv | grep forsythia | grep -v cangingivalis | awk -F "," '{print $1}' | sed 's/\"//g' | awk '{print $1, "_modified_annot_cds_from_genomic.fna"}' | sed 's/ //' | while read line; do cp ~/rna_dohmain/11-perio/02-pgap/annotations/$line ./; done
# clean up headers
sed -i 's/ .*//' *_modified_annot_cds_from_genomic.fna
sed -i 's/lcl|//' *_modified_annot_cds_from_genomic.fna
cat *_modified_annot_cds_from_genomic.fna > all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 60 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | while read line; do mafft --thread -1 $line > $line.align.fa; done
# combine alignments
sed -i 's/_.*//g' *align.fa
# test for recombination
Rscript ../../pairwise.R
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done # test for recombination
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
wc -l recomb
ls *align.fa | wc -l
cat recomb | while read line; do mv $line $line.rd; done
# combine non recombinant genes
python3 ../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tforsythia.core_genome -safe
iqtree2 -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tforsythia.core -safe
```
# 6. Expression analysis for T. forsythia
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
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-phylogenies/")
load("../04-red-diff/deseq_results_red-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_forsythia/align/tforsythia.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))

# add in annotations
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
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
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")

sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

sig_average_long <- sig_average %>%
  pivot_longer(cols = c(avg_baseMean, avg_log2FoldChange), 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# make graph
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

tr2 <- phylo4d(tree.root, res_ord_sub$species[match(tree.root$tip.label, res_ord_sub$genome)])
tr2 <- phylo4d(tree.root)

tree_plot <- ggtree(tr2)+ 
  geom_tiplab(align = TRUE, size =2.8,offset = .0033) +
  # geom_tippoint(aes(color=dt, x = .0085), alpha=1) +
  scale_color_manual(values = c("Tanne" = "#340043")) +
  theme(legend.position = "none")+
  xlim(0, .06)

pdf("tforsythia_heat.all_genes.pdf", width =10)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(3, 3))
dev.off()
system("~/.iterm2/imgcat ./tforsythia_heat.all_genes.pdf")

# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )
# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Tannerella_forsythia")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(genome, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(genome, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2_sub<- filter(sig_average_log2, species == "Tannerella_forsythia")
sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Tannerella_forsythia")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_tf & species == "Tannerella_forsythia") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_sub_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_sub$genome

tree.sub <- keep.tip(tree.root, sig_average_baseMean_sub$genome)
tr3 <- phylo4d(tree.sub, res_ord$species[match(tree.sub$tip.label, res_ord$genome)])
root_num <- getRoot(tree.sub)

tr3 <- phylo4d(tree.sub, res_ord$species[match(tree.sub$tip.label, res_ord$genome)])
tree_plot <- ggtree(tr3) + 
  geom_tiplab(align = TRUE, size =3,offset = .0009) + 
  geom_rootedge(root_num) +
  # geom_tippoint(aes(color=dt, x = .012), alpha=1) +
  # scale_color_manual(values = c("Porphyromonas_gingivalis" = "#340043", 
                               # "Tannerella_forsythia" = "#FBE51F", 
                               # "Treponema_denticola" = "#1E7F7A")) +
  theme(legend.position = "none")+
  xlim(0, .009)
tip_labels <- rev(get_taxa_name(tree_plot))
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels=tip_labels))) +
  geom_tile(data = sig_average_log2_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_log2_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, guide = guide_colorbar(title = "Average log fold change for virulence factors", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean_sub_virus, aes(fill = value), color ="black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, guide = guide_colorbar(title = "Average base mean for virulence factors", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = NULL) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)


pdf("tforsythia_heat.viru_genes.pdf", width = 12, height = 6)
tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(3, .75))
dev.off()
system("~/.iterm2/imgcat ./tforsythia_heat.viru_genes.pdf")
```






