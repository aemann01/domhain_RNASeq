# 1. Make phylogeny
```sh
cd ~/rna_dohmain/11-perio/10-p-ging-trees
grep hagA ../05-TrEMBL/combined.gff | awk -F "\t" '{print $1, $4, $5}' | sort | uniq > hagA.ids
seqtk subseq ../05-TrEMBL/combined.fna hagA.ids > hagA.fa
#rename headers to have gene id
grep hagA ../06-red-complex/red_annots.txt | awk '{print $1, $2}' | sed 's/ /_/' > combined.txt
grep '^>' hagA.fa | sed 's/:.*//' | sed 's/>//' | while read line; do grep $line combined.txt; done > headers.txt
paste -d "\t" <(grep '^>' hagA.fa) headers.txt | sed 's/>//g '> names.txt
seqkit replace -p "(.+)" -r '{kv}|$1' -k names.txt hagA.fa | sed 's/|.*//' > temp
mv temp hagA.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' hagA.fa > temp
mv temp hagA.fa
sed -i 's/.*_g/>g/' hagA.fa

#align
mafft --thread 7 hagA.fa > hagA.align.fa
sed -i 's/:.*//' hagA.align.fa
#tree
rm RAxML*
raxmlHPC-PTHREADS-SSE3 -T 7 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n hagA.tre -s hagA.align.fa

```
# 2. Run analysis
```R
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ape)
library(phytools)
library(ggtree)
library(gridExtra)
library(grid)
library(pheatmap)

setwd("~/rna_dohmain/11-perio/10-p-ging-trees")
HI <- read.csv("../06-red-complex/deseq_results_red-HIvHUU.txt",  sep = "\t")
# HEU <- read.csv("deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
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
#read in tree
tree <- read.tree("RAxML_bipartitions.hagA.tre")
md_tree <- midpoint_root(tree)
tip_order <- rev(md_tree$tip.label)
tip_order2 <- gsub("'","",tip_order)
#filter based on tip names
filt_HI <- HI[rownames(HI) %in% tip_order, ]

filt_HI_long <- filt_HI %>%
  pivot_longer(cols = c(baseMean, log2FoldChange), 
               names_to = "Type", 
               values_to = "Value")


tree_grob <- grid.grabExpr({
  plot(md_tree, main = "Phylogenetic Tree", show.tip.label = TRUE)
})

# Save the heatmap as a grob (graphical object)
log2 <- grid.grabExpr({ggplot(filt_HI, aes(x=factor(species),y=factor(tag, levels=tip_order2))) +
  geom_tile((aes(fill = log2FoldChange)))+
  scale_fill_gradient(low = "#FFC300", high = "#900C3F") +
  coord_fixed() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 2.5, vjust = 0.6, hjust = 0.0,), axis.text.y = element_text(size = 3),
        panel.background = element_rect(fill = 'grey67'), axis.title.x = element_text(size= 7), 
        axis.title.y = element_text(size= 6), plot.title = element_text(size = 10), 
        legend.title = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text = element_text(size = 4)) 
})

baseMean <- grid.grabExpr({ggplot(filt_HI, aes(x=factor(species),y=factor(tag, levels=tip_order2))) +
  geom_tile((aes(fill = baseMean)))+
  scale_fill_gradient(low = "#FFC300", high = "#900C3F") +
  coord_fixed() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 2.5, vjust = 0.6, hjust = 0.0,), axis.text.y = element_text(size = 3),
        panel.background = element_rect(fill = 'grey67'), axis.title.x = element_text(size= 7), 
        axis.title.y = element_text(size= 6), plot.title = element_text(size = 10), 
        legend.title = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text = element_text(size = 4)) 
})


# Arrange both the tree and heatmap side by side
pdf("hagA.heatmap.pdf")
grid.arrange(tree_grob, log2, ncol = 2)
dev.off()
system("~/.iterm2/imgcat ./hagA.heatmap.pdf")


#HEU