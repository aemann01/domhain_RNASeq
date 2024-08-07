library(phyloseq)
library(Biostrings)
setwd("~/rna_dohmain/rpoc/")
#rpoc
#frequency
# seqtab_rpoc <- read.table("~/rna_dohmain/rpoc/sequence_table.merged.txt", header=T, row.names=1)
# asv_rpoc <-otu_table(seqtab_rpoc, taxa_are_rows=F)
# #taxonomy
# tax_tab_rpoc <- read.table("~/rna_dohmain/rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tax_rpoc <- tax_table(as.matrix(tax_tab_rpoc))
# #ref sequences
# refseq_rpoc <- Biostrings::readDNAStringSet("~/rna_dohmain/rpoc/rep_set.fa")
# #metadata
# metadata <- read.table("~/rna_dohmain/rpoc/metadata.txt", sep="\t", header=T, row.names=1)
# metadata$age_y <- as.numeric(metadata$age_y)

# map <- sample_data(metadata)
# #merge into one phyloseq object
# rpoc.pd <- merge_phyloseq(asv_rpoc, tax_rpoc, refseq_rpoc, map)
# ps.dat <- prune_samples(sample_sums(rpoc.pd) > 4000, rpoc.pd)

# save.image("~/rna_dohmain/rpoc/ps.RData")

#start analysis
library("gridExtra")
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(corncob, verbose=F)
library(magrittr, verbose=F)
library(ggpubr, verbose=F)
library(ecole, verbose=F)
library(UpSetR)
library(smplot2)
library(RColorBrewer)
library(cowplot)

load("ps.RData")

#rarefaction curve
options(warn=-1) # suppress warnings
p <- ggrare(ps.dat, step = 1000, color = "hiv_status", se = TRUE)
p <- p + facet_wrap(~aliquot_type)
p + theme_minimal()
pdf("./rarefaction_plots.pdf")
p + theme_minimal()
dev.off()
options(warn=0) # back on
system("~/.iterm2/imgcat ./rarefaction_plots.pdf")
