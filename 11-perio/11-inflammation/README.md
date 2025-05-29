# 1. Compare boxplots
```R
library(devtools)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ape)
library(microbiome)
library(ggpubr)
library(rstatix)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-inflammation")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")

hiv_stat <- c("HI", "HEU", "HUU")
pdf("oral_score.pdf")
ggplot(metadata, aes(x=factor(hiv_status, levels=hiv_stat),y=oral_hygiene_score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./oral_score.pdf")

pdf("gingival_score.pdf")
ggplot(metadata, aes(x=factor(hiv_status, levels=hiv_stat),y=gingival_inflammation_score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./gingival_score.pdf")

```
## 1.2 CD4 count relationship
```R
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
ords <- c("HI", "HEU", "HUU")
metadata$hiv_status <- factor(metadata$hiv_status, levels = ords)

pdf("cd4.oral.corr.pdf")
ggscatter(metadata, x = "oral_hygiene_score", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor(aes(color = hiv_status), method = "spearman")+
  geom_point(aes(color = hiv_status), shape = 21, size = 3, stroke = 1, color = "black")+
  scale_color_manual(values = hivCols) +
  scale_fill_manual(values = hivCols) +
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./cd4.oral.corr.pdf")
cor.test(metadata$gingival_inflammation_score, metadata$cd4_count, method = "pearson") #not sig

pdf("cd4.ging.corr.pdf")
ggscatter(metadata, x = "gingival_inflammation_score", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./cd4.ging.corr.pdf")

pdf("cd4_count.pdf")
ggplot(metadata, aes(x=hiv_status,y=cd4_count))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./cd4_count.pdf")
```
# 3. Relationship between oral score and virulence factors of perio
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
genecounts <- read.table("../06-red-complex/red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
#get annots for genes
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="")
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
combined_sub <- combined_data %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
# oral hygein score
oral_sum <- combined_sub %>% 
  group_by(variable, hiv_status, oral_hygiene_score) %>% 
  summarise(Total = sum(value, na.rm = TRUE))

cor.test(oral_sum$oral_hygiene_score, oral_sum$Total, method = "pearson") #sig
pdf("rna.viruvoral.corr.pdf")
ggscatter(oral_sum, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.viruvoral.corr.pdf")

# gingival inflammation score
ging_sum <- combined_sub %>% 
  group_by(variable, hiv_status, gingival_inflammation_score) %>% 
  summarise(Total = sum(value, na.rm = TRUE))
pdf("rna.viruvging.corr.pdf")
ggscatter(ging_sum, x = "gingival_inflammation_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.viruvging.corr.pdf")
cor.test(ging_sum$gingival_inflammation_score, ging_sum$Total, method = "pearson") #not sig
```
# 4. Overall RNA
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
merged_data <- merge(vld_long, homd, by.x = "Var1", by.y = "tag")
merged_df2 <- merge(merged_data, metadata, by.x = "Var2", by.y = "sample_id")

#sum all genes in a transcirpt together
red_mean <- merged_df2 %>% 
  group_by(Var2, contig, species, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = sum(value, na.rm = TRUE))
#average by genome
red_mean2 <- red_mean %>% 
  group_by(Var2, species, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total_spec = mean(Total, na.rm = TRUE))
#average by species
red_mean3 <- red_mean2 %>% 
  group_by(Var2, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = mean(Total_spec, na.rm = TRUE))
#oral 
pdf("rna.redvoral.corr.pdf")
ggscatter(red_mean3, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.redvoral.corr.pdf")

hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
ords <- c("HI", "HEU", "HUU")
red_mean3$hiv_status <- factor(red_mean3$hiv_status, levels = ords)
pdf("rna.redvoral.corr.pdf", height =7, width =15)
ggscatter(red_mean3, x = "oral_hygiene_score", y = "Total", size = 3,
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
  stat_cor(aes(color = hiv_status), label.x = 3.5, size=6)+
  geom_point(aes(color = hiv_status), shape = 21, size = 3, stroke = 1, color = "black")+
  scale_color_manual(values = hivCols) +
  scale_fill_manual(values = hivCols) +
  geom_xsidedensity(
    aes(
      y    = ..count../sum(..count..),
      xfill = factor(hiv_status, levels =ords),
    ),
    alpha    = 0.5,
    size     = 1.5,
    position = "stack"
  ) +
  scale_xsidey_continuous(minor_breaks = NULL)+
  scale_xfill_manual(values = hivCols)+
  geom_ysidedensity(
    aes(
      x    = ..count../sum(..count..),
      yfill = factor(hiv_status, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1.5,
   position = "stack"
  )+
  theme(axis.text.x = element_text(angle = 90, vjust = .5))+
  scale_yfill_manual(values = hivCols)+
  theme_bw()+
  theme(
    legend.position = "none",  
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text.x = element_text(size = 16),   
    axis.text.y = element_text(size = 16),
    ggside.axis.text.y= element_text(size = 10),
    ggside.axis.text.x= element_text(size = 10)
  ) +
  scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL)+
  labs(
    x = "Oral Hygiene Score",
    y = "Averaged Transcirptome Expression" 
  )
dev.off()
system("~/.iterm2/imgcat ./rna.redvoral.corr.pdf")

# gingivals inflatmattion
pdf("rna.redvging.corr.pdf")
ggscatter(red_mean3, x = "gingival_inflammation_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.redvging.corr.pdf")

# cd4
pdf("rna.redvcd4.corr.pdf")
ggscatter(red_mean3, x = "Total", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor(aes(color = hiv_status), method = )+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.redvcd4.corr.pdf")
```
# 5. rpoC RNA
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
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/11-inflammation")
#RNA
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
proc.time() - ptm 
#pcoa diversity
vld <- varianceStabilizingTransformation(se_star)
#pull out just rpoC
homd <- read.table("../06-red-complex/red_annots.txt", header=T, sep="\t", quote="") 
rownames(homd) <- homd$tag
rpoC_tags <- homd$tag[homd$gene == "rpoC"]
vld_rna <- melt(assay(vld)[rpoC_tags,])
# filter by locus tag 
merged_data <- merge(vld_rna, homd, by.x = "Var1", by.y = "tag")
merged_df2 <- merge(merged_data, metadata, by.x = "Var2", by.y = "sample_id")

red_mean <- merged_df2 %>% 
  group_by(Var2, species, hiv_status, oral_hygiene_score, gingival_inflammation_score) %>% 
  summarise(Total = mean(value, na.rm = TRUE))
red_mean2 <- red_mean %>% 
  group_by(Var2, hiv_status, oral_hygiene_score, gingival_inflammation_score) %>% 
  summarise(Total = mean(Total, na.rm = TRUE))
#oral 
pdf("rna_rpoc.redvoral.corr.pdf")
ggscatter(red_mean2, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna_rpoc.redvoral.corr.pdf")

hivCols <- c("#40A0FA", "#FA78FA", "#8213A0")
ords <- c("HI", "HEU", "HUU")

pdf("rna_rpoc.redvoral.corr.pdf", height =7, width =15)
ggscatter(red_mean2, x = "oral_hygiene_score", y = "Total", size = 3,
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
  stat_cor(aes(color = hiv_status), label.x = 3.5, size=6)+
  geom_point(aes(color = hiv_status), shape = 21, size = 3, stroke = 1, color = "black")+
  scale_color_manual(values = hivCols) +
  scale_fill_manual(values = hivCols) +
  geom_xsidedensity(
    aes(
      y    = ..count../sum(..count..),
      xfill = factor(hiv_status, levels =ords),
    ),
    alpha    = 0.5,
    size     = 1.5,
    position = "stack"
  ) +
  scale_xsidey_continuous(minor_breaks = NULL)+
  scale_xfill_manual(values = hivCols)+
  geom_ysidedensity(
    aes(
      x    = ..count../sum(..count..),
      yfill = factor(hiv_status, levels =ords)
    ),
    alpha    = 0.5,
    size     = 1.5,
   position = "stack"
  )+
  theme(axis.text.x = element_text(angle = 90, vjust = .5))+
  scale_yfill_manual(values = hivCols)+
  theme_bw()+
  theme(
    legend.position = "none",  
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text.x = element_text(size = 16),   
    axis.text.y = element_text(size = 16),
    ggside.axis.text.y= element_text(size = 10),
    ggside.axis.text.x= element_text(size = 10)
  ) +
  scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL)+
  labs(
    x = "Oral Hygiene Score",
    y = "Averaged RNA rpoC Expression" 
  )
dev.off()
system("~/.iterm2/imgcat ./rna_rpoc.redvoral.corr.pdf")

# gingivals inflatmattion
pdf("rna_rpoc.redvging.corr.pdf")
ggscatter(red_mean, x = "gingival_inflammation_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna_rpoc.redvging.corr.pdf")
```
# 5. Normalize virulence RNA
```R
homd_sub <- homd %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA" | gene == "tfsA" | gene == "tfsB")
viru_genes <- homd_sub$tag
virus_vld <- assay(vld)[viru_genes,]
vld_long <- melt(virus_vld)
merged_data <- merge(vld_long, homd_sub, by.x = "Var1", by.y = "tag")
merged_df2 <- merge(merged_data, metadata, by.x = "Var2", by.y = "sample_id")

viru_mean <- merged_df2 %>% 
  group_by(Var2, hiv_status, oral_hygiene_score, gingival_inflammation_score) %>% 
  summarise(Total = mean(value, na.rm = TRUE))

#oral 
pdf("rna.viruvoral.corr.pdf")
ggscatter(viru_mean, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.viruvoral.corr.pdf")
# inflammation
pdf("rna.viruvging.corr.pdf")
ggscatter(viru_mean, x = "gingival_inflammation_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.viruvging.corr.pdf")
```
# 6. Look at impact with DNA
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("/home/suzanne/rna_dohmain/11-perio/11-inflammation")
load("../../rpoc/ps.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(ps.dat, "clr")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, Abundance, V8)
red_dna$nucl <- "dna"
red_dna <- red_dna %>%
  dplyr::rename(sample = Sample)
red_dna <- red_dna %>%
  dplyr::rename(species = V8)
red_dna <- red_dna %>%
  dplyr::rename(value = Abundance)

dna <- red_dna %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score) %>% 
  summarise(Total = sum(value, na.rm = TRUE), .groups = "drop")
#oral 
pdf("dna.oral.corr.pdf")
ggscatter(dna, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./dna.oral.corr.pdf")
# inflammation
pdf("dna.ging.corr.pdf")
ggscatter(dna, x = "gingival_inflammation_score", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./dna.ging.corr.pdf")


hiv_stat <- c("HI", "HEU", "HUU")
pdf("dna.red.pdf")
ggplot(dna, aes(x=factor(hiv_status, levels=hiv_stat),y=Total))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./dna.red.pdf")
