# 1. Compare oral, gingival, and CD4 for 93 samples
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
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
ords <- c("HI", "HEU", "HUU")
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
metadata$visit_num <- factor(metadata$visit_num)
pdf("oral_score.visit.pdf")
ggplot(metadata, aes(x=visit_num,y=oral_hygiene_score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # facet_wrap(~hiv_status)
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./oral_score.visit.pdf")

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


pdf("calculus_score.pdf")
ggplot(metadata, aes(x=factor(hiv_status, levels=hiv_stat),y=calculus_index))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./calculus_score.pdf")

# cd4 relationship
pdf("cd4.oral.corr.pdf")
ggscatter(metadata, x = "oral_hygiene_score", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./cd4.oral.corr.pdf")
cor.test(metadata$cd4_count, metadata$oral_hygiene_score, method = "pearson") # sig

pdf("cd4.ging.corr.pdf")
ggscatter(metadata, x = "gingival_inflammation_score", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./cd4.ging.corr.pdf")
cor.test(metadata$cd4_count, metadata$gingival_inflammation_score, method = "pearson") #not sig

pdf("cd4_count.pdf")
ggplot(metadata, aes(x=factor(hiv_status, levels=hiv_stat),y=cd4_count))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./cd4_count.pdf")
map$hiv_status <- factor(map$hiv_status, levels=c("HI", "HEU", "HUU"))
pdf("cd4_count.pdf")
ggplot(metadata, aes(x=hiv_status, y=log10(as.integer(cd4_count)))) + 
    geom_point(aes(col=log10(as.integer(cd4_count))), position="jitter") + 
    scale_color_viridis() + 
    theme_bw() +
    geom_hline(yintercept=c(log10(500), log10(1500)), linetype="dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    ylab("log10(CD4 count)") +
    xlab("HIV status") +
    theme(legend.title=element_blank())
dev.off()
system("~/.iterm2/imgcat ./cd4_count.pdf")


```
# 2. Compare oral, gingival, and CD4 for ~2000 samples
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
library(viridis)
library(rstatix)
library(effsize)
librry(lsr)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T)
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
meta$hiv_status <- factor(meta$hiv_status, levels = hiv_stat)

wilcox_test(meta, Oral_Hygiene_Score ~ sex) %>%
  add_significance()
meta$Oral_Hygiene_Score <- as.numeric(meta$Oral_Hygiene_Score)
meta$total_Ca_mg <- as.numeric(meta$total_Ca_mg)
meta$pH <- as.numeric(meta$pH)
meta$flow_rate <- as.numeric(meta$flow_rate)
meta$cd4_count <- as.numeric(meta$cd4_count)
meta$Calculus_index <- as.numeric(meta$Calculus_index)
meta$Gingival_inflammation_Score <- as.numeric(meta$Gingival_inflammation_Score)

cor.test(meta$pH, meta$Oral_Hygiene_Score, method = "spearman")
cor.test(meta$flow_rate, meta$Oral_Hygiene_Score, method = "spearman")
cor.test(meta$total_Ca_mg, meta$Oral_Hygiene_Score, method = "pearson")
cor.test(meta$total_vitA_RAE_mcg, meta$Oral_Hygiene_Score, method = "pearson")
cor.test(meta$total_vitC_mg, meta$Oral_Hygiene_Score, method = "pearson")
cor.test(meta$total_carbohydrates_g, meta$Oral_Hygiene_Score, method = "spearman")
cor.test(meta$Calculus_index, meta$Gingival_inflammation_Score, method = "pearson")

cor.test(meta$cd4_count, meta$Oral_Hygiene_Score, method = "spearman")

wilcox_test(meta, Oral_Hygiene_Score ~ hiv_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()    
aggregate(Oral_Hygiene_Score ~ hiv_status, data = meta, FUN = mean, na.rm = TRUE)
submap <- meta[meta$hiv_status == "HEU",]

cohen.d(Oral_Hygiene_Score ~ hiv_status, data = submap)
submap <- meta[meta$hiv_status == "HUU" | meta$hiv_status == "HI",]
cohen.d(Oral_Hygiene_Score ~ hiv_status, data = submap)

wilcox_test(meta, total_Ca_mg ~ hiv_status) %>%
  add_significance()
aggregate(total_Ca_mg ~ hiv_status, data = meta, FUN = mean, na.rm = TRUE)

wilcox_test(meta, total_vitA_RAE_mcg ~ hiv_status) %>%
  add_significance()
aggregate(total_vitA_RAE_mcg ~ hiv_status, data = meta, FUN = mean, na.rm = TRUE)

pdf("oral_score.long.pdf")
ggplot(meta, aes(x=factor(hiv_status, levels=hiv_stat),y=Oral_Hygiene_Score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./oral_score.long.pdf")

meta$visit_num <- factor(meta$visit_num)
pdf("oral_score.visit.long.pdf")
ggplot(meta, aes(x=visit_num,y=Oral_Hygiene_Score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  facet_wrap(~hiv_status)
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./oral_score.visit.long.pdf")

pdf("gingival_score.long.pdf")
ggplot(meta, aes(x=factor(hiv_status, levels=hiv_stat),y=Gingival_inflammation_Score))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./gingival_score.long.pdf")

pdf("calculus_score.long.pdf")
ggplot(meta, aes(x=factor(hiv_status, levels=hiv_stat),y=Calculus_index))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./calculus_score.long.pdf")

# calcium
meta$Ca_category <- ifelse(meta$total_Ca_mg <= 440, "low", ifelse(meta$total_Ca_mg >= 600, "high", "normal"))
meta$sample_visit <- paste0(meta$study_id,"V",meta$visit_num)
meta_sub <- meta[!is.na(meta$total_Ca_mg), ]
meta_unique <- meta_sub[!duplicated(meta_sub$sample_visit), ]
meta_unique <- meta_unique[meta_unique$total_Ca_mg < 2000, ]

meta$Oral_Hygiene_Score_Remark <- factor(meta$Oral_Hygiene_Score_Remark, levels=c("Good", "Fair", "Poor"))

pdf("total_Ca_mg.long.pdf")
ggplot(meta_unique, aes(x=factor(hiv_status, levels=hiv_stat),y=total_Ca_mg))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./total_Ca_mg.long.pdf")

# correlation

pdf("cd4voral.corr.long.pdf")
ggscatter(meta, x = "Oral_Hygiene_Score", y = "cd4_count",
   color = "hiv_status", shape = 21, size = 3, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = ",", size=5)
   )+
  scale_color_manual(values = hivCols)   
dev.off()
system("~/.iterm2/imgcat ./cd4voral.corr.long.pdf")

meta <- meta[meta$total_Ca_mg < 2000, ]
meta_unique$hiv_status <- factor(meta_unique$hiv_status, levels = hiv_stat)

pdf("CAvoral.corr.long.pdf")
ggscatter(meta, x = "Oral_Hygiene_Score", y = "total_Ca_mg",
   color = "hiv_status", shape = 21, size = 3, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = ",", size=5)
   )+
  scale_color_manual(values = hivCols)   
dev.off()
system("~/.iterm2/imgcat ./CAvoral.corr.long.pdf")

pdf("gingivalVcalculus.corr.long.pdf",)
ggscatter(meta, x = "Calculus_index", y = "Gingival_inflammation_Score",
   color = "hiv_status", shape = 21, size = 3, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.x = 1.5, label.sep = ",", size=5),
   point = "FALSE"
   )+
  geom_jitter(aes(color = hiv_status), width = 0.1, height = 0.1, size = 2, shape = 21, stroke =1) +  # Add jitter here
  scale_color_manual(values = hivCols)   
dev.off()
system("~/.iterm2/imgcat ./gingivalVcalculus.corr.long.pdf")


cor.test(meta$cd4_count, meta$Oral_Hygiene_Score, method = "spearman") # sig
cor.test(meta$total_Ca_mg, meta$Oral_Hygiene_Score, method = "spearman") # sig


pdf("cd4.ging.corr.long.pdf")
ggscatter(meta, x = "Gingival_inflammation_Score", y = "cd4_count",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./cd4.ging.corr.long.pdf")
cor.test(meta$cd4_count, meta$Gingival_inflammation_Score, method = "pearson") # sig

pdf("cd4_count.long.pdf")
ggplot(meta, aes(x=factor(hiv_status, levels=hiv_stat),y=cd4_count))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  # geom_jitter(aes(color=month), shape=16, position=position_jitter(0.2), size=2.5)+
  # scale_color_manual(values = month_colors)+ #color dots by sample
  # labs(x ="Bimonthly", y = "CLR Abundance")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./cd4_count.long.pdf")


meta$hiv_status <- factor(meta$hiv_status, levels=c("HI", "HEU", "HUU"))
meta$Oral_Hygiene_Score_Remark <- factor(meta$Oral_Hygiene_Score_Remark, levels=c("Good", "Fair", "Poor"))
pdf("oral_score.long.pdf")
ggplot(meta, aes(x = hiv_status, y = as.numeric(Oral_Hygiene_Score))) + 
    geom_point(aes(color = Oral_Hygiene_Score_Remark), position =  position_jitter(height = 0.05)) +
    scale_color_manual(values = c("Good" = "#006164", "Fair" = "#EDA247", "Poor" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns = TRUE, p.adjust.method = "fdr") + # Add pairwise comparisons
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    ylab("Simplified Oral Hygiene Score") +
    xlab("HIV status") +
    theme(legend.title = element_blank())
dev.off()
system("~/.iterm2/imgcat ./oral_score.long.pdf")


pdf("oral_score.violin.long.pdf")
ggplot(meta, aes(x=factor(hiv_status, levels=hiv_stat), y=Oral_Hygiene_Score, fill=hiv_status)) +
  geom_violin(trim=FALSE) +  # Create the violin plot, 'trim=FALSE' ensures the plot shows the full range of data
  geom_boxplot(width=0.1, color="black", alpha=0.5) + # Overlay a boxplot on top for additional details
  scale_fill_manual(values = hivCols) +  # Apply custom colors
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns = FALSE, p.adjust.method = "fdr") + # Add pairwise comparisons
  labs(x = "HIV Status", y = "Oral Hygiene Score") +
  theme_classic() +
  theme(legend.position="none") # Remove legend
dev.off()
system("~/.iterm2/imgcat ./oral_score.violin.long.pdf")

# Display the plot (use your terminal for image preview if using macOS)
system("~/.iterm2/imgcat ./oral_score_long_violin_plot.pdf")

meta$tooth_health <- factor(meta$tooth_health, levels=c("H", "E", "D"))
pdf("oralvtooth.long.pdf")
ggplot(meta, aes(x=tooth_health, y=as.integer(Oral_Hygiene_Score))) + 
    geom_point(aes(col=as.integer(Oral_Hygiene_Score)), position="jitter") + 
    scale_color_viridis() + 
    theme_bw() +
    geom_hline(yintercept=c(1.2, 3.1), linetype="dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    ylab("Oral Hygiene Score") +
    xlab("HIV status") +
    theme(legend.title=element_blank())
dev.off()
system("~/.iterm2/imgcat ./oralvtooth.long.pdf")
```
# 3. Overall RNA impact on oral, gingival, and CD4
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
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/red_counts.txt", header=T, sep="\t", row.names=1)
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
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
merged_data <- left_join(vld_long, homd, by = c("Var1" = "tag"))
merged_df2 <- left_join(merged_data, metadata, by = c("Var2" = "sample_id"))

#sum all genes in a transcirpt together
red_mean <- merged_df2 %>% 
  group_by(Var2, genome, species, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
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
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.redvoral.corr.pdf")
cor.test(red_mean3$oral_hygiene_score, red_mean3$Total, method = "pearson") # sig

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


hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
ords <- c("HI", "HEU", "HUU")
red_mean3$hiv_status <- factor(red_mean3$hiv_status, levels = ords)
pdf("rna.redvoral.corr.pdf", height =7, width =15)
ggscatter(red_mean3, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status", shape = 21, size = 3, stroke =2, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = ",", size=10)
   )+
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
    axis.title.x = element_text(size = 22), 
    axis.title.y = element_text(size = 22), 
    axis.text.x = element_text(size = 18),   
    axis.text.y = element_text(size = 18),
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
cor.test(red_mean3$Total, red_mean3$gingival_inflammation_score, method = "pearson") # not sig

# cd4
pdf("rna.redvcd4.corr.pdf")
ggscatter(red_mean3, x = "cd4_count", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.redvcd4.corr.pdf")
cor.test(red_mean3$Total, red_mean3$cd4_count, method = "pearson") # not sig
```
## 3.1. Normalize virulence RNA
```R
homd_sub <- homd %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA" | gene == "tfsA" | gene == "tfsB")
viru_genes <- homd_sub$tag
virus_vld <- assay(vld)[viru_genes,]
vld_long <- melt(virus_vld)
merged_data <- left_join(vld_long, homd_sub, by = c("Var1" = "tag"))
merged_df2 <- left_join(merged_data, metadata, by = c("Var2" = "sample_id"))
viru_mean <- merged_df2 %>% 
  group_by(Var2, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = median(value, na.rm = TRUE))

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
# cd4
pdf("rna.viruvcd4.corr.pdf")
ggscatter(viru_mean, x = "cd4_count", y = "Total",
   color = "hiv_status",
   add = "reg.line", conf.int = TRUE, cor.method="pearson")+
 stat_cor(aes(color = hiv_status))+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rna.viruvcd4.corr.pdf")
```
## 3.2 Composition of sample that is transcripts belong to red complex
```R
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggpubr)
library(ggside)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
counts <- read.csv("../03-global-diff/species_reads.txt", header=T,sep = "\t")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
counts$sample <- gsub(x = counts$sample, pattern = "\\.red", replacement = "") 
counts$sample <- gsub(x = counts$sample, pattern = "\\.", replacement = "-") 
data <- melt(counts)
data <- left_join(data, metadata, by = join_by(sample ==  sample_id))
data$genus <- gsub(x = data$variable, pattern = "_.*", replacement = "")

sub_data <- data[data$variable == "Porphyromonas_gingivalis" | data$variable == "Treponema_denticola" | data$variable == "Tannerella_forsythia",]
sub_data %>%
  group_by(variable) %>%
  summarise(
    cor_value = cor(value, oral_hygiene_score, method = "pearson"),  # Compute correlation directly
    p_value = cor.test(value, oral_hygiene_score, method = "pearson")$p.value,  # Extract p-value
    .groups = 'drop'
  )
cor.test(sub_data$value/100, sub_data$oral_hygiene_score, method = "pearson") #not sig
ords <- c("HI", "HEU", "HUU")
sub_data$hiv_status <- factor(sub_data$hiv_status, levels = ords)
sub_data$abundance <- sub_data$value/100
sub_data$variable <- factor(sub_data$variable, levels = c("Porphyromonas_gingivalis", "Treponema_denticola", "Tannerella_forsythia"))

pdf("redvoral.corr.species.pdf", width =20)
ggscatter(sub_data, x = "oral_hygiene_score", y = "abundance",
   color = "hiv_status", shape = 21, size = 3, stroke =2, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = ",", size=5)
   )+
  scale_color_manual(values = hivCols)   +
  facet_wrap(~variable, scales = "free") +
  labs(
    x = "Simplified Oral Hygiene Score",
    y = "Relative Abundance of Transcripts" 
  )
dev.off()
system("~/.iterm2/imgcat ./redvoral.corr.species.pdf")

# correlate between species
p_data <- sub_data[sub_data$variable == "Porphyromonas_gingivalis",]
p_data <- p_data %>%
  dplyr::rename(p_ging = abundance)
td_data <- sub_data[sub_data$variable == "Treponema_denticola",]
td_data <- td_data %>%
  dplyr::rename(t_dent = abundance)
tf_data <- sub_data[sub_data$variable == "Tannerella_forsythia",]
tf_data <- tf_data %>%
  dplyr::rename(t_for = abundance)
comb<-rbind(p_data, td_data, tf_data)
cor.test(red_mean3$p_ging, red_mean3$t_dent, method = "pearson") #not sig

#average by species
red_mean3 <- sub_data %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count, calculus_index) %>% 
  summarise(Total = sum(log10(value/100+0.0001), na.rm = TRUE))
cor.test(red_mean3$Total, red_mean3$oral_hygiene_score, method = "pearson") #not sig

test <- sub_data %>% 
  group_by(sample, hiv_status, oral_hygiene_score, oral_hygiene_score_remark, gingival_inflammation_score, cd4_count, calculus_index) %>% 
  summarise(Total = sum(value), na.rm = TRUE)
valid_data <- test[test$Total > 0.00, ]
cor.test(test$oral_hygiene_score, test$Total, method = "pearson")

test$log10_count <- log10(test$Total)
wilcox_test(as.data.frame(test), log10_count ~ hiv_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()    
aggregate(log10_count ~ hiv_status, data = test, FUN = mean, na.rm = TRUE)
# cor.test(1 / (red_mean3$Total + 1), red_mean3$oral_hygiene_score, method = "pearson")
# cor.test(red_mean3$Total, scale(red_mean3$oral_hygiene_score, center = TRUE, scale = TRUE), method = "sp")

# red_mean3$rounded <- round(as.numeric(red_mean3$Total), 2)

cor.test(red_mean3$Total, red_mean3$oral_hygiene_score, method = "pearson") #not sig
red_mean3$prop <- red_mean3$Total*100
# make corr graph
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
ords <- c("HI", "HEU", "HUU")
test$hiv_status <- factor(test$hiv_status, levels = ords)
# red_mean3$Total <- as.integer(red_mean3$Total)



pdf("rna.redvoral.relcorr.pdf", height =7, width =15)
ggscatter(red_mean3, x = "oral_hygiene_score", y = "Total",
   color = "hiv_status", shape = 21, size = 3, stroke =2, 
   add = "reg.line",  
   add.params = list(color = "black", fill = "darkgray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = ",", size=10)
   )+
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
    x = "Simplified Oral Hygiene Score",
    y = "Relative Abundance of Red Complex Transcripts" 
  )
dev.off()
system("~/.iterm2/imgcat ./rna.redvoral.relcorr.pdf")

```
## 3.3 Looking at composition in R
```R
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
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation/")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
ann <- homd[homd$tag %in% rownames(genecounts),]
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(genecounts), rownames(ann)))]
genecounts <- genecounts[sortrow, , drop=FALSE]
genecounts <- genecounts + 1
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(genecounts)==rownames(ann)) # should all return true
# if all are true, merge together
genecounts <- cbind(genecounts, ann)

# collapse by species and sum across rows
merge.count <- genecounts %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- genecounts %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
# find total transcript number
total_sums <- group_sums %>%
  summarize(across(where(is.numeric), sum))
# find the precentage 
group_percentages <- group_sums %>%
  mutate(across(where(is.numeric), ~ .x / total_sums[[cur_column()]] * 100))
# long format
collapsed_long <- group_percentages %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
red_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis" | homd$species  == "Treponema_denticola" | homd$species  == "Tannerella_forsythia",]$genome))

collapsed_long2 <- collapsed_long[collapsed_long$genome %in% red_seqs,]
sub_data <- left_join(collapsed_long2, metadata, by = c("sample" = "sample_id"))


red_mean3 <- sub_data %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = sum(count, na.rm = TRUE))
cor.test(log10(red_mean3$Total), red_mean3$oral_hygiene_score, method = "pearson") #not sig
results <- red_mean3 %>%
  group_by(hiv_status) %>%
  summarise(
    cor_value = cor(log10(Total), oral_hygiene_score, method = "pearson"),  # Compute correlation directly
    p_value = cor.test(log(Total), oral_hygiene_score, method = "pearson")$p.value,  # Extract p-value
    .groups = 'drop'
  )

# Print the results
print(results)

# look at using log 10
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
red_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis" | homd$species  == "Treponema_denticola" | homd$species  == "Tannerella_forsythia",]$genome))
collapsed_long2 <- collapsed_long[collapsed_long$genome %in% red_seqs,]
sub_data <- left_join(collapsed_long2, metadata, by = c("sample" = "sample_id"))
red_mean3 <- sub_data %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = sum(log10_count, na.rm = TRUE))
cor.test(red_mean3$Total, red_mean3$oral_hygiene_score, method = "spearman") #not sig

summary_log10_counts <- rna_abund %>%
  group_by(hiv_status, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))

# normalize by genome
# Assuming your data is in a dataframe called 'genecounts'
# Load necessary library
# Assuming 'genome' column contains species identifiers (e.g., SEQF4031.1 for the species name)
genome_counts <- genecounts %>%
  group_by(species) %>%
  summarise(genome_count = n_distinct(genome))

# Step 2: Count the number of occurrences of each species in the samples
species_counts <- genecounts %>%
  gather(key = "sample", value = "count", -gene, -tag, -product, -species) %>%
  filter(count > 0) %>%
  group_by(species) %>%
  summarise(species_count = n())

# Step 3: Merge the two counts to calculate the proportion
proportions <- merge(species_counts, genome_counts, by = "species") %>%
  mutate(proportion = species_count / genome_count)

# View the results
head(proportions)

long_data <- genecounts %>%
  gather(key = "sample", value = "count", -gene, -tag, -product, -species)

# Step 2: Filter out zero counts (species present in a sample)
long_data <- long_data %>%
  filter(count > 0)

# Step 3: For each sample, count how many species are present
species_per_sample <- long_data %>%
  group_by(sample) %>%
  summarise(total_species_in_sample = n_distinct(species))

# Step 4: For each sample and species, calculate the proportion
proportions_per_sample <- long_data %>%
  group_by(sample, species) %>%
  summarise(species_count_in_sample = n()) %>%
  left_join(species_per_sample, by = "sample") %>%
  mutate(proportion = species_count_in_sample / total_species_in_sample)

# View the results
head(proportions_per_sample)


```
## 3.4 medain for each species
```R
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
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation/")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
ann <- homd[homd$tag %in% rownames(genecounts),]
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(genecounts), rownames(ann)))]
genecounts <- genecounts[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(genecounts)==rownames(ann)) # should all return true
# if all are true, merge together
genecounts <- cbind(genecounts, ann)

# find the median for each genome
group_sums <- genecounts %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
merged_data <- as.data.frame(left_join(group_sums, unique(homd %>% select(genome, species)), by = "genome"))
species_medians <- merged_data %>%
  group_by(species) %>%
  summarise(across(starts_with("DM"), ~ median(.x, na.rm = TRUE)))

# now find the proportion of community
group_sspecies <- species_medians %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
# find total transcript number
total_sums <- group_sspecies %>%
  summarize(across(where(is.numeric), sum))
# find the precentage 
group_percentages <- group_sspecies %>%
  mutate(across(where(is.numeric), ~ .x / total_sums[[cur_column()]] * 100))
# long format
collapsed_long <- group_percentages %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
red_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis" | homd$species  == "Treponema_denticola" | homd$species  == "Tannerella_forsythia",]$genome))


collapsed_long2 <- collapsed_long[collapsed_long$species == "Porphyromonas_gingivalis" | collapsed_long$species  == "Treponema_denticola" | collapsed_long$species  == "Tannerella_forsythia",]
sub_data <- left_join(collapsed_long2, metadata, by = c("sample" = "sample_id"))
sub_data %>%
  group_by(species) %>%
  summarise(
    cor_value = cor(count, oral_hygiene_score, method = "pearson"),  # Compute correlation directly
    p_value = cor.test(count, oral_hygiene_score, method = "pearson")$p.value,  # Extract p-value
    .groups = 'drop'
  )

red_mean3 <- sub_data %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, cd4_count) %>% 
  summarise(Total = sum(count, na.rm = TRUE))
cor.test(red_mean3$Total, red_mean3$oral_hygiene_score, method = "pearson") #not sig
results <- red_mean3 %>%
  group_by(hiv_status) %>%
  summarise(
    cor_value = cor(Total, oral_hygiene_score, method = "pearson"),  # Compute correlation directly
    p_value = cor.test(Total, oral_hygiene_score, method = "pearson")$p.value,  # Extract p-value
    .groups = 'drop'
  )

# Print the results
print(results)
```
# 4. Look at impact with DNA
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("../../rpoc/ps.RData")
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
glom <- microbiome::transform(glom, method = "add", value = 1)
rel <- microbiome::transform(glom, "compositional")
actino <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
# glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(actino) # create dataframe from phyloseq object
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
red_dna <- select(data, Sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, total_Ca_mg, Abundance, V8)
red_dna$nucl <- "dna"
red_dna <- red_dna %>%
  dplyr::rename(sample = Sample)
red_dna <- red_dna %>%
  dplyr::rename(species = V8)
red_dna <- red_dna %>%
  dplyr::rename(value = Abundance)
red_dna %>%
  group_by(species) %>%
  summarise(
    cor_value = cor(value, oral_hygiene_score, method = "pearson"),  # Compute correlation directly
    p_value = cor.test(value, oral_hygiene_score, method = "pearson")$p.value,  # Extract p-value
    .groups = 'drop'
  )
dna <- red_dna %>% 
  group_by(sample, hiv_status, oral_hygiene_score, gingival_inflammation_score, total_Ca_mg) %>% 
  summarise(Total = sum(value, na.rm = TRUE), .groups = "drop")
#oral 
cor.test(log10(dna$oral_hygiene_score), dna$Total, method = "pearson") # not sig
cor.test(dna$total_Ca_mg, dna$Total, method = "spearman") # not sig

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
```
# 5. coda4microbime short data set oral hygiene score
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("../../rpoc/ps.RData")

#overall
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
#filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$oral_hygiene_score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.overall_oral.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.overall_oral.pdf")

# HI oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HI")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$oral_hygiene_score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HI_oral.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HI_oral.pdf")

# HEU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HEU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$oral_hygiene_score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HEU_oral.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HEU_oral.pdf")

# HUU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HUU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$oral_hygiene_score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HUU_oral.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUU_oral.pdf")
```
# 6. coda4microbiome long data set oral hygiene score
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
# overall t
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$Oral_Hygiene_Score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.overall_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.overall_oral.long.pdf")

# HI oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HI")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$Oral_Hygiene_Score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HI_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HI_oral.long.pdf")

# HEU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HEU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$Oral_Hygiene_Score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HEU_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HEU_oral.long.pdf")

# HUU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HUU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$cd4_count <- as.numeric(sample_data(temp)$cd4_count)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$Oral_Hygiene_Score)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HUU_oral.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUU_oral.long.pdf")
```
# 6. coda4microbiome long data set Ca
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(coda4microbiome)
#load data
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("~/long_oral/master_phyloseq.RData")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T, row.names=1)
ps.dat <- merge_phyloseq(ps.dat, sample_data(meta))
# overall t
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
temp <- subset_samples(temp, !is.na(sample_data(temp)$total_Ca_mg))
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.overall_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.overall_Ca.long.pdf")

# HI oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HI")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HI_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HI_Ca.long.pdf")

# HEU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HEU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HEU_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HEU_Ca.long.pdf")

# HUU oral hygiene score
temp <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[7])
# HI only, filter taxa that aren't found in at least 10% of all samples across visits, and at least 100 reads
temp <- subset_samples(temp, hiv_status == "HUU")
temp <- filter_taxa(temp, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_data(temp)$total_Ca_mg <- as.numeric(sample_data(temp)$total_Ca_mg)
temp
# save copy to reduce time on previous command
glom <- temp
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- as.data.frame(as.matrix(sample_data(glom))) # have to coerce to data frame
map <- tibble::rownames_to_column(map) # retain rownames for downstream processing
# get corresponding taxonomy name for each asv
taxa <- as(tax_table(glom), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- select(taxadf, V8)
orderdf <- orderdf %>%
    rownames_to_column(var = "ASV")
# rename ASV at species level
dat <- as.data.frame(dat)
dat <- dat %>% 
    rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))  
rownames(dat) <- paste(dat$V8, dat$ASV, sep="_")
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge new metadata with asv table so the response variable is in the same order
datmerge <- merge(dat, map, by.x = "row.names", by.y = "rowname")
datmerge <- datmerge[!duplicated(datmerge[c('Row.names')]), ]
row.names(datmerge) <- datmerge$Row.names
# define data and response variable
dif <- dim(datmerge)[2] - dim(map)[2]
x <- datmerge[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
dim(x)
# define response variable 
y <- as.numeric(datmerge$total_Ca_mg)
length(y)
# some stats
mean(y)
sd(y)

set.seed(852)
bal <- coda_glmnet(x, y, lambda = "lambda.min")
sum(bal$`log-contrast coefficients`)
#positive taxa
coef<-bal$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
bal$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
bal$taxa.name[negatives[on]]

pdf("./bal.HUU_Ca.long.pdf")
bal$`signature plot`
bal$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.HUU_Ca.long.pdf")
```
# 7. rpoC data for coda4microbiome species short data set (use the species identified by coda4microbiome using the long dataset)
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)

setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("~/rna_dohmain/rpoc/ps.RData")
set.seed(545433543)

# get rna data
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
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
homd <- read.table("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="", , fill = TRUE)
# remove any genes with a sum count of zero
subcount <- subcount[rowSums(subcount) != 0,]
# filter annotations by those found in our genecounts file
ann <- homd[homd$tag %in% rownames(subcount),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(subcount), rownames(ann)))]
subcount <- subcount[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(subcount)==rownames(ann)) # should all return true
# if all are true, merge together
subcount <- cbind(subcount, ann)
# get list of genera to pull from rpoC data later
subcount.gen <- unique(subcount$species)

# collapse by species and sum across rows
merge.count <- subcount %>% group_by(species) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# get total count across all rows at the species level 
group_sums <- group_sums %>% mutate(total=rowSums(across(where(is.numeric))))
# get total counts to use as denominator 
tot <- sum(group_sums$total)
thresh <- 0.00 * tot # less than 1% of total
# collapse low abundant groups
collapsed <- group_sums %>%
  mutate(species = ifelse(total < thresh, "Other", species)) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
# remove total column
collapsed <- subset(collapsed, select = -total)

# Create stacked histograms for each sample colored by species
# set color palette
# set color palette
gencols <- c(Catonella_morbi = "#F66140", 
      Porphyromonas_endodontalis = "#1FCAD4", 
      Prevotella_melaninogenica = "#97002E", 
      Treponema_vincentii = "#16AA48", 
      Porphyromonas_sp._oral_taxon_278 = "#0C5B6F", 
      Other = "#838383")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Catonella_morbi", 
                                              "Porphyromonas_endodontalis", 
                                              "Prevotella_melaninogenica", 
                                              "Treponema_vincentii", 
                                              "Porphyromonas_sp._oral_taxon_278"))) %>%
  filter(species %in% c("Catonella_morbi", 
                        "Porphyromonas_endodontalis", 
                        "Prevotella_melaninogenica", 
                        "Treponema_vincentii", 
                        "Porphyromonas_sp._oral_taxon_278"))

# Convert row names to a column
map <- metadata %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)
collapsed_long$hiv_status <- factor(collapsed_long$hiv_status, levels = c("HI", "HEU", "HUU"))

pdf("RNA_log10_balOHS.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = gencols) +
  facet_wrap(~hiv_status, scales = "free_x")
  # geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat RNA_log10_balOHS.pdf")


# # do it by oral hygiene score
# # long format
# collapsed_long <- collapsed %>%
#   pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# # Convert row names to a column
# map <- metadata %>%
#   rownames_to_column(var = "sample")
# # merge with metadata
# collapsed_long <- collapsed_long %>%
#   left_join(map, by = "sample")
# # order samples by HIV status
# collapsed_long <- collapsed_long %>%
#   mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))
# # group by oral hygiene rank
# collapsed_long <- collapsed_long %>%
#   group_by(species, hiv_status, oral_hygiene_score_remark) %>%
#   summarize(across(where(is.numeric), sum, na.rm = TRUE))


# # convert to log10 values
# collapsed_long <- collapsed_long %>%
#   mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# # reorder genera
# collapsed_long <- collapsed_long %>%
#   mutate(species = factor(species, levels = c("Catonella_morbi", 
#                                               "Porphyromonas_endodontalis", 
#                                               "Prevotella_melaninogenica", 
#                                               "Treponema_vincentii", 
#                                               "Porphyromonas_sp._oral_taxon_278"))) %>%
#   filter(species %in% c("Catonella_morbi", 
#                         "Porphyromonas_endodontalis", 
#                         "Prevotella_melaninogenica", 
#                         "Treponema_vincentii", 
#                         "Porphyromonas_sp._oral_taxon_278"))



# # get positions to add lines between samples
# sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

# pdf("RNA_log10_balOHS.pdf", width = 20, height = 5)
# ggplot(collapsed_long, aes(x = oral_hygiene_score_remark, y = log10_count, fill = species)) +
#   geom_bar(stat="identity", position="stack") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values = gencols) +
#   facet_wrap(~hiv_status, scales = "free_x")
#   # geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
# dev.off()
# system("~/.iterm2/imgcat RNA_log10_balOHS.pdf")
```
## 7. 2 Now do it for rpoC
```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
load("../../rpoc/ps.RData")
ps.dat.filt  <- subset_taxa(ps.dat, V8=="Catonella_morbi" | V8=="Porphyromonas_endodontalis" | V8=="Prevotella_melaninogenica" | V8 == "Treponema_vincentii")

# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.filt, taxrank = "V8")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))
# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V8)
asvtab <- as.data.frame(asvtab)
# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by genus (V25) and sum across rows
merge.asvtab <- asvtab %>% group_by(V94) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# rename columns
merge.asvtab <- merge.asvtab %>% rename(species = V94)

# groups to collapse by
keepgroup <- c("Catonella_morbi", "Porphyromonas_endodontalis", "Prevotella_melaninogenica", "Treponema_vincentii")
# Collapse and sum rows
collapsed <- merge.asvtab %>%
  mutate(species = ifelse(species %in% keepgroup, species, "Other")) %>%
  group_by(species) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))

# set color palette
gencols <- c(Catonella_morbi = "#F66140", 
      Porphyromonas_endodontalis = "#1FCAD4", 
      Prevotella_melaninogenica = "#97002E", 
      Treponema_vincentii = "#16AA48", 
      Porphyromonas_sp._oral_taxon_278 = "#0C5B6F", 
      Other = "#838383")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)

# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(species = factor(species, levels = c("Catonella_morbi", "Porphyromonas_endodontalis", "Prevotella_melaninogenica", "Treponema_vincentii")))

# Convert row names to a column
map <- metadata %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)
collapsed_long$hiv_status <- factor(collapsed_long$hiv_status, levels = c("HI", "HEU", "HUU"))


pdf("rpoC_log10_balOHS.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = species)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = gencols) +
  facet_wrap(~hiv_status, scales = "free_x")
  # geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("~/.iterm2/imgcat rpoC_log10_balOHS.pdf")
```
# 8. Deseq with OHS
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../02-pgap/gene_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
submap$oral_hygiene_score <- as.numeric(submap$oral_hygiene_score)
submap$oral_hygiene_score_scaled <- scale(submap$oral_hygiene_score, center = TRUE, scale = TRUE)

# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = subcount, colData = submap, design = ~oral_hygiene_score_scaled)
star_results <- star_results[rowSums(counts(star_results)) >= 50,]
star_results
# star_results$Ca_category <- factor(star_results$Ca_category, levels=c("low", "normal"))

# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm 
# normalize counts
norm_counts <- log2(counts(se_star, normalized = TRUE)+1)

res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.05:  362856"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 167304, 2.5%
# LFC < 0 (down)     : 195552, 3%
# outliers [1]       : 0, 0%
# low counts [2]     : 5355253, 81%
# (mean count < 1)
# HUU is positive, HEU cavity negative
resLFC <- lfcShrink(se_star, coef="total_Ca_scaled", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  362856"
# out of 6585836 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 201971, 3.1%
# LFC < 0 (down)     : 237241, 3.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 5355253, 81%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
write.table(resLFC, file="deseq_results_global-OHS.txt", quote=F, sep="\t")
save.image("deseq_results_global-OHS.RData")
```
# 9. Look at calculus and debris for Vince
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
library(viridis)
library(rstatix)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-inflammation")
meta <- read.csv("~/long_oral/map_domhain_long_2.txt", sep="\t", header=T)
hiv_stat <- c("HI", "HEU", "HUU")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
meta$hiv_status <- factor(meta$hiv_status, levels = hiv_stat)

wilcox_test(meta, Oral_Hygiene_Score ~ sex) %>%
  add_significance()
meta$Oral_Hygiene_Score <- as.numeric(meta$Oral_Hygiene_Score)
meta$total_Ca_mg <- as.numeric(meta$total_Ca_mg)
meta$pH <- as.numeric(meta$pH)
meta$flow_rate <- as.numeric(meta$flow_rate)
meta$cd4_count <- as.numeric(meta$cd4_count)

meta$anterior_calculus <- rowSums(meta[, c("Calculus_index_11", "Calculus_index_31")], na.rm = TRUE) / nrow(meta)
meta$posterior_calculus <- rowSums(meta[, c("Calculus_index_16", "Calculus_index_26", "Calculus_index_46", "Calculus_index_36")], na.rm = TRUE) / nrow(meta)
meta$anterior_gingival <- rowSums(meta[, c("Gingival_index_11", "Gingival_index_31")], na.rm = TRUE) / nrow(meta)
meta$posterior_gingival <- rowSums(meta[, c("Gingival_index_16", "Gingival_index_24", "Gingival_index_44", "Gingival_index_36")], na.rm = TRUE) / nrow(meta)
meta$hiv_status <- factor(meta$hiv_status, levels=c("HI", "HEU", "HUU"))
meta$Oral_Hygiene_Score_Remark <- factor(meta$Oral_Hygiene_Score_Remark, levels=c("Good", "Fair", "Poor"))

ant <- data.frame(sample_id = paste0(meta$study_id,"V",meta$visit_num))
ant$hiv_status <- meta$hiv_status
ant$location <- "anterior"
ant$calculus_index <- meta$anterior_calculus
ant$gingival_index <- meta$anterior_gingival
ant$Oral_Hygiene_Score_Remark <- meta$Oral_Hygiene_Score_Remark
ant$Oral_Hygiene_Score <- meta$Oral_Hygiene_Score

post <- data.frame(sample_id = paste0(meta$study_id,"V",meta$visit_num))
post$hiv_status <- meta$hiv_status
post$location <- "posterior"
post$calculus_index <- meta$posterior_calculus
post$gingival_index <- meta$posterior_gingival
post$Oral_Hygiene_Score_Remark <- meta$Oral_Hygiene_Score_Remark
post$Oral_Hygiene_Score <- meta$Oral_Hygiene_Score

compare_calculus <- unique(rbind(ant, post))
pdf("calculus_score.location.pdf")
ggplot(compare_calculus, aes(x = location, y = pmax(as.numeric(calculus_index), 0))) + 
    geom_point(aes(color = Oral_Hygiene_Score_Remark), position = position_jitter(width = 0.4, height = 0.0006)) + 
    scale_color_manual(values = c("Good" = "#006164", "Fair" = "#EDA247", "Poor" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns = TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    stat_summary(geom = "text", fun = "mean", aes(label = paste("Mean: ", scales::number(..y.., accuracy = 0.0001))), 
                  vjust = -40.5, size = 2, color = "black") +  
    scale_y_continuous(limits = c(0, 0.005), oob = scales::squish) +
    ylab("Calculus index") +
    xlab(NULL) +
    labs(color = "Oral Hygiene\nScore Remark") +
    facet_wrap(~hiv_status) +
    theme(legend.title = element_text(size=9), legend.title.align = 0.0) 
dev.off()
system("~/.iterm2/imgcat ./calculus_score.location.pdf")

pdf("calculus_score.hiv.pdf")
ggplot(compare_calculus, aes(x = hiv_status, y = pmax(as.numeric(calculus_index), 0))) + 
    geom_point(aes(color = Oral_Hygiene_Score_Remark), position = position_jitter(width = 0.4, height = 0.0006)) + 
    scale_color_manual(values = c("Good" = "#006164", "Fair" = "#EDA247", "Poor" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns = TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    stat_summary(geom = "text", fun = "mean", aes(label = paste("Mean: ", scales::number(..y.., accuracy = 0.0001))), 
                  vjust = -40.5, size = 2, color = "black") +  
    scale_y_continuous(limits = c(0, 0.005), oob = scales::squish) +
    ylab("Calculus index") +
    xlab(NULL) +
    labs(color = "Oral Hygiene\nScore Remark") +
    facet_wrap(~location) +
    theme(legend.title = element_text(size=9), legend.title.align = 0.0) 
dev.off()
system("~/.iterm2/imgcat ./calculus_score.hiv.pdf")

pdf("gingival_inflammation.location.pdf")
ggplot(compare_calculus, aes(x = location, y = pmax(as.numeric(gingival_index), 0))) + 
    geom_point(aes(color = Oral_Hygiene_Score_Remark), position = position_jitter(width = 0.4, height = 0.0006)) + 
    scale_color_manual(values = c("Good" = "#006164", "Fair" = "#EDA247", "Poor" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns = TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    stat_summary(geom = "text", fun = "mean", aes(label = paste("Mean: ", scales::number(..y.., accuracy = 0.0001))), 
                vjust = -40.5, size = 2, color = "black") +  
    scale_y_continuous(limits = c(0, 0.005), oob = scales::squish) +
    ylab("Gingival inflammation") +
    xlab(NULL) +
    labs(color = "Oral Hygiene\nScore Remark") +
    facet_wrap(~hiv_status) +
    theme(legend.title = element_text(size=9), legend.title.align = 0.0) 
dev.off()
system("~/.iterm2/imgcat ./gingival_inflammation.location.pdf")

pdf("gingival_inflammation.hiv_status.pdf")
ggplot(compare_calculus, aes(x = hiv_status, y = pmax(as.numeric(gingival_index), 0))) + 
    geom_point(aes(color = Oral_Hygiene_Score_Remark), position = position_jitter(width = 0.4, height = 0.0006)) + 
    scale_color_manual(values = c("Good" = "#006164", "Fair" = "#EDA247", "Poor" = "#DB4325")) +
    theme_bw() +
    # geom_hline(yintercept = c(1.2, 3.1), linetype = "dashed") +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns = TRUE, p.adjust.method = "fdr") +
    stat_summary(geom = "point", fun = "mean", size = 5, shape = 23, fill = "red") +
    stat_summary(geom = "text", fun = "mean", aes(label = paste("Mean: ", scales::number(..y.., accuracy = 0.0001))), 
                vjust = -40.5, size = 2, color = "black") +  
    scale_y_continuous(limits = c(0, 0.005), oob = scales::squish) +
    ylab("Gingival inflammation") +
    xlab(NULL) +
    labs(color = "Oral Hygiene\nScore Remark") +
    facet_wrap(~location) +
    theme(legend.title = element_text(size=9), legend.title.align = 0.0) 
dev.off()
system("~/.iterm2/imgcat ./gingival_inflammation.hiv_status.pdf")
```

