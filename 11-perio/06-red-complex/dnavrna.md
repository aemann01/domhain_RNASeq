# 3. DNA vs RNA comparison
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
setwd("~/rna_dohmain/11-perio/06-red-complex")
#get relative abundance of dna
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
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
#get relative abundance of rna
rna_counts <- read.csv("~/rna_dohmain/09-urease/09-global-distro/species_reads.txt", sep='\t')
red_rna <- select(rna_counts, sample, Tannerella_forsythia, Porphyromonas_gingivalis, Treponema_denticola)
red_rna$nucl <- "rna"
red_rna <- melt(red_rna)
red_rna <- red_rna %>%
  rename(species = variable)
map$sample <- row.names(map)
meta <- select(map, sample, hiv_status)
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

```
Porportion of samples with DNA
```R 
#proportion of samples with P. gingivalis
# proportion of different types of strep in pd vs pf
rel <- microbiome::transform(ps.dat, "compositional")
strep <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
# plot
filtered_data <- data %>% filter(Abundance > 20)

# Calculate total number of samples per HIV status
total_samples_per_hiv <- filtered_data %>%
  group_by(hiv_status) %>%
  summarize(total_samples = n_distinct(Sample))

# Calculate count of samples per species and HIV status
species_counts_by_hiv <- filtered_data %>%
  group_by(V8, hiv_status) %>%
  summarize(sample_count = n_distinct(Sample), .groups = 'drop')

# Join total counts with species counts and calculate percentages
percentages <- species_counts_by_hiv %>%
  left_join(total_samples_per_hiv, by = "hiv_status") %>%
  mutate(percentage = (sample_count / total_samples) * 100)

# Print the results
print(percentages)
percentages$hiv_status <- factor(percentages$hiv_status, levels = c("HI", "HEU", "HUU"))
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

pdf("samples.red.prop.pdf")
ggplot() + geom_bar(data = percentages, aes(x = hiv_status, y = percentage, fill = hiv_status), position = "dodge", stat = "identity")+
	facet_grid(~ V8, switch = "x", scales = "free_x")+
	scale_color_manual(values=hivCols)+
	scale_fill_manual(values=hivCols)+
	theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./samples.red.prop.pdf")
#find signficance
percentages$none <- percentages$total_samples - percentages$sample_count
sub_precent <- as.data.frame(select(percentages, V8, hiv_status, sample_count, none))
sub_precent <- filter(sub_precent, V8 == "Porphyromonas_gingivalis" & hiv_status != "HEU") %>% select(hiv_status, sample_count, none)
row.names(sub_precent) <- sub_precent$hiv_status
sub_precent <- select(sub_precent, sample_count, none)
chisq.test(sub_precent)
```
Porportion of all 1960 samples with DNA
```R 
#proportion of samples with P. gingivalis
# proportion of different types of strep in pd vs pf
load("~/long_oral/master_phyloseq.RData")
rel <- microbiome::transform(ps.dat, "compositional")
strep <- subset_taxa(ps.dat, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$Sample <- factor(data$Sample, levels=unique(data$Sample))
# plot
filtered_data <- data %>% filter(Abundance > 20)

# Calculate total number of samples per HIV status
total_samples_per_hiv <- filtered_data %>%
  group_by(hiv_status) %>%
  summarize(total_samples = n_distinct(Sample))

# Calculate count of samples per species and HIV status
species_counts_by_hiv <- filtered_data %>%
  group_by(V8, hiv_status) %>%
  summarize(sample_count = n_distinct(Sample), .groups = 'drop')

# Join total counts with species counts and calculate percentages
percentages <- species_counts_by_hiv %>%
  left_join(total_samples_per_hiv, by = "hiv_status") %>%
  mutate(percentage = (sample_count / total_samples) * 100)

# Print the results
print(percentages)
percentages$hiv_status <- factor(percentages$hiv_status, levels = c("HI", "HEU", "HUU"))
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

pdf("all_samples.red.prop.pdf")
ggplot() + geom_bar(data = percentages, aes(x = hiv_status, y = percentage, fill = hiv_status), position = "dodge", stat = "identity")+
  facet_grid(~ V8, switch = "x", scales = "free_x")+
  scale_color_manual(values=hivCols)+
  scale_fill_manual(values=hivCols)+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./all_samples.red.prop.pdf")

strep <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8 == "Treponema_denticola" | V8 == "Tannerella_forsythia")
glom <- tax_glom(strep, taxrank=rank_names(strep)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data$hiv_status <- factor(data$hiv_status, levels = c("HI", "HEU", "HUU"))
data$Abundance > 0
data_more <- data[data$Abundance > 0,]
# plot
pdf("all_samples.red.distro.pdf")
ggplot(data_more, aes(x = Abundance, fill = hiv_status)) +
  geom_density(position ="dodge", alpha = 0.4) +  # Adjust binwidth as needed
  facet_grid(~ V8, scales ="free") +
  scale_fill_manual(values = hivCols) + 
  theme_minimal() +
  labs(x = "Abundance", y = "Count", fill = "HIV Status")
dev.off()
system("~/.iterm2/imgcat ./all_samples.red.distro.pdf")

# Group by V8 and hiv_status, then summarize
# Compute summary statistics
summary_stats <- data %>%
  group_by(V8, hiv_status) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    Median = median(Abundance, na.rm = TRUE),
    # Mode = mode_function(Abundance),
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )
print(summary_stats)
summary_stats <- data_more %>%
  group_by(V8, hiv_status) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    Median = median(Abundance, na.rm = TRUE),
    # Mode = mode_function(Abundance),
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

#find signficance
percentages$none <- percentages$total_samples - percentages$sample_count
sub_precent <- as.data.frame(select(percentages, V8, hiv_status, sample_count, none))
sub_precent <- filter(sub_precent, V8 == "Porphyromonas_gingivalis" & hiv_status != "HEU") %>% select(hiv_status, sample_count, none)
row.names(sub_precent) <- sub_precent$hiv_status
sub_precent <- select(sub_precent, sample_count, none)
chisq.test(sub_precent)

```
# 2. Compare DeSeq2 outputs
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(reshape2)
library(ggpubr)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-red-complex")
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("./red_counts.txt", header=T, sep="\t", row.names=1)
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
vld <- log2(counts(se_star, normalized = TRUE)+1)
trans <- melt(vld)
homd <- read.table("red_annots.txt", header=T, sep="\t", quote="")
comb_rna <- left_join(trans, homd, by = c("Var1"="tag"))
#combine info
red_rna <- comb_rna %>% filter(species == "Tannerella_forsythia" | species == "Porphyromonas_gingivalis" | species == "Treponema_denticola")
rna_change1 <- red_rna %>%
  group_by(Var2, species) %>%
  summarize(log2fold = mean(value), .groups = 'drop')
rna_change1
```
```R
#dna
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1)) # get rid of weird empty column in genecounts
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HUU" | metadata$hiv_status == "HEU" | metadata$hiv_status == "HI",]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI","HEU", "HUU"))
# run deseq
ptm <- proc.time()
se_star <- DESeq(star_results, fitType="local")
proc.time() - ptm #pcoa diversity
# vld <- varianceStabilizingTransformation(se_star)
# trans <- melt(assay(vld))
vld <- log2(counts(se_star, normalized = TRUE)+1)
trans <- melt(vld)
homd <- read.table("../../rpoc/taxonomy_bac.txt", header=F, sep="\t")
comb_dna <- left_join(trans, homd, by = c("Var1"="V1"))
#combine info
red_dna <- comb_dna %>% filter(V8 == "Tannerella_forsythia" | V8 == "Porphyromonas_gingivalis" | V8 == "Treponema_denticola")
dna_change <- red_dna %>%
  group_by(Var2, V8) %>%
  summarize(log2fold = mean(value), .groups = 'drop')
dna_change

comb1 <- left_join(dna_change, rna_change1, by = c("Var2"="Var2", "V8"="species"))
hiv <- select(metadata, sample_id, hiv_status)
comb1 <- left_join(comb1, hiv, by = c("Var2"="sample_id"))
comb1
pdf("red.log2.corr.pdf")
ggscatter(comb1, x = "log2fold.x", y = "log2fold.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 facet_wrap(V8 ~hiv_status)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.log2.corr.pdf")

dna_change$nucl <- "dna"
colnames(dna_change) <- c("sample", "species", "Log2Fold", "nucl")
rna_change <- red_rna %>%
  group_by(Var2, species) %>%
  summarize(log2fold = mean(value), .groups = 'drop')
rna_change1$nucl <- "rna"
colnames(rna_change1) <- c("sample", "species", "Log2Fold", "nucl")

combined_df <- rbind(dna_change,rna_change1)
combined_df <- left_join(combined_df, hiv, by = c("sample"="sample_id"))

pdf("red.RNAvDNA.bar.pdf")
ggplot() + geom_bar(data = combined_df, aes(x = sample, y = Log2Fold, fill = nucl), position = "dodge", stat = "identity")+
  facet_grid(species ~ hiv_status, switch = "x", scales = "free_x")+
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.RNAvDNA.bar.pdf")
```
RNA/DNA ration
```R
viru_rna <- red_rna %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")

rna_change1 <- viru_rna %>%
  group_by(Var2, species) %>%
  summarize(log2fold = mean(value), .groups = 'drop')
rna_change1

comb1 <- left_join(dna_change, rna_change1, by = c("sample"="Var2", "species"="species"))
hiv <- select(metadata, sample_id, hiv_status)
comb1 <- left_join(comb1, hiv, by = c("sample"="sample_id"))
comb1
colnames(comb1) <- c("sample", "species", "dna", "nucl", "rna","hiv_status")

comb1$ratio <- comb1$rna / comb1$dna
comb1$hiv_status <- factor(comb1$hiv_status, levels = c("HI", "HEU", "HUU"))
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")

pdf("red.ratio.box.pdf")
ggplot(comb1, aes(x=factor(hiv_status),y=ratio))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") +
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+
  labs(x ="Sample", y = "RNA/DNA")+
  facet_wrap(~species)+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./red.ratio.box.pdf")

metadata$hiv_status <- factor(metadata$hiv_status, levels = c("HI", "HEU", "HUU"))
pdf("hiv.age.box.pdf")
ggplot(metadata, aes(x=factor(hiv_status),y=age_y))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", hide.ns =FALSE, p.adjust.method = "fdr") +
  geom_boxplot() +
  geom_jitter(aes(color=hiv_status), shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = hivCols)+
  labs(x ="Sample", y = "Age")+
  theme_classic()
dev.off()
system("~/.iterm2/imgcat ./hiv.age.box.pdf")
```
RNA vs Species using also orange complex
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(reshape2)
library(ggpubr)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/07-orange-complex")
resdf <- read.table("./annotated_HIvHUU_deseq_results.txt", header=T, sep="\t", quote="")

resdf$species <- paste0(resdf$Genus,"_", resdf$Species)
red_rna <- resdf %>% filter(species == "Tannerella_forsythia" | species == "Porphyromonas_gingivalis" | species == "Treponema_denticola" | species=="Campylobacter_gracilis" | species=="Campylobacter _rectus" | species=="Campylobacter_showae" | species=="Eubacterium_nodatum" | species=="Fusobacterium_nucleatum" | species=="Fusobacterium_periodonticum" | species=="Prevotella_micros" | species=="Prevotella_intermedia" | species=="Prevotella_nigrescens"| species=="Streptococcus_constellatus")

rna_change <- red_rna %>%
  group_by(species) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
rna_change

#dna
resLFC <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/deseq_results_rpoC-HUUvHI.txt", sep="\t", header=T, row.names=1 )
# add in annotations
homd <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/rpoC.annotations.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]
# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
#combine name 
resdf$species <- paste0(resdf$Genus,"_",resdf$Species)

dna_change <- resdf %>%
  group_by(species) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
dna_change <- as.data.frame(dna_change)
red_dna <- dna_change %>% filter(species == "Tannerella_forsythia" | species == "Porphyromonas_gingivalis" | species == "Treponema_denticola" | species=="Campylobacter_gracilis" | species=="Campylobacter _rectus" | species=="Campylobacter_showae" | species=="Eubacterium_nodatum" | species=="Fusobacterium_nucleatum" | species=="Fusobacterium_periodonticum" | species=="Prevotella_micros" | species=="Prevotella_intermedia" | species=="Prevotella_nigrescens"| species=="Streptococcus_constellatus")

#combine data
comb1 <- left_join(red_dna, rna_change, by = c("species"="species"))
pdf("red.log2.corr.pdf")
ggscatter(comb1, x = "log2fold.x", y = "log2fold.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.log2.corr.pdf")
```






# 3. Corncob for all rpoc
```R
library(phyloseq)
library(corncob)
library(magrittr)

set.seed(12349)
setwd("/hone/suzannerna_dohmain/11-perio/01-rpoC-diff-abundance")
load("~/long_oral/master_phyloseq.RData")


ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HI")
da_analysis <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.HUUvHI.all_rpoc.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHI.all_rpoc.pdf")

ps.dat.sub <- subset_samples(ps.dat, hiv_status == "HUU" | hiv_status == "HEU")
da_analysis1 <- differentialTest(formula = ~ hiv_status,
                               phi.formula = ~ hiv_status,
                               formula_null = ~ 1,
                               phi.formula_null = ~ hiv_status,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sub,
                               fdr_cutoff = 0.05)
da_analysis1
#look at sign taxa
da_analysis1$significant_taxa
pdf("./diffab.HUUvHEU.all_rpoc.pdf", width = 20)
plot(da_analysis1, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.HUUvHEU.all_rpoc.pdf")
save.image("ASV_corncob.RData")
```


```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
resLFC <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/deseq_results_rpoC-HUUvHI.txt", sep="\t", header=T, row.names=1 )
# add in annotations
homd <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/rpoC.annotations.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
#combine name 
resdf$GeneInfo <- paste0(resdf$Genus,"_",resdf$Species)

dna_change <- resdf %>%
  group_by(GeneInfo) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
dna_change <- as.data.frame(dna_change)
red_dna <- dna_change %>% filter(GeneInfo == "Tannerella_forsythia" | GeneInfo == "Porphyromonas_gingivalis" | GeneInfo == "Treponema_denticola")

#get rna changes for groups of interest
HI <- read.csv("../06-red-complex/deseq_results_red-HIvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
# add in annotations
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
HI <- HI %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")

rna_change <- HI %>%
  group_by(species) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
rna_change

#combine data
comb1 <- left_join(red_dna, rna_change, by = c("GeneInfo"="species"))
pdf("red.log2.corr.pdf")
ggscatter(comb1, x = "log2fold.x", y = "log2fold.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.log2.corr.pdf")
```
HEU vs HUU
```R
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
resLFC <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/deseq_results_rpoC-HUUvHI.txt", sep="\t", header=T, row.names=1 )
# add in annotations
homd <- read.table("~/rna_dohmain/11-perio/01-rpoC-diff-abundance/rpoC.annotations.txt", header=T, sep="\t", quote="")
# filter by locus tag 
resdf <- as.data.frame(resLFC)
ann <- homd[homd$locus_tag %in% rownames(resdf),]

# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(resdf), rownames(ann)))]

resdf <- resdf[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(resdf)==rownames(ann)) # should all return true
# if all are true, merge together
resdf <- cbind(resdf, ann)
#combine name 
resdf$GeneInfo <- paste0(resdf$Genus,"_",resdf$Species)

dna_change <- resdf %>%
  group_by(GeneInfo) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
dna_change <- as.data.frame(dna_change)
red_dna <- dna_change %>% filter(GeneInfo == "Tannerella_forsythia" | GeneInfo == "Porphyromonas_gingivalis" | GeneInfo == "Treponema_denticola")

#get rna changes for groups of interest
HI <- read.csv("../06-red-complex/deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
# add in annotations
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
HI <- HI %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")

rna_change <- HI %>%
  group_by(species) %>%
  summarize(log2fold = mean(log2FoldChange), .groups = 'drop')
rna_change

#combine data
comb2 <- left_join(red_dna, rna_change, by = c("GeneInfo"="species"))
pdf("red.log2.corr.pdf")
ggscatter(comb2, x = "log2fold.x", y = "log2fold.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.log2.corr.pdf")
```
Combine
```R
comb <- rbind(comb1, comb2)
pdf("red.log2.corr.pdf")
ggscatter(comb, x = "log2fold.x", y = "log2fold.y",
   add = "reg.line", conf.int = TRUE, cor.method="spearman")+
 stat_cor()+
 # facet_wrap(~species)+
 theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./red.log2.corr.pdf")
