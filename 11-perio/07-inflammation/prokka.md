## 1.2 Composition of sample that is transcripts belong to red complex
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
counts <- read.csv("../03-global-diff/species_reads.prokka.txt", header=T,sep = "\t")
hivCols <- c("#8213A0", "#FA78FA", "#40A0FA")
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

pdf("redvoral.corr.species.prokka.pdf", width =20)
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
system("~/.iterm2/imgcat ./redvoral.corr.species.prokka.pdf")

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
  summarise(Total = sum((value), na.rm = TRUE))
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