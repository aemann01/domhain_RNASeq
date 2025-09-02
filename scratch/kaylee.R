library(ggplot)
library(ggpubr)
setwd("~/rna_dohmain/scratch")

df <- read.csv("SA_plastic.csv")

pdf("kaylee_surfaces.pdf")
ggplot(df, aes(x=factor(Surfacetype),y=Average))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", method ="t_test", 
  hide.ns =TRUE, p.adjust.method = "fdr") + #adds signficance between the categories
  geom_boxplot() +
  theme_classic()
dev.off()

system("~/.iterm2/imgcat ./kaylee_surfaces.pdf")

pdf("kaylee_surfaces.no_correction.pdf")
ggplot(df, aes(x=factor(Surfacetype),y=Average))+
  geom_pwc(label = "{p.adj.format}{p.adj.signif}", method ="t_test", 
  hide.ns =TRUE, p.adjust.method = "none") + #adds signficance between the categories
  geom_boxplot() +
  theme_classic()
dev.off()

system("~/.iterm2/imgcat ./kaylee_surfaces.no_correction.pdf")

permanova_pairwise(df$average, grp=sample_data(ps.dat.clr)$study_group, method="euclidean") # check signficane
