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

load("~/rna_dohmain/rpoc/ps.RData")
#clr transformation     
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
#beta diversity
metadata <- as(sample_data(ps.dat.clr), "data.frame")
clr.dist <- dist(otu_table(ps.dat.clr), method="euclidean")
adonis2(clr.dist ~ hiv_status * aliquot_type * sex * age_y * delivery_type, data=metadata,  padj= "fdr")

permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$aliquot_type, method="euclidean", padj= "fdr")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$hiv_status, method="euclidean", padj= "fdr")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$sex, method="euclidean", padj= "fdr")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$age_y, method="euclidean", padj= "fdr")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$delivery_type, method="euclidean", padj= "fdr")
#make plots
hivCols <- c("#FA78FA", "#8213A0", "#40A0FA")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~hiv_status)
hiv <- plot_ordination(ps.dat.clr, ordcap, "samples", color="hiv_status") + 
    theme_minimal() + 
    geom_point(size=1.5) +
    scale_color_manual(values=hivCols)
dev.off()
pdf("./bdiv_cap.hiv_status.pdf")
print(hiv)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.hiv_status.pdf")

healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~aliquot_type)
aliq <- plot_ordination(ps.dat.clr, ordcap, "samples", color="aliquot_type") + 
    theme_minimal() + 
    geom_point(size=1.5) +
    scale_color_manual(values=healthCols)
pdf("./bdiv_cap.aliquot_type.pdf")
print(aliq)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.aliquot_type.pdf")

ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~age_y)
age <- plot_ordination(ps.dat.clr, ordcap, "samples", color="age_y") + 
	geom_point(size=1.5) + theme_minimal() +
	scale_colour_gradientn(colours = brewer.pal(9,"Spectral"))
pdf("bdiv_cap.age.pdf")
print(age)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.age.pdf")
[1] "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#E6F598" "#ABDDA4"
[8] "#66C2A5" "#3288BD"

sexCols <- c("#F8509C", "#3591E7")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~sex)
sex <- plot_ordination(ps.dat.clr, ordcap, "samples", color="sex") + 
	geom_point(size=1.5) + theme_minimal() +
	scale_color_manual(values=sexCols)
pdf("bdiv_cap.sex.pdf")
print(sex)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.sex.pdf")

birthCols <- c("#30BAF2", "#98002E")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~delivery_type)
birth <- plot_ordination(ps.dat.clr, ordcap, "samples", color="delivery_type") +
	 geom_point(size=1.5) + theme_minimal() +
	scale_color_manual(values=birthCols)
pdf("bdiv_cap.delivery_type.pdf")
print(birth)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.delivery_type.pdf")

#remove NAs
ps.sub <- subset_samples(ps.dat , membrane_rupture != "NA")
ps.dat.clr <- microbiome::transform(ps.sub, transform="clr", target="OTU")
metadata <- as(sample_data(ps.dat.clr), "data.frame")
clr.dist <- dist(otu_table(ps.dat.clr), method="euclidean")
adonis2(clr.dist ~ hiv_status * aliquot_type * sex * age_y * delivery_type* membrane_rupture, data=metadata)
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$membrane_rupture, method="euclidean", padj= "fdr")

memCols <- c("#86599B", "#4BA24F")
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~membrane_rupture)
mem <- plot_ordination(ps.dat.clr, ordcap, "samples", color="membrane_rupture")+
	 geom_point(size=1.5) + theme_minimal() +
	 scale_color_manual(values=memCols)
pdf("bdiv_cap.membrane_rupture.pdf")
print(mem)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.membrane_rupture.pdf")

#combine plots
combined <- plot_grid(aliq + theme(legend.position = c(0.8, 0.8)), 
	hiv +theme(legend.position = c(0.8, 0.8)), 
	age+theme(legend.position = c(0.8, 0.8)), 
	sex+theme(legend.position = c(0.8, 0.8)), 
	birth+theme(legend.position = c(0.8, 0.8)), 
	mem+theme(legend.position = c(0.8, 0.8)), 
	labels=c("a", "b", "c", "d", "e", "f"), ncol = 3, nrow = 2)
pdf("combine.pdf", width = 9.3, height =7)
print(combined)
dev.off()
system("~/.iterm2/imgcat ./combine.pdf")
svg("combine.svg")
print(combined)
dev.off()

#look at proportion of health
library(MASS)
library(vcd)
library(ggstatsplot)

meta_cats <- c("aliquot_type", "hiv_status", "sex", "delivery_type", "membrane_rupture")
for(var in meta_cats) {
print(var)
pdf(paste0("tooth.", var, ".pdf"))
print(ggbarstats(data = metadata,
               x = !!var, 
               y = tooth_health, 
               title = paste('Mosaic plot using ggstatsplot'),
               type = 'parametric',
               conf.level = 0.95,
               proportion.test = TRUE,
               ggtheme = ggplot2::theme_classic()))
dev.off()
}

#high strep mutans

ps.sub <- subset_taxa(ps.dat, V8 == "Streptococcus_mutans")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.sub, n=50)
ps_50 <- subset_taxa(ps.sub, rownames(tax_table(ps.sub)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./asv_dis.strep_mutans.pdf")
ggplot(data,aes(x=factor(Sample),y=Abundance,fill=factor(names))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
  #scale_fill_manual(values=healthCols)
dev.off()
system("~/.iterm2/imgcat ./asv_dis.strep_mutans.pdf")


#rarefied
rare_rpoc <- rarefy_even_depth(otu_table(ps.dat), rngseed = TRUE, replace = FALSE)
ps.rar.rpoc <- phyloseq(rare_rpoc, tax_rpoc, map) # create a phyloseq object
ps.sub <- subset_taxa(ps.rar.rpoc, V8 == "Streptococcus_mutans")
#more than 10 reads in 5% or more samples
#ps.rar.candida.its <- filter_taxa(ps.rar.candida.its, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
top50 <- top_taxa(ps.rar.rpoc, n=50)
ps_50 <- subset_taxa(ps.sub, rownames(tax_table(ps.sub)) %in% top50)
data <- psmelt(ps_50) # create dataframe from phyloseq object
#getting order for ggplot2
data$names <- paste(data$OTU, data$V8, sep="_")
data[c('cat1', 'cat2')] <- str_split_fixed(data$OTU, 'SV', 2)
data2 <- select(data, names, cat2)
data3 <- unique(data2)
order <- data3[order(as.numeric(as.character(data3$cat2))), ] #sort by ASV number
order_50 <- order$names
healthCols <- c("#AA0A3B", "#F87850", "#F0F032", "#2F5AC8", "#3FD2DC", "#24B45A")

pdf("./rare.strep_mutans.pdf")
ggplot(data,aes(x=factor(Sample),y=Abundance,fill=factor(names))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_minimal()+
  theme(axis.text.x=element_text(size=rel(0.5)))
  #scale_fill_manual(values=healthCols)
dev.off()
system("~/.iterm2/imgcat ./rare.strep_mutans.pdf")



#give order
sub_meta <- filter(metadata, !ketones == "NA")
sub_meta <- filter(sub_meta, !menopause == "pre")

urine_col <- c("colorless", "light_yellow", "yellow", "dark_yellow", "amber", "red")
urine_rank <- c(1,2,3,4,5,6)
sub_meta$rank_col <- replace(sub_meta$color_urine, sub_meta$color_urine %in% urine_col, urine_rank)

uti_ans <- c("y", "n")
uti_rep <- c(1,0)
sub_meta$uti_ans <- replace(sub_meta$uti, sub_meta$uti %in% uti_ans, uti_rep)
cor.test(as.numeric(sub_meta$rank_col),as.numeric(sub_meta$uti_ans))

#clairty
urine_col <- c("clear", "slightlycoudy", "cloudy", "turbid")
urine_rank <- c(1,2,3,4)
sub_meta$rank_clar <- replace(sub_meta$clarity_urine, sub_meta$clarity_urine %in% urine_col, urine_rank)
#blood
current <- c("negative", "trace", "small", "moderate", "large")
rep <- c(1,2,3,4,5)
sub_meta$rank_blood <- replace(sub_meta$blood, sub_meta$blood %in% current, rep)
cor.test(as.numeric(sub_meta$rank_blood),as.numeric(sub_meta$uti_ans))

#leuks
current <- c("negative", "trace", "small", "moderate", "large")
rep <- c(1,2,3,4,5)
sub_meta$rank_leuk <- replace(sub_meta$leukocyte_esterase, sub_meta$leukocyte_esterase %in% current, rep)
test <- cor.test(as.numeric(sub_meta$rank_leuk),as.numeric(sub_meta$uti_ans))
p.adjust(test, method = "fdr", n = length(test))

#rcorr
rank_meta <- as.matrix(dplyr::select(sub_meta, uti_ans, rank_col, rank_clar, rank_blood, rank_leuk, pH, estrogen))
cor_meta <- rcorr(rank_meta,type = c("pearson"))
cor_meta.adj <- rcorr_padjust(cor_meta, method = "BH")
cor_meta.adj


#balance of taxa
library(grid)
library(coda4microbiome)
library(tidyverse)
# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.dat, "V8")
sub <- subset_samples(glom, mouth_health == "CF" | mouth_health == "CD")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
sub <- filter_taxa(sub, function(x) sum(x > 15) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(sub)))
map <- sample_data(sub)

# get corresponding taxonomy name for each asv
taxa <- as(tax_table(sub), "matrix")
taxadf <- as.data.frame(taxa)
orderdf <- dplyr::select(taxadf, V8)
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

# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,2:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
x# response variable
y <- as.factor(dat$mouth_health)
y
# z <- data.frame(Tooth_Classification = as.factor(dat$sex)) #possible cofound

geo_its <- coda_glmnet(x=x,y=y)
geo_its
sum(geo_its$`log-contrast coefficients`)

#positive taxa
coef<-geo_its$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
geo_its$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
geo_its$taxa.name[negatives[on]]

pdf("./bal.tooth.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.tooth.pdf")

pdf("./bal.hiv.pdf")
geo_its$`signature plot`
geo_its$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./bal.hiv.pdf")


#diff abundance
library(corncob)
library(magrittr)
#glom to species level
its.dat.sp <- ps.dat %>%
                tax_glom("V8") %>% subset_samples(hiv_status != "HEU") 
its.dat.sp

#diff abundance
da_analysis_its <- differentialTest(formula = ~ mouth_health,
                               phi.formula = ~ mouth_health,
                               formula_null = ~ 1,
                               phi.formula_null = ~ mouth_health,
                               test = "Wald",
                               boot = FALSE,
                               data = its.dat.sp,
                               fdr_cutoff = 0.05)
da_analysis_its
#look at sign taxa
da_analysis_its$significant_taxa
pdf("./diffab.tooth_health.pdf", width = 20)
plot(da_analysis_its, level=c("V8"))
dev.off()
