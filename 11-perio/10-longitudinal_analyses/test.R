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
setwd("/home/suzanne/rna_dohmain/11-perio/10-longitudinal_analyses/")
virulence_genes_pg <- c("rgpA", "rgpB")
for (virulence2 in virulence_genes_pg) {
load("./deseq_results_V1-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
# add in annotations
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
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
# now make it for virulence factors
# get average for virulence
average_pg <- res_ord %>%
  filter(gene ==virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V1_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V1_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V2
load("./deseq_results_V2-HIvHUU.RData")
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
# for p. gingivalis
# # get avg_log2FoldChange
# sig_average_log2 <- sig_average %>%
#   select(genome, avg_log2FoldChange, species) %>%
#   pivot_longer(cols = "avg_log2FoldChange", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_baseMean
# sig_average_baseMean <- sig_average %>%
#   select(genome, avg_baseMean, species) %>%
#   pivot_longer(cols = "avg_baseMean", 
#                names_to = "metric", 
#                values_to = "value")

# # get avg_log2FoldChange
# sig_average_log2_sub<- filter(sig_average_log2, species == "Porphyromonas_gingivalis")
# sig_average_baseMean_sub<- filter(sig_average_baseMean, species == "Porphyromonas_gingivalis")

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V2_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V2_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")
# V3
load("./deseq_results_V3-HIvHUU.RData")
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
# now make it for virulence factors

# for p. gingivalis
res_ord_sub <- filter(res_ord, species == "Porphyromonas_gingivalis")
sig_average <- res_ord_sub %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# get average for virulence
average_pg <- res_ord %>%
  filter(gene %in% virulence2 & species == "Porphyromonas_gingivalis") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_V3_virus <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_V3_virus <- average_pg %>%
  select(genome, avg_baseMean_factor, species) %>%
  pivot_longer(cols = "avg_baseMean_factor", 
               names_to = "metric", 
               values_to = "value")


tr3 <- phylo4d(tree.root)
root_num <- getRoot(tree.root)

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

#base mean
combined_metrics<- bind_rows(
  sig_average_baseMean_V1_virus %>% mutate(type = "baseMean_V1"),
  sig_average_baseMean_V2_virus %>% mutate(type = "baseMean_V2"),
  sig_average_baseMean_V3_virus %>% mutate(type = "baseMean_V3"),
  sig_average_log2_V1_virus %>% mutate(type = "log2FoldChange_V1"),
  sig_average_log2_V2_virus %>% mutate(type = "log2FoldChange_V2"),
  sig_average_log2_V3_virus %>% mutate(type = "log2FoldChange_V3")
)
# fill in any missing genmoes
tips <- gsub("'","",rev(tree.root$tip.label))
complete_genomes <- expand.grid(
  genome = tip_labels,
  type = unique(combined_metrics$type)
)
combined_metrics_complete <- complete_genomes %>%
  left_join(combined_metrics, by = c("genome", "type")) %>%
  mutate(value = replace_na(value, 0))  # Replace NAs in 'value' with 0

combined_metrics_complete <- combined_metrics_complete %>%
  mutate(
    metric = case_when(
      grepl("^baseMean", type) ~ replace_na(metric, "avg_baseMean_factor"),
      TRUE ~ replace_na(metric, "avg_log2FoldChange_factor")
    )
  )


# Create the ggplot with the correct filtering for 'baseMean_V1', 'baseMean_V2', etc.
facet_labels <- c(
  "baseMean_V1" = "V1",
  "baseMean_V2" = "V2",
  "baseMean_V3" = "V3",
  "log2FoldChange_V1" = "V1",
  "log2FoldChange_V2" = "V2",
  "log2FoldChange_V3" = "V3"
)

heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V1", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "baseMean_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean V3", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V1"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Start a new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V2"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V2", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(type == "log2FoldChange_V3"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = -1, 
                       guide = guide_colorbar(title = "Average log fold change V3", title.position = "top")) +
  facet_wrap(~type, ncol =6, labeller = labeller(type = facet_labels) ) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-3,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_heat_", virulence2, ".pdf"), width = 8, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .45))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_heat_", virulence2 , ".pdf", sep=""))
}
