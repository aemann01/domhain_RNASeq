# 1. Heatmap
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(microbiome)
#load data
setwd("/home/suzanne/rna_dohmain/11-perio/05-virulence")
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
pdf("pca.hiv.red.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca.hiv.red.pdf")

#make dendogram
#get virulence factors ids
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
homd_sub <- homd %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA" | gene == "tfsA" | gene == "tfsB")
viru_genes <- homd_sub$tag
#Get top varying genes
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
# topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 50)
 
#make a subset of the log transformed counts for just the top 25 varying genes
# topCounts <- assay(vld)[topVarGenes,]
virus_vld <- assay(vld)[viru_genes,]
top_viru <- head(order(rowVars(virus_vld), decreasing=TRUE), 22580)
topCounts <- virus_vld[top_viru,]

df <- as.data.frame(colData(vld)[,c("hiv_status","age_y")])

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#give colors
mat_colors <- list(hiv_status = c("#8213A0","#FA78FA", "#40A0FA"))
names(mat_colors$hiv_status) <- c("HI", "HEU","HUU")
# mat_colors$tooth_health <- c("#F35E5A","#25AE2B", "#4F87FF")
# names(mat_colors$tooth_health) <- c("D", "E", "H")
#get relative abundance of red complex
seqtab <- t(read.table("../../rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("../../rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("../../homd_map/map.txt", sep="\t", header=T, row.names=1)
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
rel <- microbiome::transform(glom, "log10")
sub <- subset_taxa(rel, V8=="Porphyromonas_gingivalis" | V8=="Tannerella_forsythia" | V8=="Treponema_denticola")
abund <- t(as.data.frame(otu_table(sub)))
abund_df <- as.data.frame(abund)
abund_with_sum <- abund_df %>%
  rowwise() %>%
  mutate(Sum = sum(c_across(everything()))) %>%
  ungroup()
abund <-as.data.frame(abund_with_sum$Sum)
row.names(abund) <- row.names(t(as.data.frame(otu_table(sub))))

merged_df <- merge(df, abund, by = "row.names", all = FALSE)
colnames(merged_df) <- c("sample","hiv_status", "age", "rel_abundance")
row.names(merged_df) <- merged_df$sample
merged_df <- merged_df[, -1]
merged_df$rel_abundance <- as.numeric(merged_df$rel_abundance)
#get row annotaions
rows <- as.data.frame(row.names(topCounts))
colnames(rows) <- c('seq')
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  tag))
rows2$name <- paste0(rows2$species, "_", rows2$Gene)
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(topCounts) <- rows$name

x <- pheatmap(topCounts, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),fontsize_row= 4, fontsize_col =4)
save_pheatmap_pdf(x, "heatmap.hiv.viru.red.pdf")
system("~/.iterm2/imgcat ./heatmap.hiv.viru.red.pdf")

#collapse same gene
library(reshape2)
library(tidyr)
# gene_counts <- assay(vld)
gene_counts <- assay(vld)[viru_genes,]

#get row annotaions
rows <- as.data.frame(row.names(gene_counts))
colnames(rows) <- c('seq')
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
#get the names
homd$Gene <- gsub(x = homd$gene, pattern = "_[0-9]", replacement = "")
homd$Gene <- gsub(x = homd$Gene, pattern = "[0-9]", replacement = "")
rows2 <-left_join(rows, homd, by = join_by(seq ==  tag))
# rows2$name <- paste0(rows2$Genus, "_", rows2$Species, "_",rows2$Gene)
rows2$name <- paste0(rows2$species, "_", rows2$Gene)
rows2 <- rows2 %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
row <- rows2[, c("seq", "name")]
rows <- as.data.frame(row[, !(colnames(row) %in% "seq")])
rownames(rows) <- rows2$seq
colnames(rows) <- c('name')
row.names(gene_counts) <- rows$name
long_df <- melt(gene_counts)
group_all <- long_df %>% group_by(Var2, Var1) %>% summarise(average_baseMean = mean(value))
wide_df <- pivot_wider(group_all, id_cols = Var1,  names_from = Var2,  values_from = average_baseMean)
wide_df$sum <-rowSums(wide_df[, -1])
wide_df <- wide_df[order(-wide_df$sum), ]
df <- head(wide_df, 100)
df <- df[, -ncol(df)]
df <- as.data.frame(df)
row.names(df) <- df$Var1
df <- df[, -1]
x <- pheatmap(df, annotation_col = merged_df, annotation_colors = mat_colors,color = brewer.pal(9, "Greys"),)
save_pheatmap_pdf(x, "heatmap.hiv.viru_gene.pdf")
system("~/.iterm2/imgcat ./heatmap.hiv.viru_gene.pdf")
```
# 2. Donut Plot
Reads belonging to each HIV group
```R
library(tidyverse)
library(reshape2)
library(webr)

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/05-virulence")
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
#get annots for genes
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(genecounts),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(genecounts), rownames(ann)))]
genecounts <- genecounts[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(genecounts)==rownames(ann)) # should all return true
# if all are true, merge together
genecounts <- cbind(genecounts, ann)
gene_long <- melt(genecounts)
combined_data <- left_join(gene_long, metadata, by = c("variable" = "sample_id"))

combined_sub <- combined_data %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
row <- rows2[, c("seq", "name")]
#make donut plots
gene_sum <- combined_sub %>% 
  group_by(hiv_status, genome, species) %>% 
  summarise(Total_genome = sum(value, na.rm = TRUE))
gene_sum2 <- gene_sum %>% 
  group_by(hiv_status, species) %>% 
  summarise(Total = mean(Total_genome, na.rm = TRUE))

pdf("viru.donut.pdf")
PieDonut(gene_sum2, aes("hiv_status", "species", count="Total"), showRatioThreshold = F)
dev.off()
system("~/.iterm2/imgcat ./viru.donut.pdf")
```
# 3. Upset plot
```R
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ggplot2)

setwd("~/rna_dohmain/11-perio/05-virulence")
HI <- read.csv("../04-red-diff/deseq_results_red-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("../04-red-diff/deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
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
#HEU
HEU$hiv_status <- "HEU"
#get annots for HEEU
homd <- read.table("../02-pgap/red_annots.txt", header=T, sep="\t", quote="")
# filter by locus tag 
ann <- homd[homd$tag %in% rownames(HEU),]
# reorder
rownames(ann) <- ann$tag
sortrow <- rownames(ann)[order(match(rownames(HEU), rownames(ann)))]
HEU <- HEU[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]
# check that locus tags match between the two dataframes
table(rownames(HEU)==rownames(ann)) # should all return true
# if all are true, merge together
HEU <- cbind(HEU, ann)

#bind together HI and HEU
resdf <- rbind(HI, HEU)

#get only signficant genes
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of -2 since that is the side from their respective deseq analysis
resdf <- resdf %>%
  mutate(mark = if_else(log2FoldChange <= -lfc & padj <= pval, "1", "0"))
resdf$mark <- as.numeric(resdf$mark)
resdf_filtered <- select(resdf, mark, hiv_status, gene_tag, species)
resdf_wide <- resdf_filtered %>%
  pivot_wider(names_from = hiv_status, values_from = mark, values_fill = list(mark = 0)) %>%
  arrange(gene_tag)
#make upset plot
resdf_df <- as.data.frame(resdf_wide)
resdf_clean_df <- resdf_df %>%
  filter(!is.na(HI) & !is.na(HEU))
resdf_clean_df$species <- as.factor(resdf_clean_df$species)
#get just virulence genes
homd_sub <- homd %>% filter(gene == "kgp" | gene == "rgpB" | gene == "rgpA" | gene == "hagA" | gene == "fimA" | gene== "serC" | gene =="folD" | gene== "fhs" | gene =="sda" | gene == "susB" | gene == "kly" | gene == "eno" | gene == "hagA" | gene == "fimA" | gene == "oppA" | gene == "prtP"  | gene == "flaA" | gene == "flaB"| gene == "fliE" | gene == "cheX" | gene == "cheY" | gene == "hbpA" | gene == "hbpB" | gene == "troA" | gene =="bspA")
viru_genes <- homd_sub$tag
filtered_df <- resdf_clean_df[resdf_clean_df$gene_tag %in% viru_genes, ]

pdf("HI_HEU.upset.viru.pdf")
upset(filtered_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially Expressed Genes", sets.x.label="Total Number of Differentially Expressed Genes",
	queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()
system("~/.iterm2/imgcat ./HI_HEU.upset.viru.pdf")
```
# 4. Virulence db
```sh
mkdir databases && cd databases
wget https://github.com/Wanting-Dong/MetaVF_toolkit/raw/refs/heads/main/databases/eVFGC.fasta.gz
wget https://github.com/Wanting-Dong/MetaVF_toolkit/raw/refs/heads/main/databases/pVFGC.fasta.gz
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/databases/VF_info_file
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/databases/VFs_complete_num

gzip -d *gz
cd ../
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/metaVF.py
mkdir workflows && cd workflows
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/workflows/Snakefile_draft
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/workflows/Snakefile_PE
cd ../
mkdir scripts && cd scripts
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/filter_blast_result.py
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/select_best_hit.py
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/cal_VF_genes_draft.py
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/cal_VF_genes.py
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/cal_VF_draft.py
wget https://raw.githubusercontent.com/Wanting-Dong/MetaVF_toolkit/refs/heads/main/scripts/cal_VF.py
cd ../
metaVF.py -h

# copy red complex files
cp $(grep -i "Treponema_denticola\|Porphyromonas_gingivalis\|Tannerella_forsythia" ~/rna_dohmain/11-perio/02-pgap/gene_annots.txt | awk '{print $1}' | sed 's/_.*//' | sort | uniq | sed 's/SEQ/\/home\/suzanne\/rna_dohmain\/11-perio\/02-pgap\/annotations\/SEQ/' | sed 's/\.1/\.1\.fna/') ./test

metaVF.py -p ~/rna_dohmain/11-perio/05-virulence -pjn draft_test -id ./test -o ./test_draft -m draft -c 20 -ti 90 -tc 80
```
Pathofact2
```sh
git clone https://gitlab.lcsb.uni.lu/ESB/PathoFact2
cd PathoFact2
wget https://zenodo.org/records/14192463/files/DATABASES.tar.gz
tar -xvf DATABASES.tar.gz    
rm DATABASES.tar.gz 
bash install_key_envs.sh -d=~/rna_dohmain/11-perio/05-virulence/PathoFact2
sed -i 's/\/abs\/path\/to/\/home\/suzanne\/rna_dohmain\/11-perio\/05-virulence/' ./Config.yaml
conda activate PathoFact_env
cp $(grep -i "Treponema_denticola\|Porphyromonas_gingivalis\|Tannerella_forsythia" ~/rna_dohmain/11-perio/02-pgap/gene_annots.txt | awk '{print $1}' | sed 's/_.*//' | sort | uniq | sed 's/SEQ/\/home\/suzanne\/rna_dohmain\/11-perio\/02-pgap\/annotations\/SEQ/' | sed 's/\.1/\.1\.fna/') ./red_files
ls *fna | sed 's/.fna//' | while read line; do mv $line.fna $line.fasta; done
# give unique identifier
for f in *.fasta; do awk -v file="$f" '{if ($0 ~ /^>/) {print ">" file "_" NR} else {print}}' "$f" > "${f%.fasta}_unique_ids.fasta"; done
mkdir old
ls *fasta | grep -v unique | while read line; do mv $line ./old; done

nano list_of_samples.csv
bash run_PathoFact.sh "-p --printshellcmds"
bash run_PathoFact.sh "-n -p"
bash run_PathoFact.sh "-"
```
# 5. P. ging tree
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
setwd("/home/suzanne/rna_dohmain/11-perio/05-virulence/")
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
load("../03-global-diff/deseq_results-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/p_gingivalis/align/pging.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
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
  xlim(0, .035)
tip_labels <- rev(get_taxa_name(tree_plot))
# make log10 for RNA
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
merge.count <- subcount %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# long format
collapsed_long <- group_sums %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
pgings_seqs <- unique(sort(homd[homd$species  == "Porphyromonas_gingivalis",]$genome))

collapsed_long2 <- collapsed_long[collapsed_long$genome %in% pgings_seqs,]
rna_abund <- left_join(collapsed_long2, submap, by = c("sample" = "sample_id"))

summary_log10_counts <- rna_abund %>%
  group_by(hiv_status, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))

#loop
virulence_genes_pg <- c("rgpA", "rgpB", "hagA", "rpoC")
for (virulence2 in virulence_genes_pg) {

# add in annotations
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
sig_average_log2_200days <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_200days <- average_pg %>%
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
  xlim(0, .055)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
summary_log10_counts$metric <- "log10_count"
combined_metrics<- bind_rows(
  sig_average_baseMean_200days %>% mutate(type = "baseMean"),
  sig_average_log2_200days %>% mutate(type = "log2FoldChange"),
  summary_log10_counts %>% mutate(type = "log10(count)")
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
  "baseMean" = "baseMean",
  "log2FoldChange" = "log2FoldChange",
  "log10(count)" = "log10(count)"
)
combined_metrics_complete[combined_metrics_complete$genome  == "SEQF3212.1",]
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_gradient2(
    low = "#5AB4AC", mid = "white", high = "#D8B365", 
    midpoint = 0,
    guide = guide_colorbar(title = "Average log fold change", title.position = "top")
  ) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "log10_count"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, 
                       guide = guide_colorbar(title = "Average log10(count)", title.position = "top")) +
  facet_wrap(~ interaction(metric, ifelse(metric == "log10_count", as.character(hiv_status), as.character(type))),
             scales = "free", ncol = 4,
             labeller = labeller(
               `type` = facet_labels,
               `hiv_status` = label_value
             )) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./pging_viru", virulence2, ".pdf"), width = 10, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.1, .5))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./pging_viru", virulence2 , ".pdf", sep=""))
}
```
# 5. P. ging tree
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
setwd("/home/suzanne/rna_dohmain/11-perio/05-virulence/")
homd <- read.csv("../02-pgap/gene_annots.txt", header=T, sep="\t", quote="")
homd$genome <- gsub(x = homd$genome, pattern = "\\_.*", replacement = "")
load("../03-global-diff/deseq_results-HIvHUU.RData")
tree <- read.tree("~/rna_dohmain/11-perio/06-phylogenies/t_denticola/align/tdent.core.treefile")
tree$node.label <- NULL
tree.root <- midpoint_root(tree)
tip_labels <- gsub("'","",rev(tree.root$tip.label))
gsub("'","",rev(tree$tip.label))
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
  xlim(0, .095)
tip_labels <- rev(get_taxa_name(tree_plot))
# make log10 for RNA
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
merge.count <- subcount %>% group_by(genome) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by species and calculate group sums
group_sums <- merge.count %>%
  group_by(genome) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# long format
collapsed_long <- group_sums %>%
  pivot_longer(cols = -genome, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# genomes to get 
pgings_seqs <- unique(sort(homd[homd$species  == "Treponema_denticola",]$genome))

collapsed_long2 <- collapsed_long[collapsed_long$genome %in% pgings_seqs,]
rna_abund <- left_join(collapsed_long2, submap, by = c("sample" = "sample_id"))

summary_log10_counts <- rna_abund %>%
  group_by(hiv_status, genome) %>% 
  summarise(value = mean(log10_count, na.rm = TRUE))

#loop
virulence_genes_pg <- c("prtP", "prcB", "prcA", "oppA", "rpoC")
for (virulence2 in virulence_genes_pg) {

# add in annotations
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
  filter(gene ==virulence2 & species == "Treponema_denticola") %>%
  group_by(genome, species) %>%
  summarise(
    avg_baseMean_factor = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange_factor = mean(log2FoldChange, na.rm = TRUE)
  )

# get avg_log2FoldChange for virus
sig_average_log2_200days <- average_pg %>%
  select(genome, avg_log2FoldChange_factor, species) %>%
  pivot_longer(cols = "avg_log2FoldChange_factor", 
               names_to = "metric", 
               values_to = "value")
sig_average_baseMean_200days <- average_pg %>%
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
  xlim(0, .055)
tip_labels <- rev(get_taxa_name(tree_plot))

#base mean
summary_log10_counts$metric <- "log10_count"
combined_metrics<- bind_rows(
  sig_average_baseMean_200days %>% mutate(type = "baseMean"),
  sig_average_log2_200days %>% mutate(type = "log2FoldChange"),
  summary_log10_counts %>% mutate(type = "log10(count)")
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
  "baseMean" = "baseMean",
  "log2FoldChange" = "log2FoldChange",
  "log10(count)" = "log10(count)"
)
combined_metrics_complete[combined_metrics_complete$genome  == "SEQF3212.1",]
heatmap_plot <- ggplot(mapping = aes(x = metric, y = factor(genome, levels = tip_labels))) +
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_baseMean_factor"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, 
                       guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "avg_log2FoldChange_factor"), aes(fill = value), color = "black") +
  scale_fill_gradient2(
    low = "#5AB4AC", mid = "white", high = "#D8B365", 
    midpoint = 0,
    guide = guide_colorbar(title = "Average log fold change", title.position = "top")
  ) +
  ggnewscale::new_scale_fill() +  # Another new scale for the next plot
  geom_tile(data = combined_metrics_complete %>% filter(metric == "log10_count"), aes(fill = value), color = "black") +
  scale_fill_distiller(type = "seq", palette = "Purples", direction = 1, 
                       guide = guide_colorbar(title = "Average log10(count)", title.position = "top")) +
  facet_wrap(~ interaction(metric, ifelse(metric == "log10_count", as.character(hiv_status), as.character(type))),
             scales = "free", ncol = 4,
             labeller = labeller(
               `type` = facet_labels,
               `hiv_status` = label_value
             )) +
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  theme(
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(-1,'lines'),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title = element_blank() # Remove axis titles
  ) +
  labs(x = NULL, y = NULL)

pdf(paste0("./tdent_viru", virulence2, ".pdf"), width = 10, height = 8)
print(tree_plot + heatmap_plot + plot_layout(ncol = 2, widths = c(.45, .5))) 
dev.off()
system(paste0("~/.iterm2/imgcat ./tdent_viru", virulence2 , ".pdf", sep=""))
}
```
