# 1. Heatmap
```R
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
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
#pcoa diversity
vld <- varianceStabilizingTransformation(se_star)
pdf("pca.hiv.red.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pca.hiv.red.pdf")

#make dendogram
#get virulence factors ids
homd <- read.table("red_annots.txt", header=T, sep="\t", quote="")
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
homd <- read.table("red_annots.txt", header=T, sep="\t", quote="")
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
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
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
save_pheatmap_pdf(x, "heatmap.hiv.red.pdf")
system("~/.iterm2/imgcat ./heatmap.hiv.red.pdf")
```
# 2. Donut Plot
Reads belonging to each group
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
genecounts <- read.table("./red_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
# genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
#get annots for genes
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
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
row <- rows2[, c("seq", "name")]
#make donut plots
gene_sum <- combined_sub %>% 
  group_by(hiv_status, contig, species) %>% 
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

setwd("~/rna_dohmain/11-perio/06-red-complex")
HI <- read.csv("deseq_results_red-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("deseq_results_red-HEUvHUU.txt",  sep = "\t")
HI$gene_tag <- row.names(HI)
HEU$gene_tag <- row.names(HEU)
#HI
HI$hiv_status <- "HI"
#get annots for HI
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
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
homd <- read.table("./red_annots.txt", header=T, sep="\t", quote="")
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
# 4. Tree 
Make core gene tree
```sh
mkdir tree && cd tree
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep -v cangingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_cds_from_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query
rm query
cat *query > query
wget -i query 

#get core genes
python3 core_genes.py > core_genes
zcat *fna.gz > all_genomes.fna
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' all_genomes.fna > temp
mv temp all_genomes.fna
cat core_genes | while read line; do echo "$line: $(grep -cw "$line" all_genomes.fna)"; done # make sure only 53 are there

# get genes for alignment
cat core_genes | while read line; do grep -w $line all_genomes.fna | awk '{print $1}' | sed 's/>//' > $line.ids; done

ls *ids | sed 's/.ids//' | while read line; do seqtk subseq all_genomes.fna $line.ids > $line.fa; done

sed -i 's/lcl|//' *fa
sed -i 's/].*//' *fa
sed -i 's/ \[/_/' *fa

# align
ls *fa | sed 's/.fa//'| while read line; do mafft --thread -1 $line.fa > $line.align.fa; done
cat *align.fa > core_genome.align.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' core_genome.align.fa > temp
mv temp core_genome.align.fa
python3 combine_core.py
grep ">" combined_core.align.fna -c

# make trees
fasttree -nt combined_core.align.fna > combined_core.tre
raxmlHPC-PTHREADS-SSE3 -T 190 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s combined_core.align.fna

# make annotation file
cat rpoC.ids | sed 's/\..*//' | sed 's/lcl|//' | sort | uniq | while read line; do grep -m 1 $line ../red_annots.txt | awk -F "\t" '{print $1, $6}' | sed 's/ /\t/' ; done | sed '1s/^/genome\tspecies\n/'> combined_core.annots.txt
```
Examin expression of genes in R
```R
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
library("cowplot")

#load data
setwd("/home/suzanne/rna_dohmain/11-perio/06-red-complex/tree")
load("../deseq_results_red-HIvHUU.RData")
tree <- read.tree("./RAxML_bestTree.tre")
tree.root <- midpoint_root(tree)

# add in annotations
homd <- read.table("../red_annots.txt", header=T, sep="\t", quote="")
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

# get tip order
tip_labels <- gsub("'","",rev(tree.root$tip.label))

sig_average <- res_ord %>%
  group_by(contig, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

sig_average_long <- sig_average %>%
  pivot_longer(cols = c(avg_baseMean, avg_log2FoldChange), 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2 <- sig_average %>%
  select(contig, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- sig_average %>%
  select(contig, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

# combine the heatmaps 
pdf("combined_heat.all_genes.pdf", width =10)
ggplot(mapping = aes(x = metric, y = factor(contig, levels=tip_labels))) +
  geom_tile(data = sig_average_log2, aes(fill = value)) +
  scale_fill_distiller(type = "div", palette = "BrBG", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = function(x) sig_average_log2$species[match(x, sig_average_log2$contig)]) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()
system("~/.iterm2/imgcat ./combined_heat.all_genes.pdf")


# now make it for virulence factors
virulence_genes_pg <- c("kgp", "rgpB", "rgpA", "hagA", "fimA")
# Tannerella forsythia
virulence_genes_tf <- c("susB", "kly", "eno", "hagA", "fimA")
# Treponema denticola
virulence_genes_td <- c("oppA", "flaA", "flaB", "fliE", "cheX", "cheY", "hbpA", "hbpB", "troA" )

# calculate the average baseMean and log2fold
average_pg <- res_ord %>%
  filter(gene %in% virulence_genes_pg & species == "Porphyromonas_gingivalis") %>%
  group_by(contig, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )
average_tf <- res_ord %>%
  filter(gene %in% virulence_genes_tf & species == "Tannerella_forsythia") %>%
  group_by(contig, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )
average_td <- res_ord %>%
  filter(gene %in% virulence_genes_td & species == "Treponema_denticola") %>%
  group_by(contig, species) %>%
  summarise(
    avg_baseMean = mean(baseMean, na.rm = TRUE),
    avg_log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
  )

# combine
average_values <- rbind(average_pg, average_tf, average_td)

# overall
sig_average_long <- average_values %>%
  pivot_longer(cols = c(avg_baseMean, avg_log2FoldChange), 
               names_to = "metric", 
               values_to = "value")

# get avg_log2FoldChange
sig_average_log2 <- average_values %>%
  select(contig, avg_log2FoldChange, species) %>%
  pivot_longer(cols = "avg_log2FoldChange", 
               names_to = "metric", 
               values_to = "value")

# get avg_baseMean
sig_average_baseMean <- average_values %>%
  select(contig, avg_baseMean, species) %>%
  pivot_longer(cols = "avg_baseMean", 
               names_to = "metric", 
               values_to = "value")

pdf("combined_heat.viru_genes.pdf", width =6)
ggplot(mapping = aes(x = metric, y = factor(contig, levels=tip_labels))) +
  geom_tile(data = sig_average_log2, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "GnBu", guide = guide_colorbar(title = "Average log fold change", title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = sig_average_baseMean, aes(fill = value)) +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1, guide = guide_colorbar(title = "Average base mean", title.position = "top")) +
  facet_grid(.~metric, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(labels = function(x) sig_average_log2$species[match(x, sig_average_log2$contig)]) +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()
system("~/.iterm2/imgcat ./combined_heat.viru_genes.pdf")


