# 1. Download files
```sh
#make gbk file
cd ~/rna_dohmain/11-perio/09-GO_terms
mkdir annotations
cp ../02-pgap/annotations/*gbk ./annotations
ls ./annotations/*gbk | sed 's/.*SEQ/SEQ/' | sed 's/\..*//' | while read line; do
  sed -i 's/locus_tag="pgap/locus_tag="'$line'-gene-pgap/g' ./annotations/$line.[12345]_modified_annot.gbk
done
# sed -i 's/locus_tag="pgap/locus_tag="SEQF9999_pgap/g' SEQF9999.1_modified_annot.gbk
# sed 's/\(SEQF[0-9]*\.[0-9]*\)[^ ]*/\1 /g' SEQF9668.1_modified_annot.gbk | grep LOCUS
for file in ./annotations/SEQF*.gbk; do
    sed -i 's/\(SEQF[0-9]*\.[0-9]*\)[^ ]*/\1     /g' "$file"
done

for file in ./annotations/SEQF*.gbk; do
    sed -i 's/linear[[:space:]]*$/linear   BCT/' "$file"
done

# reformat locus header
python3 combine_gbk.py
# Get GO process
python3 locus_tag_to_go.py
sed -i 's/ -.*//' locus2go.txt
```
Get refseq IDs for each Uniprot ID, then gene ID from each accession
```sh
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
# unzip
gzip -d *.gz
# only want to include specific columns from gene2accession to reduce memory load
awk -F "\t" '{print $2}' locus2go.txt | sed '1d' | sort | uniq > go_terms
parallel -a go_terms -j 80 -k  "grep '{}' gene2go -w" > gene2go.sub
cat <(head -n 1 gene2go) gene2go.sub > temp
mv temp gene2go.sub
sed -i 's/ /_/g' gene2go.sub
# merge files together 
python3 merge_tables1.py
# parallel -a <(awk '{print $2}' go_reads.txt | grep -v no_term) -j 90 -k "grep -wm 1 '{}' gene2go.sub"> go_info 
parallel -a <(awk '{print $2}' go_reads.txt | grep -v no_term) -j 90 -k "grep -wm 1 '{}' gene2go.sub || echo 'no_term'" > go_info
paste -d "\t" <(cat go_reads.txt | grep -v no_term) go_info | grep -v no_term > read_counts_go.txt
paste -d "\t" <(head -n 1 go_reads.txt) <(head -n 1 gene2go.sub) > headers
cat headers read_counts_go.txt > temp
mv temp read_counts_go.txt
# python3 merge_tables.py 
# clean up header
sed -i 's/#//g' read_counts_go.txt
parallel -a <(awk '{print $1}' read_counts_go.txt) -j 90 -k "grep -wm 1 '{}'  <(sed 's/_[12345]//g' ../02-pgap/gene_annots.txt) | awk '{print \$6}'" > species
# awk '{print $1}' read_counts_go.txt | while read line; do grep -wm 1 $line <(sed 's/_[12345]//g' ../02-pgap/gene_annots.txt) | awk '{print $6}'; done > species
sed -i '1s/^/species\n/' species
paste -d "\t" read_counts_go.txt species > temp
mv temp read_counts_go.txt
grep -i "Treponema_denticola\|Porphyromonas_gingivalis\|Tannerella_forsythia\|Locus_Tag" read_counts_go.txt > red_counts_go.txt

```
# 2. Run deseq for go process
## HI vs HUU
```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
genecounts <- unique(genecounts)
genecounts$Geneid <- make.unique(as.character(genecounts$Geneid))
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Group.1")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  766"
# out of 10331 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 409, 4%
# LFC < 0 (down)     : 357, 3.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 7317, 71%
# (mean count < 1)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  766"
# out of 10331 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 481, 4.7%
# LFC < 0 (down)     : 457, 4.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 7317, 71%
# (mean count < 1)
write.table(resLFC, file="deseq_results_go-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_go-HIvHUU.RData")
```
Make Valcona Plot
```R
load("deseq_results_go-HIvHUU.RData")
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
# filter by locus tag 
resdf <- as.data.frame(resLFC)
resdf <- resdf %>%
  rownames_to_column(var = "original_name") %>%  # Convert row names to a column
  mutate(genus = str_extract(original_name, "(Tannerella|Porphyromonas|Treponema)"))  # Extract specific genera
row.names(resdf) <- resdf$original_name
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# get new dataframe with concatenated genus species and locus id


sigsp <- paste("x", sigloc$genus, sep="_")
sigloc$tag <- row.names(sigloc)
sigdf <- as.data.frame(cbind(sigloc$original_name, sigsp))

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
# res_ord <- res_ord %>% filter(gene != "none")
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
res_ord$GeneInfo <- row.names(res_ord)
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$GeneInfo, 10)
# negative top 10
top <- tail(sortdf$GeneInfo, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
# res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$GeneInfo %in% labgenes, res_ord$GeneInfo, ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap ,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HIvHUU.go.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.go.pdf")
```
## HUU vs HEU
```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
genecounts <- unique(genecounts)
genecounts$Geneid <- make.unique(as.character(genecounts$Geneid))
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Group.1")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HEU", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  506"
# out of 10331 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 168, 1.6%
# LFC < 0 (down)     : 338, 3.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 7205, 70%
# (mean count < 1)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  506"
# out of 10331 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 256, 2.5%
# LFC < 0 (down)     : 448, 4.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 7205, 70%
# (mean count < 1)
write.table(resLFC, file="deseq_results_go-HEUvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_go-HEUvHUU.RData")
```
```R
load("deseq_results_go-HEUvHUU.RData")
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
# filter by locus tag 
resdf <- as.data.frame(resLFC)
resdf <- resdf %>%
  rownames_to_column(var = "original_name") %>%  # Convert row names to a column
  mutate(genus = str_extract(original_name, "(Tannerella|Porphyromonas|Treponema)"))  # Extract specific genera
row.names(resdf) <- resdf$original_name
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# get new dataframe with concatenated genus species and locus id


sigsp <- paste("x", sigloc$genus, sep="_")
sigloc$tag <- row.names(sigloc)
sigdf <- as.data.frame(cbind(sigloc$original_name, sigsp))

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
# res_ord <- res_ord %>% filter(gene != "none")
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
res_ord$GeneInfo <- row.names(res_ord)
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$GeneInfo, 10)
# negative top 10
top <- tail(sortdf$GeneInfo, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
# res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$GeneInfo %in% labgenes, res_ord$GeneInfo, ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap ,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HEUvHUU.go.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.go.pdf")
```
Working cluster HEU vs HUU
```R
library(clusterProfiler)
# library(org.Hs.eg.db)
library(AnnotationDbi)
library(Orthology.eg.db)
library(biomaRt)
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
setwd("/home/suzanne/rna_dohmain/11-perio/09-GO_terms")
set.seed(2020)
metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
genecounts <- unique(genecounts)
genecounts$Geneid <- make.unique(as.character(genecounts$Geneid))
# genecounts$Geneid_unique <- paste0(genecounts$Geneid, "_", seq_along(genecounts$Geneid))
genecounts$Locus_tag <- paste0(genecounts$Locus_Tag, "_", seq_along(genecounts$Locus_Tag))
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
# genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
# gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <-  distinct(genecounts, Locus_tag, .keep_all = TRUE)
gogroup$Geneid <- gsub("gene-", "", gogroup$Geneid)

gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Geneid")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HEU", "HUU"))

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

resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# extract a named vector of all terms
# goterms <- Term(GOTERM)
# term2name <- data.frame("term"=names(goterms),"name"=goterms )
#get gene id to go term
genes2go <- read.csv("read_counts_go.txt", header=T, sep="\t")
genes2go <- unique(genes2go)
genes2go$Geneid <- make.unique(as.character(genes2go$Geneid))
genes2go_added <- genes2go %>%
  mutate(GO_ID = case_when(
    species == "Porphyromonas_gingivalis" ~ paste0(GO_ID, "_1"),
    species == "Treponema_denticola" ~ paste0(GO_ID, "_2"),
    species == "Tannerella_forsythia" ~ paste0(GO_ID, "_3"),
    TRUE ~ GO_ID 
  ))
genes2go <- genes2go_added[, c("GO_ID", "Geneid")]
names(genes2go) <- c("gs_name", "gene")
genes2go$gene <- gsub("gene-", "", genes2go$gene)
unique(sub(".*_", "", genes2go$gs_name))
#get go term to function
genes2go_added <- genes2go_added %>%
  mutate(GO_term = case_when(
    species == "Porphyromonas_gingivalis" ~ paste0(GO_term, "_Porphyromonas_gingivalis"),
    species == "Treponema_denticola" ~ paste0(GO_term, "_Treponema_denticola"),
    species == "Tannerella_forsythia" ~ paste0(GO_term, "_Tannerella_forsythia"),
    TRUE ~ GO_term 
  ))

term2name <- genes2go_added[, c("GO_ID", "GO_term")]
names(term2name) <- c("term", "name")
unique(sub(".*_", "", term2name$term))

#prepare lefsc
df_resLFC <- as.data.frame(resLFC)
filtered_resLFC <- df_resLFC %>%
  mutate(transformed_value = -log10(padj) * sign(log2FoldChange)) %>%
  arrange(desc(transformed_value))
lfc_vector <- filtered_resLFC$transformed_value
names(lfc_vector) <- row.names(filtered_resLFC)

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
head(lfc_vector)

intersect(row.names(filtered_resLFC), genes2go$gene)
# remove duplicate values
set.seed(123) # To make it reproducible
lfc_vector_with_noise <- lfc_vector + runif(length(lfc_vector), 0, 1e-6)
lfc_vector_with_noise <- lfc_vector_with_noise[order(lfc_vector_with_noise, decreasing = TRUE)]

gsea_results <- GSEA(
  geneList = lfc_vector_with_noise, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.5, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = genes2go,
  TERM2NAME =term2name
)
head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)

# This returns the 3 rows with the largest NES values
gsea_result_df %>%
  dplyr::slice_max(NES, n = 3)
# This returns the 3 rows with the largest NES values
gsea_result_df %>%
  dplyr::slice_min(NES, n = 3)

pdf("dot-HEUvHUU.pdf", width =10, height =20)
dotplot(gsea_results, showCategory=100, color= "p.adjust", split = ".sign" )+
facet_grid(.~.sign)
dev.off()
system("~/.iterm2/imgcat ./dot-HEUvHUU.pdf")

pdf("ridge-HEUvHUU.pdf", width =15, height =30)
ridgeplot(gsea_results, orderBy = "NES", showCategory=15)
dev.off()
system("~/.iterm2/imgcat ./ridge-HEUvHUU.pdf")


gsea_result_df$enrichmentScore
write.table(gsea_result_df, file="gsea_results_go-HEUvHUU.txt", quote=F, sep="\t")
save.image("gsea_results_go-HEUvHUU.RData")






pdf("test.pdf")
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GO:0005524_3",
  title = "TEST",
  color.line = "#0d76ff"
)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf", width =10, height =30)
dotplot(gsea_results, showCategory=100, color= "p.adjust" )
dev.off()
system("~/.iterm2/imgcat ./test.pdf")
gsea_result_df$NES

pdf("test.pdf")
cnetplot(gsea_results, categorySize="p.adjust", node_label = NULL)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

library(enrichplot)
gsea_results_with_sim <- pairwise_termsim(gsea_results)
pdf("test.pdf")
emapplot(gsea_results_with_sim)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf", width =20)
treeplot(gsea_results_with_sim, hclust_method = "average")
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf", width =40)
upsetplot(gsea_results)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf")
heatplot(gsea_results_with_sim, showCategory=5)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")



```
Working cluster HI vs HUU
```R
library(clusterProfiler)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# library(Orthology.eg.db)
library(biomaRt)
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
setwd("/home/suzanne/rna_dohmain/11-perio/09-KEGG")
set.seed(2020)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
genecounts <- unique(genecounts)
genecounts$Geneid <- make.unique(as.character(genecounts$Geneid))
# genecounts$Geneid_unique <- paste0(genecounts$Geneid, "_", seq_along(genecounts$Geneid))
genecounts$Locus_tag <- paste0(genecounts$Locus_Tag, "_", seq_along(genecounts$Locus_Tag))
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
# genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
# gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <-  distinct(genecounts, Locus_tag, .keep_all = TRUE)
gogroup$Geneid <- gsub("gene-", "", gogroup$Geneid)

gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Geneid")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

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

resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# extract a named vector of all terms
# goterms <- Term(GOTERM)
# term2name <- data.frame("term"=names(goterms),"name"=goterms )
#get gene id to go term
genes2go <- read.csv("read_counts_go.txt", header=T, sep="\t")
genes2go <- unique(genes2go)
genes2go$Geneid <- make.unique(as.character(genes2go$Geneid))
genes2go_added <- genes2go %>%
  mutate(GO_ID = case_when(
    species == "Porphyromonas_gingivalis" ~ paste0(GO_ID, "_1"),
    species == "Treponema_denticola" ~ paste0(GO_ID, "_2"),
    species == "Tannerella_forsythia" ~ paste0(GO_ID, "_3"),
    TRUE ~ GO_ID 
  ))
genes2go <- genes2go_added[, c("GO_ID", "Geneid")]
names(genes2go) <- c("gs_name", "gene")
genes2go$gene <- gsub("gene-", "", genes2go$gene)
unique(sub(".*_", "", genes2go$gs_name))
#get go term to function
genes2go_added <- genes2go_added %>%
  mutate(GO_term = case_when(
    species == "Porphyromonas_gingivalis" ~ paste0(GO_term, "_Porphyromonas_gingivalis"),
    species == "Treponema_denticola" ~ paste0(GO_term, "_Treponema_denticola"),
    species == "Tannerella_forsythia" ~ paste0(GO_term, "_Tannerella_forsythia"),
    TRUE ~ GO_term 
  ))

term2name <- genes2go_added[, c("GO_ID", "GO_term")]
names(term2name) <- c("term", "name")
unique(sub(".*_", "", term2name$term))

#prepare lefsc
df_resLFC <- as.data.frame(resLFC)
filtered_resLFC <- df_resLFC %>%
  mutate(transformed_value = -log10(padj) * sign(log2FoldChange)) %>%
  arrange(desc(transformed_value))
lfc_vector <- filtered_resLFC$transformed_value
names(lfc_vector) <- row.names(filtered_resLFC)

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
head(lfc_vector)

intersect(row.names(filtered_resLFC), genes2go$gene)
# remove duplicate values
set.seed(123) # To make it reproducible
lfc_vector_with_noise <- lfc_vector + runif(length(lfc_vector), 0, 1e-6)
lfc_vector_with_noise <- lfc_vector_with_noise[order(lfc_vector_with_noise, decreasing = TRUE)]

gsea_results <- GSEA(
  geneList = lfc_vector_with_noise, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.5, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = genes2go,
  TERM2NAME =term2name
)
head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)

# This returns the 3 rows with the largest NES values
gsea_result_df %>%
  dplyr::slice_max(NES, n = 3)
# This returns the 3 rows with the largest NES values
gsea_result_df %>%
  dplyr::slice_min(NES, n = 3)

pdf("dot-HIvHUU.pdf", width =10, height =20)
dotplot(gsea_results, showCategory=100, color= "p.adjust", split = ".sign" )+
facet_grid(.~.sign)
dev.off()
system("~/.iterm2/imgcat ./dot-HIvHUU.pdf")

pdf("ridge-HIvHUU.pdf", width =15, height =30)
ridgeplot(gsea_results, orderBy = "NES", showCategory=15)
dev.off()
system("~/.iterm2/imgcat ./ridge-HIvHUU.pdf")

gsea_result_df$enrichmentScore

write.table(gsea_result_df, file="gsea_results_go-HIvHUU.txt", quote=F, sep="\t")
save.image("gsea_results_go-HIvHUU.RData")












pdf("test.pdf")
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GO:0016639",
  title = "TEST",
  color.line = "#0d76ff"
)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf", width =10, height =30)
dotplot(gsea_results, showCategory=100, color= "p.adjust" )
dev.off()
system("~/.iterm2/imgcat ./test.pdf")
gsea_result_df$NES

pdf("test.pdf")
cnetplot(gsea_results, categorySize="p.adjust", node_label = NULL)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

library(enrichplot)
gsea_results_with_sim <- pairwise_termsim(gsea_results)
pdf("test.pdf")
emapplot(gsea_results_with_sim)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")

pdf("test.pdf")
treeplot(gsea_results_with_sim, hclust_method = "average")
dev.off()
system("~/.iterm2/imgcat ./test.pdf")


pdf("test.pdf")
heatplot(gsea_results_with_sim, showCategory=5)
dev.off()
system("~/.iterm2/imgcat ./test.pdf")
```
Make upset plot comparing HEU vs HI go terms
```R
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ggplot2)
library(dplyr)
# library(ComplexUpset)
keep_after_second_last_underscore <- function(description) {
  # Split by underscores
  parts <- strsplit(description, "_")[[1]]
  # Extract everything after the second-to-last underscore (everything after the second last part)
  second_last_part <- paste(parts[(length(parts) - 1):length(parts)], collapse = "_")
  return(second_last_part)
}

# Apply the function to the Description column
df$last_part <- sapply(df$Description, keep_after_last_underscore)


setwd("~/rna_dohmain/11-perio/09-KEGG")
HI <- read.csv("gsea_results_go-HIvHUU.txt",  sep = "\t")
HEU <- read.csv("gsea_results_go-HEUvHUU.txt",  sep = "\t")
HI$hiv_status <- "HI"
HEU$hiv_status <- "HEU"
intersect(row.names(HEU), row.names(HI))

resdf <- rbind(HI, HEU)
resdf_filtered <- resdf %>%
  dplyr::select(ID, Description, hiv_status, p.adjust)
resdf_filtered$species <- sapply(resdf_filtered$Description, keep_after_second_last_underscore)
head(resdf_filtered$last_part)


head(resdf_filtered)
pval = 0.05
resdf_filtered <- resdf_filtered %>%
  mutate(mark = if_else(p.adjust <= pval, "1", "0"))
resdf_filtered$mark <- as.numeric(resdf_filtered$mark)
resdf_filtered <- dplyr::select(resdf_filtered, mark, hiv_status, ID, species)

resdf_wide <- resdf_filtered %>%
  pivot_wider(names_from = hiv_status, values_from = mark, values_fill = list(mark = 0)) %>%
  arrange(ID)

#make upset plot
resdf_df <- as.data.frame(resdf_wide)
resdf_clean_df <- resdf_df %>%
  filter(!is.na(HI) & !is.na(HEU))

pdf("go-HI_HEU.upset.pdf")
upset(resdf_clean_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially GO terms", sets.x.label="Total Number of Differentially GO terms",
    queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola", "Tannerella_forsythia")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola","Tannerella_forsythia")), color = "#FBE51F", active = T),
        list(query = elements, 
             params = list("species", "Tannerella_forsythia"), color = "#1E7F7A", active = T)))
dev.off()

system("~/.iterm2/imgcat ./go-HI_HEU.upset.pdf")

pdf("go-HI_HEU.upset.pdf")
upset(resdf_clean_df, order.by="freq", sets = c("HEU", "HI"), mainbar.y.label="Number of Differentially GO terms", sets.x.label="Total Number of Differentially GO terms",
    queries = list(
        list(query = elements, 
             params = list("species", c("Porphyromonas_gingivalis","Treponema_denticola")), color = "#340043", active = T),
        list(query = elements, 
             params = list("species", c("Treponema_denticola")), color = "#FBE51F", active = T)))
dev.off()

system("~/.iterm2/imgcat ./go-HI_HEU.upset.pdf")
```
# 4. Deseq proccess go (on hillary)
## HI vs HUU
```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
genecounts <- unique(genecounts)
genecounts$Geneid <- make.unique(as.character(genecounts$Geneid))
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Group.1")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  273"

# out of 539 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 75, 14%
# LFC < 0 (down)     : 198, 37%
# outliers [1]       : 0, 0%
# low counts [2]     : 21, 3.9%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  266"
# out of 539 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 78, 14%
# LFC < 0 (down)     : 229, 42%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
write.table(resLFC, file="deseq_results_gop-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_gop-HIvHUU.RData")
```
Make Valcona Plot
```R
load("deseq_results_gop-HIvHUU.RData")
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
# filter by locus tag 
resdf <- as.data.frame(resLFC)
resdf <- resdf %>%
  rownames_to_column(var = "original_name") %>%  # Convert row names to a column
  mutate(genus = str_extract(original_name, "(Tannerella|Porphyromonas|Treponema)"))  # Extract specific genera
row.names(resdf) <- resdf$original_name
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# get new dataframe with concatenated genus species and locus id


sigsp <- paste("x", sigloc$genus, sep="_")
sigloc$tag <- row.names(sigloc)
sigdf <- as.data.frame(cbind(sigloc$original_name, sigsp))

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
# res_ord <- res_ord %>% filter(gene != "none")
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
res_ord$GeneInfo <- row.names(res_ord)
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$GeneInfo, 10)
# negative top 10
top <- tail(sortdf$GeneInfo, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
# res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$GeneInfo %in% labgenes, res_ord$GeneInfo, ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap ,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HIvHUU.gop.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.gop.pdf")
```
## HUU vs HEU
```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Group.1")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HEU", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  111"
# out of 539 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 32, 5.9%
# LFC < 0 (down)     : 79, 15%
# outliers [1]       : 0, 0%
# low counts [2]     : 167, 31%
# (mean count < 3)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  100"
# out of 539 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 69, 13%
# LFC < 0 (down)     : 76, 14%
# outliers [1]       : 0, 0%
# low counts [2]     : 53, 9.8%
# (mean count < 2)
write.table(resLFC, file="deseq_results_gop-HEUvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_gop-HEUvHUU.RData")
```
```R
load("deseq_results_gop-HEUvHUU.RData")
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
# filter by locus tag 
resdf <- as.data.frame(resLFC)
resdf <- resdf %>%
  rownames_to_column(var = "original_name") %>%  # Convert row names to a column
  mutate(genus = str_extract(original_name, "(Tannerella|Porphyromonas|Treponema)"))  # Extract specific genera
row.names(resdf) <- resdf$original_name
# set LFC and P value filters
lfc = 2
pval = 0.05

# get list of loci that are significant and have a log FC of 2+
sigloc <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# get new dataframe with concatenated genus species and locus id


sigsp <- paste("x", sigloc$genus, sep="_")
sigloc$tag <- row.names(sigloc)
sigdf <- as.data.frame(cbind(sigloc$original_name, sigsp))

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
# res_ord <- res_ord %>% filter(gene != "none")
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
res_ord$GeneInfo <- row.names(res_ord)
tempdf <- res_ord[res_ord$target_gene == TRUE,]
# sort the dataframe by logfold change
sortdf <- tempdf[order(tempdf$log2FoldChange),]
# get top 10 genes for disease
low <- head(sortdf$GeneInfo, 10)
# negative top 10
top <- tail(sortdf$GeneInfo, 10)
# concatenate
labgenes <- c(top, low)
#combine species and gene
# res_ord$gene_name <- gsub(x = res_ord$gene, pattern = "_.", replacement = "") 

#Create volcano plot
overall_plot <- EnhancedVolcano(res_ord,
  lab = ifelse(res_ord$GeneInfo %in% labgenes, res_ord$GeneInfo, ""),
  x = 'log2FoldChange',
  y = 'padj',
  FCcutoff = lfc,
  pCutoff = pval,
  colCustom = colormap ,
  title = "",
  subtitle = "",
  caption = "",
  shape = 19,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  max.overlaps = Inf,
  pointSize = (ifelse(rownames(res_ord) %in% all_genes == T, 3, 3)),
  colAlpha = (ifelse(rownames(res_ord) %in% all_genes == F, 0.5, 0.75)),
) 
pdf("volcano-HEUvHUU.gop.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.gop.pdf")
```







# 5. KEGG
Get geneids
```sh
# wget wget https://github.com/SilentGene/Bio-py/raw/master/prokka2kegg/idmapping_KO.tab.gz
# wget https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz
# gzip -d gene_refseq_uniprotkb_collab.gz
awk '{print $1}' ../06-red-complex/red_counts.txt | sed '1d' | sed 's/^\([^_]*\_[^_]*\)_.*/\1/' > red_genes

parallel -a red_genes -j 7 -k  "grep Parent= ../05-TrEMBL/combined.gff | grep -wm 1 '{}' | sed 's/;Parent.*//' | sed 's/.*cds-//'" > protein_ids
wc -l red_genes #22580 red_genes
wc -l protein_ids #22580 red_genes
parallel -a protein_ids -j 7 -k "grep -wm 1 '{}' gene_refseq_uniprotkb_collab || echo 'no_match'" > uniprotkb
wc -l uniprotkb #22580
sed -i 's/no_match/no_match\tno_match/' uniprotkb
parallel -a <(awk -F "\t" '{print $2}' uniprotkb) -j 7 -k "grep -w '{}' idmapping_KO.tab || echo 'no_kegg'" > ko_ids
sed -i 's/no_kegg/no_kegg\tno_kegg/' ko_ids
wc -l ko_ids #22580
paste -d "\t" red_genes protein_ids <(awk '{print $2}' uniprotkb) <(awk '{print $2}' ko_ids) > kegg_ids

rm kegg2map
awk '{print $4}' kegg_ids | sort | uniq | grep -v no_kegg | while read line; do curl -s "https://rest.kegg.jp/link/pathway/$line" >> kegg2map; done

wget https://rest.kegg.jp/list/pathway

python3 gene2map.py #output ko_read_counts.tsv
sed -i 's/ /_/g' ko_read_counts.tsv
awk '{print $1}' ko_read_counts.tsv | while read line; do grep -wm 1 $line <(sed 's/_[12]//g' ../06-red-complex/red_annots.txt) || echo $line "no_species" |  awk '{print $6}'; done > species
# sed -i '1s/^/species\n/' species
paste -d "\t" ko_read_counts.tsv species > temp
mv temp read_counts_ko.txt
```
# 3. Run deseq2
```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("ko_read_counts.tsv", header=T, sep="\t")
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Group.1")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HI" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HI", "HUU"))

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
# [1] "number of genes with adjusted p value lower than 0.05:  298"

# out of 767 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 124, 16%
# LFC < 0 (down)     : 174, 23%
# outliers [1]       : 0, 0%
# low counts [2]     : 90, 12%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  290"
# out of 767 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 136, 18%
# LFC < 0 (down)     : 205, 27%
# outliers [1]       : 0, 0%
# low counts [2]     : 15, 2%
# (mean count < 2)
write.table(resLFC, file="deseq_results_red-HIvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_red-HIvHUU.RData")
```
scratch
```R
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Orthology.eg.db)
library(biomaRt)
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(tidyverse, warn.conflicts = F, quietly = T)

metadata <- read.table("../../homd_map/map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.csv("read_counts_go.txt", header=T, sep="\t")
row.names(genecounts) <- genecounts$Geneid
# fix sample names gene count file
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.red", replacement = "") 
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")
# genecounts$GO_LOCUS  <- paste0(genecounts$GO_term,"_", genecounts$species)
# gogroup <- aggregate(genecounts[, 4:96], by=list(genecounts$GO_LOCUS), FUN=sum)
gogroup <-  distinct(genecounts, Locus_tag, .keep_all = TRUE)
gogroup$Locus_tag <- gsub("gene-", "", gogroup$Locus_tag)

gogroup <- gogroup %>% remove_rownames %>% column_to_rownames("Locus_tag")
# filter metadata so that we only compare H to D
submap <- metadata[metadata$hiv_status == "HEU" | metadata$hiv_status == "HUU",]
# submap <- submap[submap$sample_id %in% sample_list, ]
subcount <- gogroup[, colnames(gogroup) %in% row.names(submap)]
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
star_results$hiv_status <- factor(star_results$hiv_status, levels=c("HEU", "HUU"))

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

resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
#genest o test for 
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
genes_to_test <- rownames(resLFC[resLFC$log2FoldChange > 0.5,])

go_results <- enrichGO(gene = gene_list,
                       OrgDb = NULL,  # No OrgDb as we're using custom TERM2GENE
                       TERM2GENE = TERM2GENE,  # Custom mapping of GO terms to genes
                       pvalueCutoff = 0.05)



GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENTERZ", ont = "BP")








library(DESeq2)
Counts <- read.delim("./count_table.csv", header = TRUE, row.names = 1, sep = ",")
Counts <- Counts[which(rowSums(Counts) > 0),]
condition <- factor(c("C","C","C","C", "S","S","S","S"))
coldata <- data.frame(row.names = colnames(Counts), condition)

dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "S", "C"))
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05 & sigs$baseMean > 50,]


library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5,])
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 15))

png("out.pdf", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()
system("~/.iterm2/imgcat ./out.pdf")

test <- ("WP_004583429.1")
GO_results <- enrichGO(gene = test, OrgDb = "org.Hs.eg.db", keyType = "REFSEQ", ont = "BP")



library(AnnotationForge)
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "431947",
                       genus = "Porphyromonas",
                       species = "gingivalis")
org.Pgingivalis.eg.sqlite
GO_results <- enrichGO(gene = test, OrgDb = "org.Pgingivalis.eg.sqlite", keyType = "REFSEQ", ont = "BP")











```
