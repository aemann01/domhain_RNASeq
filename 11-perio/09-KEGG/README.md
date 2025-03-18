# 1. Download files
```sh
# download genbank file from HOMD and prokka2kegg files
wget https://github.com/SilentGene/Bio-py/raw/master/prokka2kegg/idmapping_KO.tab.gz
wget https://github.com/SilentGene/Bio-py/raw/master/prokka2kegg/prokka2kegg.py
chmod +x prokka2kegg.py
#make gbk file
cd ~/rna_dohmain/11-perio/05-TrEMBL
```
```py
from Bio import SeqIO

def combine_genbank_files(input_files, output_file):
    records = []

    # Read each GenBank file and add the records to the list
    for input_file in input_files:
        print(f"Processing {input_file}...")
        with open(input_file, 'r') as infile:
            for record in SeqIO.parse(infile, 'genbank'):
                records.append(record)

    # Write all records to the output GenBank file
    with open(output_file, 'w') as outfile:
        SeqIO.write(records, outfile, 'genbank')
    print(f"Combined GenBank file saved as {output_file}")

# List of GenBank files to combine
input_files = [
    'GCF_000007585.1_ASM758v1_genomic.gbff', 'GCF_018141805.1_ASM1814180v1_genomic.gbff',
    'GCF_000008185.1_ASM818v1_genomic.gbff', 'GCF_023822425.1_ASM2382242v1_genomic.gbff',
    'GCF_000010505.1_ASM1050v1_genomic.gbff', 'GCF_023822445.1_ASM2382244v1_genomic.gbff',
    'GCF_000238215.1_ASM23821v1_genomic.gbff', 'GCF_023822465.1_ASM2382246v1_genomic.gbff',
    'GCF_000270225.1_ASM27022v1_genomic.gbff', 'GCF_024181405.1_ASM2418140v1_genomic.gbff',
    'GCF_001263815.1_ASM126381v1_genomic.gbff', 'GCF_024181425.1_ASM2418142v1_genomic.gbff',
    'GCF_001274615.1_ASM127461v1_genomic.gbff', 'GCF_024181445.1_ASM2418144v1_genomic.gbff',
    'GCF_001314265.1_ASM131426v1_genomic.gbff', 'GCF_024181465.1_ASM2418146v1_genomic.gbff',
    'GCF_001444325.1_ASM144432v1_genomic.gbff', 'GCF_024181485.1_ASM2418148v1_genomic.gbff',
    'GCF_001547855.1_ASM154785v1_genomic.gbff', 'GCF_024181505.1_ASM2418150v1_genomic.gbff',
    'GCF_001547875.1_ASM154787v1_genomic.gbff', 'GCF_024181545.1_ASM2418154v1_genomic.gbff',
    'GCF_002753935.1_ASM275393v1_genomic.gbff', 'GCF_024181565.1_ASM2418156v1_genomic.gbff',
    'GCF_002753955.1_ASM275395v1_genomic.gbff', 'GCF_024181605.1_ASM2418160v1_genomic.gbff',
    'GCF_002753975.1_ASM275397v1_genomic.gbff', 'GCF_024181625.1_ASM2418162v1_genomic.gbff',
    'GCF_002754015.1_ASM275401v1_genomic.gbff', 'GCF_024181645.1_ASM2418164v1_genomic.gbff',
    'GCF_002754035.1_ASM275403v1_genomic.gbff', 'GCF_024400535.1_ASM2440053v1_genomic.gbff',
    'GCF_002754055.1_ASM275405v1_genomic.gbff', 'GCF_024400725.1_ASM2440072v1_genomic.gbff',
    'GCF_002754075.1_ASM275407v1_genomic.gbff', 'GCF_025905525.1_ASM2590552v1_genomic.gbff',
    'GCF_002754095.1_ASM275409v1_genomic.gbff', 'GCF_028335085.1_ASM2833508v1_genomic.gbff',
    'GCF_002754115.1_ASM275411v1_genomic.gbff', 'GCF_028335105.1_ASM2833510v1_genomic.gbff',
    'GCF_002754135.1_ASM275413v1_genomic.gbff', 'GCF_028335125.1_ASM2833512v1_genomic.gbff',
    'GCF_002754155.1_ASM275415v1_genomic.gbff', 'GCF_030144345.1_ASM3014434v1_genomic.gbff',
    'GCF_002892555.1_ASM289255v1_genomic.gbff', 'GCF_030252365.1_ASM3025236v1_genomic.gbff',
    'GCF_002892575.1_ASM289257v1_genomic.gbff', 'GCF_030440475.1_ASM3044047v1_genomic.gbff',
    'GCF_002892595.1_ASM289259v1_genomic.gbff', 'GCF_030440495.1_ASM3044049v1_genomic.gbff',
    'GCF_018141745.1_ASM1814174v1_genomic.gbff', 'GCF_900638305.1_57043_C01_genomic.gbff',
    'GCF_018141765.1_ASM1814176v1_genomic.gbff'
]

output_file = 'combined_genbank_file.gbk'

combine_genbank_files(input_files, output_file)
```
```sh
cp combined_genbank_file.gbk ~/rna_dohmain/11-perio/09-KEGG
cd ~/rna_dohmain/11-perio/09-KEGG
# run prokka2kegg on genbank file
./prokka2kegg.py -i ../05-TrEMBL/combined_genbank_file.gbk -d idmapping_KO.tab.gz -o ALL_genomes_KO.txt
```
Get corresponding uniprot ids
```sh
python3.8 locus_tag_to_go.py
```
Get refseq IDs for each Uniprot ID, then gene ID from each accession
```sh
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
# unzip
gzip -d *.gz
# only want to include specific columns from gene2accession to reduce memory load
awk -F "\t" '{print $2}' locus2go.txt | sed '1d' | sort | uniq > go_terms
parallel -a go_terms -j 7 -k  "grep '{}' gene2go -w " > gene2go.sub
cat <(head -n 1 gene2go) gene2go.sub > temp
mv temp gene2go.sub
sed -i 's/ /_/g' gene2go.sub
# merge files together 
python3 merge_tables1.py
parallel -a <(awk '{print $2}' go_reads.txt | grep -v no_term) -j 50 -k "grep -wm 1 '{}' gene2go.sub"> go_info #get cords of all rpoC genes
parallel -a <(awk '{print $2}' go_reads.txt | grep -v no_term) -j 50 -k "grep -wm 1 '{}' gene2go.sub || echo 'no_term'" > go_info
paste -d "\t" <(cat go_reads.txt | grep -v no_term) go_info | grep -v no_term > read_counts_go.txt
paste -d "\t" <(head -n 1 go_reads.txt) <(head -n 1 gene2go.sub) > headers
cat headers read_counts_go.txt > temp
mv temp read_counts_go.txt
# python3 merge_tables.py 
# clean up header
sed -i 's/#//g' read_counts_go.txt
awk '{print $1}' read_counts_go.txt | while read line; do grep -wm 1 $line <(sed 's/_[12]//g' ../06-red-complex/red_annots.txt) | awk '{print $6}'; done > species
sed -i '1s/^/species\n/' species
paste -d "\t" read_counts_go.txt species > temp
mv temp read_counts_go.txt
```
# 2. Run deseq
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
Make Valcona Plot
```R
load("deseq_results_red-HIvHUU.RData")
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
pdf("volcano-HIvHUU.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HIvHUU.red.pdf")
```
# 3. HUU vs HEU
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
# [1] "number of genes with adjusted p value lower than 0.05:  172"
# out of 767 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 82, 11%
# LFC < 0 (down)     : 90, 12%
# outliers [1]       : 0, 0%
# low counts [2]     : 119, 16%
# (mean count < 2)
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# [1] "number of genes with adjusted p value lower than 0.05:  168"
# out of 767 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 135, 18%
# LFC < 0 (down)     : 86, 11%
# outliers [1]       : 0, 0%
# low counts [2]     : 104, 14%
# (mean count < 2)
write.table(resLFC, file="deseq_results_red-HEUvHUU.txt", quote=F, sep="\t")
save.image("deseq_results_red-HEUvHUU.RData")
```
```R
load("deseq_results_red-HEUvHUU.RData")
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
pdf("volcano-HEUvHUU.red.pdf", width=15, height=10)
overall_plot
dev.off()
system("~/.iterm2/imgcat ./volcano-HEUvHUU.red.pdf")
```
# 2. KEGG
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

python3 gene2map.py #output merged_output.tsv
sed -i 's/ /_/g' merged_output.tsv
