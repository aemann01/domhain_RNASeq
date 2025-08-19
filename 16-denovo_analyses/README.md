# Denovo assembly analyses

### 1. Set up environment

```bash
cd ~/domhain_RNAseq/12-denovo_analyses
# install singularity version of prokka
singularity build prokka.sif docker://staphb/prokka:latest
singularity exec prokka.sif prokka -h # test to see if install worked
```

### 2. Annotate unique transcripts using prokka

```bash
# make headers prokka compliant
sed -i 's/_/ .*/2' ../02-denovo_assembly/all_transcripts.uniq.fa
sed -i 's/ .*//' ../02-denovo_assembly/all_transcripts.uniq.fa 
# add unique numbers to sequence identifiers
awk '/^>/ { $0 = $0 "_" ++count } 1' ../02-denovo_assembly/all_transcripts.uniq.fa > all_transcripts.fasta

# run prokka for gene annotations on all unique transcripts generated during denovo assembly
singularity exec prokka.sif prokka \
		--cpus 40 \
        --outdir prokka_out \
        --addgenes \
        --norrna \
        --notrna \
        --force \
        all_transcripts.fasta >> prokka.out 2>> prokka.err &
```

### 3. Identify any transcripts with arcABC

```bash
cd ~/domhain_RNAseq/12-denovo_analyses
# pull records for all arc genes annotated by prokka
grep -E "eC_number=3.5.3.6|eC_number=2.1.3.3|eC_number=2.7.2.2" prokka_out/PROKKA_04022025.gff > arc_genes.gff
# clean up to only include the data that we need and format for parsing
awk -F"\t" '{print $1 "\t" $4 "\t" $5 "\t" $9}' arc_genes.gff | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > arc_genes.clean.txt
# add in header line
sed -i '1i homdID\tstart\tend\tlocustag\tEC' arc_genes.clean.txt
# load environment
conda activate 2024-HIV_RNASeq
ipython
```

### 4. Identify ADS operon (modified from /10-ads_DE/operon_slice.md)

```python
# import libraries
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

arcgenes = pd.read_csv("arc_genes.clean.txt", header=0, sep="\t")
# get length of each concatenated gene
arcgenes["length"] = abs(arcgenes["end"] - arcgenes["start"])
# add gene names
arcgenes.loc[arcgenes['EC'] == "eC_number=3.5.3.6", 'EC'] = 'arcA'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.1.3.3", 'EC'] = 'arcB'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.7.2.2", 'EC'] = 'arcC'

# first need to check that each transcript ID has at lest the three genes required for the ADS pathway
ads_genes = {'arcA', 'arcB', 'arcC'}

def check_genes(homdID): # using homdID because it was the indicator in the last script
    return ads_genes.issubset(set(homdID['EC']))

genespres = arcgenes.groupby('homdID').filter(check_genes)

# check to see if the leptos are still present
genespres["homdID"].str.contains('SEQF2444').any() # should return true

# slice out locus tag number and add to dataframe
genespres["locusnum"] = genespres.locustag.str.split('_', expand=True)[1].astype("int")

# define the allowable deviation from synteny
allowdev = 2

def find_consecutive_groups(sub_df):
    sub_df = sub_df.sort_values('locusnum').reset_index(drop=True)
    sub_df['diff'] = sub_df['locusnum'].diff().fillna(0)
    sub_df['group'] = (sub_df['diff'] > allowdev).cumsum()
    consecutive_groups = sub_df.groupby('group')['locustag'].apply(list).tolist()
    return [group for group in consecutive_groups if len(group) >= 3]

# group by genome ID
result = genespres.groupby('homdID').apply(find_consecutive_groups).tolist()

# flatten result
flatres = [group for sublist in result for group in sublist]

# make a new dataframe and merge with previous data
result_df = pd.DataFrame({'syntenous_loci': flatres})

# explode lists into indivudal rows
exploded_df = result_df.explode('syntenous_loci').reset_index(drop=True)

# merge with previous data
filtarcgenes = pd.merge(arcgenes, exploded_df, left_on="locustag", right_on="syntenous_loci")

# save this to file for troubleshooting purposes
filtarcgenes.to_csv('filtered_arc_genes.txt', index=False, sep="\t")

# next, we know that stutter annotations (or fragmented annotations) are an issue with prokka, need to identify these genes so we can fix downstream after defining the genes that are in full syntenous operons

# to give it some wiggle room (i.e., length distribution) first define a function to calculate the interquartile range
# Function to calculate IQR (x2)
def iqr(series):
    iqr = series.quantile(0.75) - series.quantile(0.25)
    return iqr * 2

genlen_stats = filtarcgenes.groupby('EC').agg(
    median_length=('length', 'median'),
    IQR_length=('length', iqr)
).reset_index()

genlen_stats.columns = ['EC', 'med_len', 'IQR_len']

# find genes that may be stutter annotated by prokka and flag for further investigation
dups = filtarcgenes[filtarcgenes.duplicated(subset=['homdID', 'EC'], keep=False)]
# calculate the sum of the length for each duplicated gene
sumlen_dups = dups.groupby(['homdID', 'EC'])['length'].sum().reset_index()
sumlen_dups.columns = ['homdID', 'EC', 'sum_length']
# merge with the median length dataframe
comparison = pd.merge(sumlen_dups, genlen_stats, on='EC')

# what I want to know if the genes in this list really look like a duplicated gene (i.e., double the length of the expected median) or if they are less than that (and are therefore probably a stutter annotation)
# calculate lower bound for the gene to be a duplicate
comparison['lower_bound'] = 2 * comparison['med_len'] - comparison['IQR_len']

# is summed length of suspected stutter genes is less than or equal to a suspected duplicate length?
comparison['within_range'] = comparison.apply(
    lambda row: row['sum_length'] <= row['lower_bound'], axis=1
)

# how many suspected true duplicate genes (False) compared to suspected stutter annotated genes (True) do we have?
comparison['within_range'].value_counts()
# using the true values query our filtarcgenes data frame to print out a list of locus tags that we will need to sum together before analysis
# get list of true values
comptrue = comparison[comparison['within_range']]
stutter_genes = pd.merge(comptrue, filtarcgenes, on=['homdID', 'EC'])
# write to file
stutter_genes.to_csv('stutter_genes_to_fix.txt', index=False, sep="\t")

# now group by homdID and get min/max coordinates for the operon
filt_df = filtarcgenes.groupby("homdID")["length"].sum().reset_index().merge(filtarcgenes.groupby("homdID")["start"].min().reset_index()).merge(filtarcgenes.groupby("homdID")["end"].max().reset_index())
# get actual operon length from start and end coordinates
filt_df["operon_len"] = filt_df["end"] - filt_df["start"]

# first remove extreme outliers (longer than 1.5x the median)
filt_df = filt_df[(filt_df['operon_len'] <= (filt_df['operon_len'].median() * 1.5))]

# calculate summary statistics across operon length column after removing outliers
summary_stats = {
    'mean': filt_df['operon_len'].mean(),
    'median': filt_df['operon_len'].median(),
    'std': filt_df['operon_len'].std(),
    'min': filt_df['operon_len'].min(),
    'max': filt_df['operon_len'].max(),
    '25%': filt_df['operon_len'].quantile(0.25),
    '50%': filt_df['operon_len'].quantile(0.50),  # same as median
    '75%': filt_df['operon_len'].quantile(0.75),
    'count': filt_df['operon_len'].count(),
}
print(summary_stats)

# {'mean': 3311.656314699793, 'median': 3414.0, 'std': 330.4715132268441, 'min': 2071, 'max': 4796, '25%': 3294.0, '50%': 3414.0, '75%': 3437.0, 'count': 483}
# since data seems to be generally normally distributed (mean ~= median), will use IQR to define length distribution threshold 
iqr = (filt_df["operon_len"].quantile(0.75) - filt_df["operon_len"].quantile(0.25)) * 2

# finally, filter the dataframe based on the 25th and 75th percentiles to remove outliers (plus minus iqr)
q25 = filt_df["operon_len"].quantile(0.25) - iqr
q75 = filt_df["operon_len"].quantile(0.75) + iqr

lenfilt_df = filt_df[(filt_df['operon_len'] >= q25) & (filt_df['operon_len'] <= q75)]

# check that you didn't remove leptos
# lenfilt_df['homdID'].str.contains('SEQF2444').any()

# finally, save to file
lenfilt_df.to_csv('filtered_ads_operon.txt', index=False, sep="\t")
# so using these parameters, we are left with 407 transcripts
```

### 5. Filter full ads operons from your denovo assembly fna and check for any odd operons

```bash
# slice these out from each transcript -- should only be left with arcABC (and any intermediate genes)
awk -F"\t" '{print $1 "\t" $3 "\t" $4}' filtered_ads_operon.txt > filtcoord
seqtk subseq prokka_out/PROKKA_04022025.fna filtcoord > arcABC_operons.fna
# dereplicate full length (if there are any)
vsearch --derep_fulllength arcABC_operons.fna --output arcABC_operons.uniq.fna # 371 unique sequences
# align 
 mafft  --adjustdirectionaccurately \
   arcABC_operons.uniq.fna > arcABC_operons.align.fna
# trim alignment 
trimal -in arcABC_operons.align.fna \
	-out arcABC_operons.trim.fna \
	-htmlout arcABC_operons.trim.html \
	-gt 0.5 \
	-resoverlap 0.5 \
    -seqoverlap 50
# build into ML tree to check for oddities
rm *tre
raxmlHPC-PTHREADS-SSE3 -T 40 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31514 \
	-f a \
	-N 100 \
	-x 02938 \
	-n tre -s arcABC_operons.trim.fna
# these look pretty good, do they cluster by taxonomy?
```

### 6. Revert alignment to get a cleaned up reference database

```python
from Bio.Seq import Seq
from Bio import SeqIO
with open("arcABC_operons.clean.fna", "w") as o:
    for record in SeqIO.parse("arcABC_operons.align.fna", "fasta"):
        record.seq = record.seq.ungap("-")
        SeqIO.write(record, o, "fasta")
```

### 6. Get taxonomy for each operon to check the tree 

```bash
# first need to download kraken2 database 
# uncomment out below to only download bacterial database
# kraken2-build --download-library bacteria \
# 			  --db kraken2-bac \
# 			  --threads 25 \
# 			  1>krakenbuild.out 2>krakenbuild.err & # only need to run this once
kraken2-build --standard --db kraken_nt # only need to run this once
# get taxonomy for each operon with kraken2 -- ran this on pickles since nt database already downloaded
kraken2 --db kraken_nt \
		--use-names \
		--output arcABC_operons.kraken.out \
		--unclassified arcABC_operons.unclassified.kraken.out \
		--confidence 0.01 \
		--threads 40 \
		arcABC_operons.clean.fna
# 370 of the sequences were classified, one unclassified -- what does blast think this is??
# create annotations file for tree
awk -F"\t" '{print $2 "\t" $3}' arcABC_operons.kraken.out | sed 's/ (taxid.*//' | sed 's/ /_/g' | sed '1i seqid\ttaxonomy' | sed 's/:.*\t/\t/' > annotations.txt
```

### 7. Rebuild GFF file from contiguous ADS operon sequences, generate star reference database, remap filtered reads to new ADS database

```bash
singularity exec prokka.sif prokka \
		--cpus 40 \
        --outdir prokka_out_ADS \
        --addgenes \
        --norrna \
        --notrna \
        --force \
        arcABC_operons.clean.fna >> prokka.out 2>> prokka.err &
# convert gff to gtf for star
conda activate 2024-HIV_RNASeq
# conda install bioconda::gffread=0.12.7
gffread -T prokka_out_ADS/PROKKA_04112025.gff -o arcABC_assemblies.gtf
# conda install bioconda::star=2.7.11b
# generate reference database for STAR
STAR --runMode genomeGenerate \
   --genomeFastaFiles arcABC_operons.clean.fna  \
   --runThreadN 10 \
   --limitGenomeGenerateRAM 228697571594 \
   --sjdbGTFfile arcABC_assemblies.gtf \
   --genomeChrBinNbits 10 \
   --limitSjdbInsertNsj 3204994 \
   --genomeSAindexNbases 10\
   --sjdbGTFfeatureExon CDS

# now need to remap all of our samples to this smaller database
cd ~/domhain_RNAseq/01-processing/filtered
# run star mapping on cleaned up ADS operons
ls *1.fastq.gz | sed 's/.1.fastq.gz//' | while read line; do
	STAR --runThreadN 40 \
   --genomeDir ~/domhain_RNAseq/12-denovo_analyses/GenomeDir \
   --readFilesIn <(gunzip -c $line.1.fastq.gz) <(gunzip -c $line.2.fastq.gz) \
   --outFileNamePrefix ~/domhain_RNAseq/12-denovo_analyses/$line.denovo \
   --outSAMtype BAM Unsorted \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold ;
done &
```

### 8. Run feature counts on resulting bam file

```bash
cd /home/allie/domhain_RNAseq/12-denovo_analyses
mkdir featurecounts

# run feature counts
ls *denovoAligned.out.bam | sed 's/.denovoAligned.out.bam//' | while read line; do
    featureCounts -f -p -C -B -M -a arcABC_assemblies.gtf \
    -o featurecounts/$line.out \
    -T 40 $line.denovoAligned.out.bam \
    -t CDS \
    -g transcript_id; done

# combine output
cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt

# clean up sample names
sed 's/.denovoAligned.out.bam//g' read_counts.txt | sed 's/_S..\t/\t/'g | sed 's/_S.\t/\t/g' > temp
mv temp read_counts.txt
```

### 9. Load denovo assembly data for deseq analysis

```R
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("apeglm")
library(ggplot2, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
# Load data 
setwd("~/domhain_RNAseq/12-denovo_analyses")
metadata <- read.table("~/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id

# read in gene counts file
genecounts <- read.table("~/domhain_RNAseq/12-denovo_analyses/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# also fix sample name error in genecounts file
genecounts <- genecounts %>% 
  rename(DM00428V2PQ75 = DM00428V2PQ84)
# remove empty last column (if you run this more than once it will start removing actual samples, make sure your dim after is expected)
genecounts <- genecounts[, -ncol(genecounts)]
dim(genecounts)    
```

### 10. First comparing healthy teeth only HI vs HUU

```R
# filter genecounts and metadata to only include healthy teeth
submap <- metadata[metadata$tooth_health == "H",]
# HI vs HUU
submap <- submap[submap$hiv_status == "HI" | submap$hiv_status == "HUU",]
dim(submap)
# filter genecounts by submap samples
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
dim(subcount)
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
# out of 938 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 99, 11%
# LFC < 0 (down)     : 180, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 19, 2%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resultsNames(se_star)
# HUU is positive, HI negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# out of 938 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 141, 15%
# LFC < 0 (down)     : 204, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 19, 2%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-H-HIvHUU.txt", quote=F, sep="\t")
# save.image()

## PCA Plot
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-H-HIvHUU.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-H-HIvHUU.pdf")
# like before, no real clustering patterns
```

### 11. Tree figure

```bash
# First need to get file that links locus tag to transcript of origin to collapse our deseq results
grep "NODE" PROKKA_04112025.gff | grep -v "^>" | awk -F"\t" '{print $1 "\t" $9}' | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $2}' | grep "gene" | sed 's/ID=//' | sed 's/:.*\t/\t/' > locustag_transcript.txt
```


```R
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ggtreeExtra")
# install.packages("ggstar")
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
library(dplyr)
library(stringr)

# read in newick tree file (this was midpoint rooted in figtree first)
tree <- read.tree("RAxML_bestTree-root.tre")

# read in transcript/locus tag annotation file
annot <- read.table("prokka_out_ADS/locustag_transcript.txt", header=F)
# read back in deseq results
deres <- read.table("deseq_results_ADS-H-HIvHUU.txt", header=T)
# merge based on rownames and V2
deres <- deres %>% 
  tibble::rownames_to_column("RowNames")
# merge
mergedat <- annot %>% 
  left_join(deres, by = c("V2" = "RowNames"))
# group by transcript ID, get the average of all numerical values 
mean_sum <- mergedat %>%
  group_by(V1) %>%
  summarize(
    across(c(baseMean, log2FoldChange, lfcSE, pvalue, padj), 
           ~ mean(., na.rm = TRUE),
           .names = "mean_{.col}")
  )
# now add in taxonomic identifiers for annotating the tree
tax <- read.table("annotations.txt", header=T)
mean_sum <- mean_sum %>% 
  left_join(tax, by = c("V1" = "seqid"))
# convert to dataframe 
mean_sum <- as.data.frame(mean_sum)
# get genus only column 
mean_sum <- mean_sum %>%
  mutate(genus = str_extract(taxonomy, "^[^_ ]+"))
# make genus a factor
mean_sum$genus <- as.factor(mean_sum$genus)
# rename column
colnames(mean_sum)[colnames(mean_sum) == "V1"] <- "label"

# sanity check -- do the tip labels match up with our mean_sum file?
tips <- tree$tip.label
nodes <- mean_sum$label
setdiff(tips, nodes) # should come back as character(0)

# base tree
p <- ggtree(tree, layout="fan", open.angle=10, size=0.25)

# genus colors
gencol <- c(
    "Actinomyces" = "#F04932",
    "Bacteria" = "black",
    "Bulleidia" = "#F8E0B1",
    "Clostridium" = "#97002E",
    "Cutibacterium" = "#93FE70",
    "Kingella" = "#23C0CA",
    "Streptococcus" = "#1C9D38",
    "Treponema" = "#0C5B6F", 
    "unclassified" = "#6C6C6C")

# add genus ring
p2 <- p + 
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = genus),  # Assuming 'label' matches tip labels
    width = 0.30,  # Adjust ring width
    offset = 0.05  # Adjust distance from tree
  ) +
  scale_fill_manual(
    name = "Genus",
    values = gencol,
    na.value = "grey90"  # Color for NA values
  ) +
  guides(fill = guide_legend(
    ncol = 2,  # Adjust legend columns
    keywidth = 0.5,
    keyheight = 0.5
  ))

# add heatmap ring
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = mean_log2FoldChange),
    width = 0.5,
    offset = 0.15  # Place after genus ring
  ) +
  scale_fill_gradient2(
    name = "log2FC",
    low = "#D73027",  # low is HI 
    mid = "#f7f7f7",  # White
    high = "#4575B4",  # high is HUU 
    midpoint = 0,      # Center at zero
    na.value = "grey90",
    limits = c(-5, 5)  # Symmetrical limits
  )

# finally add barchart showing base mean
p4 <- p3 +
  new_scale_fill() +  # Create new scale for this layer
  geom_fruit(
    data = mean_sum,
    geom = geom_bar,
    mapping = aes(
      y = label, 
      x = mean_baseMean,  
      fill = genus   # color bars by genus
    ),
    stat = "identity",
    orientation = "y",
    width = 0.5,     # Width of the bar chart ring
    offset = 0.15      # Position after previous rings
  ) +
  # Option 1: Use genus colors consistent with first ring
  scale_fill_manual(
    values = gencol,
    guide = "none"    # Hide legend since we already have genus colors
  ) +
  # Add space between rings
  theme(legend.position = "right")

p_final <- p4 +
  # Add significance markers
  geom_fruit(
    data = subset(mean_sum, mean_padj < 0.05),
    geom = geom_point,
    mapping = aes(y = label),
    shape = 8,          # Star shape
    size = 0.5,         # Adjust size
    color = "black",    # Marker color
    fill = "black",     # For shapes with fill
    stroke = 0.5,       # Border thickness
    offset = 0.005       # Position outside other rings
  )

# how many up or down regulated taxa?
mean_sum %>%
  filter(mean_padj < 0.05) %>%
  summarize(
    total_sig = n(),
    upregulated = sum(mean_log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(mean_log2FoldChange < 0, na.rm = TRUE),
    .groups = 'drop'
  )

#   total_sig upregulated downregulated
# 1        49          19            30
# more upregulated taxa (mostly in Treponema but all at very low (comparatively) base mean)

pdf("denovo_tree.H-HIvHUU.pdf")
p_final
dev.off()
system("/home/allie/.iterm2/imgcat denovo_tree.H-HIvHUU.pdf")

# I want to have a companion plot of the distribution of genera up or down regulated in these categories 
mean_sum_enriched <- mean_sum %>%
  mutate(
    HIV_status = case_when(
      mean_log2FoldChange > 0 & mean_padj < 0.05 ~ "HUU_upregulated",
      mean_log2FoldChange < 0 & mean_padj < 0.05 ~ "HI_upregulated",
      TRUE ~ "Non_sig"  # Not significant
    ),
    Genus = str_extract(taxonomy, "^[^_]+")  # Extract genus (e.g., "Streptococcus" from "Streptococcus_sp.")
  ) %>%
  filter(HIV_status != "Non_sig")  # Keep only significant genes

genus_counts <- mean_sum_enriched %>%
  count(HIV_status, Genus) %>%
  group_by(HIV_status) %>%
  mutate(Percent = 100 * n / sum(n)) %>%  # Convert to percentage
  ungroup()

# Plot 
p <- ggplot(genus_counts, aes(x = HIV_status, y = Percent, fill = Genus)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  geom_text(
    aes(label = ifelse(Percent > 5, sprintf("%.1f%%", Percent), "")),  # Label only >5%
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = gencol) +  
  labs(
    title = "Taxonomic Distribution of Upregulated Genes",
    x = "HIV Status",
    y = "Percentage of Upregulated Genes",
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

pdf("genera_prop.H-HIvHUU.pdf")
p
dev.off()
system("/home/allie/.iterm2/imgcat genera_prop.H-HIvHUU.pdf")
```

### 12. Comparing healthy teeth HUU vs HEU

```R
# filter genecounts and metadata to only include healthy teeth
submap <- metadata[metadata$tooth_health == "H",]
# HI vs HUU
submap <- submap[submap$hiv_status == "HEU" | submap$hiv_status == "HUU",]
dim(submap)
# filter genecounts by submap samples
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
dim(subcount)
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
# out of 923 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 65, 7%
# LFC < 0 (down)     : 93, 10%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resultsNames(se_star)
# HUU is positive, HEU negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# out of 923 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 102, 11%
# LFC < 0 (down)     : 120, 13%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-H-HEUvHUU.txt", quote=F, sep="\t")
# save.image()

## PCA Plot
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-H-HEUvHUU.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-H-HEUvHUU.pdf")

### Tree figure

# read back in deseq results
deres <- read.table("deseq_results_ADS-H-HEUvHUU.txt", header=T)
# merge based on rownames and V2
deres <- deres %>% 
  tibble::rownames_to_column("RowNames")
# merge
mergedat <- annot %>% 
  left_join(deres, by = c("V2" = "RowNames"))
# group by transcript ID, get the average of all numerical values 
mean_sum <- mergedat %>%
  group_by(V1) %>%
  summarize(
    across(c(baseMean, log2FoldChange, lfcSE, pvalue, padj), 
           ~ mean(., na.rm = TRUE),
           .names = "mean_{.col}")
  )
# now add in taxonomic identifiers for annotating the tree
tax <- read.table("annotations.txt", header=T)
mean_sum <- mean_sum %>% 
  left_join(tax, by = c("V1" = "seqid"))
# convert to dataframe 
mean_sum <- as.data.frame(mean_sum)
# get genus only column 
mean_sum <- mean_sum %>%
  mutate(genus = str_extract(taxonomy, "^[^_ ]+"))
# make genus a factor
mean_sum$genus <- as.factor(mean_sum$genus)
# rename column
colnames(mean_sum)[colnames(mean_sum) == "V1"] <- "label"

# sanity check -- do the tip labels match up with our mean_sum file?
tips <- tree$tip.label
nodes <- mean_sum$label
setdiff(tips, nodes) # should come back as character(0)

# base tree
p <- ggtree(tree, layout="fan", open.angle=10, size=0.25)

# genus colors
gencol <- c(
    "Actinomyces" = "#F04932",
    "Bacteria" = "black",
    "Bulleidia" = "#F8E0B1",
    "Clostridium" = "#97002E",
    "Cutibacterium" = "#93FE70",
    "Kingella" = "#23C0CA",
    "Streptococcus" = "#1C9D38",
    "Treponema" = "#0C5B6F", 
    "unclassified" = "#6C6C6C")

# add genus ring
p2 <- p + 
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = genus),  # Assuming 'label' matches tip labels
    width = 0.30,  # Adjust ring width
    offset = 0.05  # Adjust distance from tree
  ) +
  scale_fill_manual(
    name = "Genus",
    values = gencol,
    na.value = "grey90"  # Color for NA values
  ) +
  guides(fill = guide_legend(
    ncol = 2,  # Adjust legend columns
    keywidth = 0.5,
    keyheight = 0.5
  ))

# add heatmap ring
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = mean_log2FoldChange),
    width = 0.5,
    offset = 0.15  # Place after genus ring
  ) +
  scale_fill_gradient2(
    name = "log2FC",
    low = "#D73027",  # low is HEU
    mid = "#f7f7f7",  # White
    high = "#4575B4",  # high is HUU color
    midpoint = 0,      # Center at zero
    na.value = "grey90",
    limits = c(-5, 5)  # Symmetrical limits
  )

# finally add barchart showing base mean
p4 <- p3 +
  new_scale_fill() +  # Create new scale for this layer
  geom_fruit(
    data = mean_sum,
    geom = geom_bar,
    mapping = aes(
      y = label, 
      x = mean_baseMean,  
      fill = genus   # color bars by genus
    ),
    stat = "identity",
    orientation = "y",
    width = 0.5,     # Width of the bar chart ring
    offset = 0.15      # Position after previous rings
  ) +
  # Option 1: Use genus colors consistent with first ring
  scale_fill_manual(
    values = gencol,
    guide = "none"    # Hide legend since we already have genus colors
  ) +
  # Add space between rings
  theme(legend.position = "right")

p_final <- p4 +
  # Add significance markers
  geom_fruit(
    data = subset(mean_sum, mean_padj < 0.05),
    geom = geom_point,
    mapping = aes(y = label),
    shape = 8,          # Star shape
    size = 0.5,         # Adjust size
    color = "black",    # Marker color
    fill = "black",     # For shapes with fill
    stroke = 0.5,       # Border thickness
    offset = 0.005       # Position outside other rings
  )

# how many up or down regulated taxa?
mean_sum %>%
  filter(mean_padj < 0.05) %>%
  summarize(
    total_sig = n(),
    upregulated = sum(mean_log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(mean_log2FoldChange < 0, na.rm = TRUE),
    .groups = 'drop'
  )

#   total_sig upregulated downregulated
# 1        19           9            10

pdf("denovo_tree.H-HEUvHUU.pdf")
p_final
dev.off()
system("/home/allie/.iterm2/imgcat denovo_tree.H-HEUvHUU.pdf")

# taxa stacked plot
mean_sum_enriched <- mean_sum %>%
  mutate(
    HIV_status = case_when(
      mean_log2FoldChange > 0 & mean_padj < 0.05 ~ "HUU_upregulated",
      mean_log2FoldChange < 0 & mean_padj < 0.05 ~ "HEU_upregulated",
      TRUE ~ "Non_sig"  # Not significant
    ),
    Genus = str_extract(taxonomy, "^[^_]+")  # Extract genus (e.g., "Streptococcus" from "Streptococcus_sp.")
  ) %>%
  filter(HIV_status != "Non_sig")  # Keep only significant genes

genus_counts <- mean_sum_enriched %>%
  count(HIV_status, Genus) %>%
  group_by(HIV_status) %>%
  mutate(Percent = 100 * n / sum(n)) %>%  # Convert to percentage
  ungroup()

# Plot 
p <- ggplot(genus_counts, aes(x = HIV_status, y = Percent, fill = Genus)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  geom_text(
    aes(label = ifelse(Percent > 5, sprintf("%.1f%%", Percent), "")),  # Label only >5%
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = gencol) +  
  labs(
    title = "Taxonomic Distribution of Upregulated Genes",
    x = "HIV Status",
    y = "Percentage of Upregulated Genes",
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

pdf("genera_prop.H-HEUvHUU.pdf")
p
dev.off()
system("/home/allie/.iterm2/imgcat genera_prop.H-HEUvHUU.pdf")
```

### 13. Comparing diseased teeth HI vs HUU

```R
# filter genecounts and metadata to only include diseased teeth
submap <- metadata[metadata$tooth_health == "D",]
# HI vs HUU
submap <- submap[submap$hiv_status == "HI" | submap$hiv_status == "HUU",]
dim(submap)
# filter genecounts by submap samples
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
dim(subcount)
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
# out of 858 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 36, 4.2%
# LFC < 0 (down)     : 38, 4.4%
# outliers [1]       : 335, 39%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resultsNames(se_star)
# HUU is positive, HI negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HI", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# out of 858 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 49, 5.7%
# LFC < 0 (down)     : 44, 5.1%
# outliers [1]       : 335, 39%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-D-HIvHUU.txt", quote=F, sep="\t")
# save.image()

## PCA Plot
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-D-HIvHUU.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-D-HIvHUU.pdf")

### Tree figure

# read back in deseq results
deres <- read.table("deseq_results_ADS-D-HIvHUU.txt", header=T)
# merge based on rownames and V2
deres <- deres %>% 
  tibble::rownames_to_column("RowNames")
# merge
mergedat <- annot %>% 
  left_join(deres, by = c("V2" = "RowNames"))
# group by transcript ID, get the average of all numerical values 
mean_sum <- mergedat %>%
  group_by(V1) %>%
  summarize(
    across(c(baseMean, log2FoldChange, lfcSE, pvalue, padj), 
           ~ mean(., na.rm = TRUE),
           .names = "mean_{.col}")
  )
# now add in taxonomic identifiers for annotating the tree
tax <- read.table("annotations.txt", header=T)
mean_sum <- mean_sum %>% 
  left_join(tax, by = c("V1" = "seqid"))
# convert to dataframe 
mean_sum <- as.data.frame(mean_sum)
# get genus only column 
mean_sum <- mean_sum %>%
  mutate(genus = str_extract(taxonomy, "^[^_ ]+"))
# make genus a factor
mean_sum$genus <- as.factor(mean_sum$genus)
# rename column
colnames(mean_sum)[colnames(mean_sum) == "V1"] <- "label"

# sanity check -- do the tip labels match up with our mean_sum file?
tips <- tree$tip.label
nodes <- mean_sum$label
setdiff(tips, nodes) # should come back as character(0)

# base tree
p <- ggtree(tree, layout="fan", open.angle=10, size=0.25)

# genus colors
gencol <- c(
    "Actinomyces" = "#F04932",
    "Bacteria" = "black",
    "Bulleidia" = "#F8E0B1",
    "Clostridium" = "#97002E",
    "Cutibacterium" = "#93FE70",
    "Kingella" = "#23C0CA",
    "Streptococcus" = "#1C9D38",
    "Treponema" = "#0C5B6F", 
    "unclassified" = "#6C6C6C")

# add genus ring
p2 <- p + 
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = genus),  # Assuming 'label' matches tip labels
    width = 0.30,  # Adjust ring width
    offset = 0.05  # Adjust distance from tree
  ) +
  scale_fill_manual(
    name = "Genus",
    values = gencol,
    na.value = "grey90"  # Color for NA values
  ) +
  guides(fill = guide_legend(
    ncol = 2,  # Adjust legend columns
    keywidth = 0.5,
    keyheight = 0.5
  ))

# add heatmap ring
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = mean_log2FoldChange),
    width = 0.5,
    offset = 0.15  # Place after genus ring
  ) +
  scale_fill_gradient2(
    name = "log2FC",
    low = "#D73027",  # low is HI
    mid = "#f7f7f7",  # White
    high = "#4575B4",  # high is HUU color
    midpoint = 0,      # Center at zero
    na.value = "grey90",
    limits = c(-5, 5)  # Symmetrical limits
  )

# finally add barchart showing base mean
p4 <- p3 +
  new_scale_fill() +  # Create new scale for this layer
  geom_fruit(
    data = mean_sum,
    geom = geom_bar,
    mapping = aes(
      y = label, 
      x = mean_baseMean,  
      fill = genus   # color bars by genus
    ),
    stat = "identity",
    orientation = "y",
    width = 0.5,     # Width of the bar chart ring
    offset = 0.15      # Position after previous rings
  ) +
  # Option 1: Use genus colors consistent with first ring
  scale_fill_manual(
    values = gencol,
    guide = "none"    # Hide legend since we already have genus colors
  ) +
  # Add space between rings
  theme(legend.position = "right")

p_final <- p4 +
  # Add significance markers
  geom_fruit(
    data = subset(mean_sum, mean_padj < 0.05),
    geom = geom_point,
    mapping = aes(y = label),
    shape = 8,          # Star shape
    size = 0.5,         # Adjust size
    color = "black",    # Marker color
    fill = "black",     # For shapes with fill
    stroke = 0.5,       # Border thickness
    offset = 0.005       # Position outside other rings
  )

# how many up or down regulated taxa?
mean_sum %>%
  filter(mean_padj < 0.05) %>%
  summarize(
    total_sig = n(),
    upregulated = sum(mean_log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(mean_log2FoldChange < 0, na.rm = TRUE),
    .groups = 'drop'
  )

#   total_sig upregulated downregulated
# 1        26          14            12

pdf("denovo_tree.D-HIvHUU.pdf")
p_final
dev.off()
system("/home/allie/.iterm2/imgcat denovo_tree.D-HIvHUU.pdf")

# taxa stacked plot
mean_sum_enriched <- mean_sum %>%
  mutate(
    HIV_status = case_when(
      mean_log2FoldChange > 0 & mean_padj < 0.05 ~ "HUU_upregulated",
      mean_log2FoldChange < 0 & mean_padj < 0.05 ~ "HI_upregulated",
      TRUE ~ "Non_sig"  # Not significant
    ),
    Genus = str_extract(taxonomy, "^[^_]+")  # Extract genus (e.g., "Streptococcus" from "Streptococcus_sp.")
  ) %>%
  filter(HIV_status != "Non_sig")  # Keep only significant genes

genus_counts <- mean_sum_enriched %>%
  count(HIV_status, Genus) %>%
  group_by(HIV_status) %>%
  mutate(Percent = 100 * n / sum(n)) %>%  # Convert to percentage
  ungroup()

# Plot 
p <- ggplot(genus_counts, aes(x = HIV_status, y = Percent, fill = Genus)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  geom_text(
    aes(label = ifelse(Percent > 5, sprintf("%.1f%%", Percent), "")),  # Label only >5%
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = gencol) +  
  labs(
    title = "Taxonomic Distribution of Upregulated Genes",
    x = "HIV Status",
    y = "Percentage of Upregulated Genes",
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

pdf("genera_prop.D-HIvHUU.pdf")
p
dev.off()
system("/home/allie/.iterm2/imgcat genera_prop.D-HIvHUU.pdf")
```

### 14. Comparing diseased teeth HEU vs HUU

```R
# filter genecounts and metadata to only include diseased teeth
submap <- metadata[metadata$tooth_health == "D",]
# HEU vs HUU
submap <- submap[submap$hiv_status == "HEU" | submap$hiv_status == "HUU",]
dim(submap)
# filter genecounts by submap samples
subcount <- genecounts[, colnames(genecounts) %in% row.names(submap)]
dim(subcount)
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
# out of 856 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 125, 15%
# LFC < 0 (down)     : 66, 7.7%
# outliers [1]       : 335, 39%
# low counts [2]     : 33, 3.9%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resultsNames(se_star)
# HUU is positive, HEU negative
resLFC <- lfcShrink(se_star, coef="hiv_status_HUU_vs_HEU", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# out of 856 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 160, 19%
# LFC < 0 (down)     : 78, 9.1%
# outliers [1]       : 335, 39%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write results to file
write.table(resLFC, file="deseq_results_ADS-D-HEUvHUU.txt", quote=F, sep="\t")
# save.image()

## PCA Plot
# transform for visualizations
vld <- varianceStabilizingTransformation(se_star, fitType="local")
pdf("pca_pdvpf_ADS-D-HEUvHUU.pdf")
plotPCA(vld, intgroup=c("hiv_status")) + theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat pca_pdvpf_ADS-D-HEUvHUU.pdf")

### Tree figure

# read back in deseq results
deres <- read.table("deseq_results_ADS-D-HEUvHUU.txt", header=T)
# merge based on rownames and V2
deres <- deres %>% 
  tibble::rownames_to_column("RowNames")
# merge
mergedat <- annot %>% 
  left_join(deres, by = c("V2" = "RowNames"))
# group by transcript ID, get the average of all numerical values 
mean_sum <- mergedat %>%
  group_by(V1) %>%
  summarize(
    across(c(baseMean, log2FoldChange, lfcSE, pvalue, padj), 
           ~ mean(., na.rm = TRUE),
           .names = "mean_{.col}")
  )
# now add in taxonomic identifiers for annotating the tree
tax <- read.table("annotations.txt", header=T)
mean_sum <- mean_sum %>% 
  left_join(tax, by = c("V1" = "seqid"))
# convert to dataframe 
mean_sum <- as.data.frame(mean_sum)
# get genus only column 
mean_sum <- mean_sum %>%
  mutate(genus = str_extract(taxonomy, "^[^_ ]+"))
# make genus a factor
mean_sum$genus <- as.factor(mean_sum$genus)
# rename column
colnames(mean_sum)[colnames(mean_sum) == "V1"] <- "label"

# sanity check -- do the tip labels match up with our mean_sum file?
tips <- tree$tip.label
nodes <- mean_sum$label
setdiff(tips, nodes) # should come back as character(0)

# base tree
p <- ggtree(tree, layout="fan", open.angle=10, size=0.25)

# genus colors
gencol <- c(
    "Actinomyces" = "#F04932",
    "Bacteria" = "black",
    "Bulleidia" = "#F8E0B1",
    "Clostridium" = "#97002E",
    "Cutibacterium" = "#93FE70",
    "Kingella" = "#23C0CA",
    "Streptococcus" = "#1C9D38",
    "Treponema" = "#0C5B6F", 
    "unclassified" = "#6C6C6C")

# add genus ring
p2 <- p + 
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = genus),  # Assuming 'label' matches tip labels
    width = 0.30,  # Adjust ring width
    offset = 0.05  # Adjust distance from tree
  ) +
  scale_fill_manual(
    name = "Genus",
    values = gencol,
    na.value = "grey90"  # Color for NA values
  ) +
  guides(fill = guide_legend(
    ncol = 2,  # Adjust legend columns
    keywidth = 0.5,
    keyheight = 0.5
  ))

# add heatmap ring
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = mean_sum,
    geom = geom_tile,
    mapping = aes(y = label, fill = mean_log2FoldChange),
    width = 0.5,
    offset = 0.15  # Place after genus ring
  ) +
  scale_fill_gradient2(
    name = "log2FC",
    low = "#D73027",  # low is HI
    mid = "#f7f7f7",  # White
    high = "#4575B4",  # high is HUU color
    midpoint = 0,      # Center at zero
    na.value = "grey90",
    limits = c(-5, 5)  # Symmetrical limits
  )

# finally add barchart showing base mean
p4 <- p3 +
  new_scale_fill() +  # Create new scale for this layer
  geom_fruit(
    data = mean_sum,
    geom = geom_bar,
    mapping = aes(
      y = label, 
      x = mean_baseMean,  
      fill = genus   # color bars by genus
    ),
    stat = "identity",
    orientation = "y",
    width = 0.5,     # Width of the bar chart ring
    offset = 0.15      # Position after previous rings
  ) +
  # Option 1: Use genus colors consistent with first ring
  scale_fill_manual(
    values = gencol,
    guide = "none"    # Hide legend since we already have genus colors
  ) +
  # Add space between rings
  theme(legend.position = "right")

p_final <- p4 +
  # Add significance markers
  geom_fruit(
    data = subset(mean_sum, mean_padj < 0.05),
    geom = geom_point,
    mapping = aes(y = label),
    shape = 8,          # Star shape
    size = 0.5,         # Adjust size
    color = "black",    # Marker color
    fill = "black",     # For shapes with fill
    stroke = 0.5,       # Border thickness
    offset = 0.005       # Position outside other rings
  )

# how many up or down regulated taxa?
mean_sum %>%
  filter(mean_padj < 0.05) %>%
  summarize(
    total_sig = n(),
    upregulated = sum(mean_log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(mean_log2FoldChange < 0, na.rm = TRUE),
    .groups = 'drop'
  )

#   total_sig upregulated downregulated
# 1        68          44            24

pdf("denovo_tree.D-HEUvHUU.pdf")
p_final
dev.off()
system("/home/allie/.iterm2/imgcat denovo_tree.D-HEUvHUU.pdf")

# taxa stacked plot
mean_sum_enriched <- mean_sum %>%
  mutate(
    HIV_status = case_when(
      mean_log2FoldChange > 0 & mean_padj < 0.05 ~ "HUU_upregulated",
      mean_log2FoldChange < 0 & mean_padj < 0.05 ~ "HEU_upregulated",
      TRUE ~ "Non_sig"  # Not significant
    ),
    Genus = str_extract(taxonomy, "^[^_]+")  # Extract genus (e.g., "Streptococcus" from "Streptococcus_sp.")
  ) %>%
  filter(HIV_status != "Non_sig")  # Keep only significant genes

genus_counts <- mean_sum_enriched %>%
  count(HIV_status, Genus) %>%
  group_by(HIV_status) %>%
  mutate(Percent = 100 * n / sum(n)) %>%  # Convert to percentage
  ungroup()

# Plot 
p <- ggplot(genus_counts, aes(x = HIV_status, y = Percent, fill = Genus)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  geom_text(
    aes(label = ifelse(Percent > 5, sprintf("%.1f%%", Percent), "")),  # Label only >5%
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = gencol) +  
  labs(
    title = "Taxonomic Distribution of Upregulated Genes",
    x = "HIV Status",
    y = "Percentage of Upregulated Genes",
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

pdf("genera_prop.D-HEUvHUU.pdf")
p
dev.off()
system("/home/allie/.iterm2/imgcat genera_prop.D-HEUvHUU.pdf")
```