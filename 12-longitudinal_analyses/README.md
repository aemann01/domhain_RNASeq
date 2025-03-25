## Longitudinal changes in ADS activity and abundance of ADS competent bacteria on single teeth 

### 1. Load environment

```bash
cd ~/domhain_RNAseq/09-longitudinal_analyses
conda activate 2024-HIV_RNASeq
```

### 2. Load libraries

```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
```

### 3. Load data

```R
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("../07-ads_expression/arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")  
# reorder columns by metadata 
metadata <- metadata[order(colnames(genecounts)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(genecounts)==metadata$sample_id) # should return all true

# list samples of interest
twice <- c("DM00101V2PQ16","DM00101V3PQ16","DM00305V2PQ16","DM00305V3PQ16","DM00338V2PQ16","DM00338V3PQ16","DM00008V1PQ16-2","DM00008V2PQ16","DM00023V1PQ16-1","DM00023V2PQ16","DM00035V1PQ16","DM00035V2PQ16","DM00074V1PQ55","DM00074V2PQ55","DM00142V1PQ16","DM00142V2PQ16","DM00255V1PQ26","DM00255V2PQ26","DM00266V1PQ16","DM00266V2PQ16","DM00464V1PQ65","DM00464V2PQ65","DM00563V1PQ55","DM00563V2PQ55")
triple <- c("DM00044V1PQ16","DM00044V2PQ16","DM00044V3PQ16","DM00103V1PQ16","DM00103V2PQ16","DM00103V3PQ16","DM00254V1PQ55","DM00254V2PQ55","DM00254V3PQ55","DM00519V1PQ16","DM00519V2PQ16","DM00519V3PQ16")
```

### 4. Load annotations (this takes a bit)

```R
# load annotations
homd <- read.table("~/domhain_RNAseq/03-star_map/homd_db/annotations.merge.txt", header=T, sep="\t", quote="")
```

### 5. First create log2 expression of ADS genes dataframe at the genus level 

```R
# filter metadata
map.twice <- metadata[rownames(metadata) %in% twice,]
count.twice <- genecounts[,colnames(genecounts) %in% twice]
# remove any genes with a sum count of zero
count.twice <- count.twice[rowSums(count.twice) != 0,]
# filter annotations by those found in our genecounts file
ann <- homd[homd$locus_tag %in% rownames(count.twice),]
# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(count.twice), rownames(ann)))]
count.twice <- count.twice[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.twice)==rownames(ann)) # should all return true
# if all are true, merge together
count.twice <- cbind(count.twice, ann)
# get list of Genera to pull from rpoC data later
count.twice.gen <- unique(count.twice$Genus)

# collapse by genus and sum across rows
merge.count <- count.twice %>% group_by(Genus) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by Genus and calculate group sums
group_sums <- merge.count %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# get total count across all rows at the genus level 
group_sums <- group_sums %>% mutate(total=rowSums(across(where(is.numeric))))
# get total counts to use as denominator 
tot <- sum(group_sums$total)
thresh <- 0.01 * tot # less than 1% of total

# collapse low abundant groups
collapsed <- group_sums %>%
  mutate(Genus = ifelse(total < thresh, "Other", Genus)) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
# remove total column
collapsed <- subset(collapsed, select = -total)
```

### 6. Create stacked histograms for each sample colored by genus

```R
# set color palette
gencols <- c(Actinomyces = "#F66140", 
			Kingella = "#1FCAD4", 
			Leptotrichia = "#97002E", 
			Streptococcus = "#16AA48", 
			Treponema = "#F658F8", 
			Other = "#838383")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(Genus = factor(Genus, levels = c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("RNA_log10_ADS_200days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Genus)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat RNA_log10_ADS_200days.pdf")
```

### 7. Do the same with teeth sampled at three time points (~400 days apart)

```R
# filter metadata
map.triple <- metadata[rownames(metadata) %in% triple,]
count.triple <- genecounts[,colnames(genecounts) %in% triple]
# remove any genes with a sum count of zero
count.triple <- count.triple[rowSums(count.triple) != 0,]
# filter annotations by those found in our genecounts file
ann <- homd[homd$locus_tag %in% rownames(count.triple),]
# reorder
rownames(ann) <- ann$locus_tag
sortrow <- rownames(ann)[order(match(rownames(count.triple), rownames(ann)))]
count.triple <- count.triple[sortrow, , drop=FALSE]
ann <- ann[sortrow, , drop=FALSE]

# check that locus tags match between the two dataframes
table(rownames(count.triple)==rownames(ann)) # should all return true
# if all are true, merge together
count.triple <- cbind(count.triple, ann)
# get list of genera to pull from rpoC data later
count.triple.gen <- unique(count.triple$Genus)

# collapse by genus and sum across rows
merge.count <- count.triple %>% group_by(Genus) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# Group by Genus and calculate group sums
group_sums <- merge.count %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
# get total count across all rows at the genus level 
group_sums <- group_sums %>% mutate(total=rowSums(across(where(is.numeric))))
# get total counts to use as denominator 
tot <- sum(group_sums$total)
thresh <- 0.01 * tot # less than 1% of total
# collapse low abundant groups
collapsed <- group_sums %>%
  mutate(Genus = ifelse(total < thresh, "Other", Genus)) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
# remove total column
collapsed <- subset(collapsed, select = -total)
```

### 8. Create stacked histograms for each sample colored by genus

```R
# set color palette
gencols <- c(Actinomyces = "#F66140", 
			Kingella = "#1FCAD4", 
			Leptotrichia = "#97002E", 
			Streptococcus = "#16AA48", 
			Treponema = "#F658F8", 
			Other = "#838383")
# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(Genus = factor(Genus, levels = c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.triple %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(3.5, length(unique(collapsed_long$sample)) -0.5, by = 3)

pdf("RNA_log10_ADS_400days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Genus)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat RNA_log10_ADS_400days.pdf")
```

### 9. Load in rpoC ASV count data and filter out species of interest (from RNA seq data)

```R
# load and clean up rpoC data
map <- read.table("~/domhain_RNAseq/map.txt", header=T, sep="\t")
map$aliquot_type <- sub("-", "", map$aliquot_type)
row.names(map) <- map$sample_id
# sequence table
seqtab <- read.table("~/domhain_RNAseq/homd_rpoc_suzanne-7.2.24/sequence_table.merged.txt", header=T, sep="\t", row.names=1)
tax <- read.table("~/domhain_RNAseq/homd_rpoc_suzanne-7.2.24/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
notinmeta <- setdiff(row.names(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), row.names(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw
# should both come back as character(0) 
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)))
ps.dat

# first need to filter our dataframe by genera of interest (going off of RNA seq data)
filter_gen <- unique(c(count.twice.gen, count.triple.gen))
# filter genera of interest at the V8 level
ps.dat.filt <- prune_taxa(tax_table(ps.dat)[, "V7"] %in% filter_gen, ps.dat)
ps.dat.filt 
```

### 10. Run histograms for teeth sampled ~200 days apart first

```R
# filter samples to only include our twice samples
ps.dat.twice <- prune_samples(twice, ps.dat.filt)
ps.dat.twice
# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.twice, taxrank = "V7")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))
# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V7)
asvtab <- as.data.frame(asvtab)
# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by genus (V25) and sum across rows
merge.asvtab <- asvtab %>% group_by(V25) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# rename columns
merge.asvtab <- merge.asvtab %>% rename(Genus = V25)

# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")
# Collapse and sum rows
collapsed <- merge.asvtab %>%
  mutate(Genus = ifelse(Genus %in% keepgroup, Genus, "Other")) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
```

### 11. Stacked barchart for 200 day sampled teeth

```R
# set color palette
gencols <- c(Actinomyces = "#F66140", 
			Kingella = "#1FCAD4", 
			Leptotrichia = "#97002E", 
			Streptococcus = "#16AA48", 
			Treponema = "#F658F8", 
			Other = "#838383")

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)

# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(Genus = factor(Genus, levels = c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.twice %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("rpoC_log10_ADS_200days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Genus)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat rpoC_log10_ADS_200days.pdf")
```

### 12. Run histograms for teeth sampled ~400 days apart 

```R
# filter samples to only include our twice samples
ps.dat.triple <- prune_samples(triple, ps.dat.filt)
ps.dat.triple
# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.triple, taxrank = "V7")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))
# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V7)
asvtab <- as.data.frame(asvtab)
# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by genus (V13) and sum across rows
merge.asvtab <- asvtab %>% group_by(V13) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))
# rename columns
merge.asvtab <- merge.asvtab %>% rename(Genus = V13)

# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")
# Collapse and sum rows
collapsed <- merge.asvtab %>%
  mutate(Genus = ifelse(Genus %in% keepgroup, Genus, "Other")) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))
```

### 13. Stacked barchart for 400 day sampled teeth

```R
# set color palette
gencols <- c(Actinomyces = "#F66140", 
			Kingella = "#1FCAD4", 
			Leptotrichia = "#97002E", 
			Streptococcus = "#16AA48", 
			Treponema = "#F658F8", 
			Other = "#838383")

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)

# reorder genera
collapsed_long <- collapsed_long %>%
  mutate(Genus = factor(Genus, levels = c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema", "Other")))

# Convert row names to a column
map <- map.triple %>%
  rownames_to_column(var = "sample")
# merge with metadata
collapsed_long <- collapsed_long %>%
  left_join(map, by = "sample")
# order samples by HIV status
collapsed_long <- collapsed_long %>%
  mutate(sample = factor(sample, levels = unique(sample[order(hiv_status)])))

# get positions to add lines between samples
sample_positions <- seq(3.5, length(unique(collapsed_long$sample)) -0.5, by = 3)

pdf("rpoC_log10_ADS_400days.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Genus)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = gencols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat rpoC_log10_ADS_400days.pdf")
```







