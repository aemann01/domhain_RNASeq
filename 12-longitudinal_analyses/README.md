## Longitudinal changes in ADS activity and abundance of ADS competent bacteria on single teeth 

### 1. Load environment

```bash
cd ~/domhain_RNAseq/09-longitudinal_analyses
conda activate 2024-HIV_RNASeq
R
```

### 2. Load libraries

```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(stringr)
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
# get list of species to pull from rpoC data later
keepsp.twice <- unique(paste(count.twice$Genus, count.twice$Species, sep="_"))

# collapse by genus and sum across rows
merge.count <- count.twice %>% group_by(Genus) %>% summarize(across(where(is.numeric), \(x) sum(x, na.rm=TRUE)))
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

### 5a. Set color palette for genera

```R
# set color palette
gencols <- c(Actinomyces = "#F66140", 
			Kingella = "#1FCAD4", 
			Leptotrichia = "#97002E", 
			Streptococcus = "#16AA48", 
			Treponema = "#0C5B6F", 
			Other = "#838383")
```

### 5b. Set color palette for species

```R
# set species colors
spcols <- c(Actinomyces = "#800026",
	Actinomyces_graevenitzii = "#bd0026",
	Actinomyces_johnsonii = "#e31a1c",
	Actinomyces_naeslundii = "#fc4e2a",
	Actinomyces_oris = "#fd8d3c",
	Actinomyces_sp._oral_taxon_169 = "#feb24c",
	Actinomyces_sp._oral_taxon_170 = "#fed976",
	Actinomyces_sp._oral_taxon_171 = "#ffeda0",
	Actinomyces_sp._oral_taxon_175 = "#ffffcc",
	Actinomyces_viscosus = "#f16913",
	Kingella = "#08306b",
	Kingella_oralis = "#08519c",
	Leptotrichia = "#67001f",
	Leptotrichia_sp._oral_taxon_212 = "#980043",
	Leptotrichia_sp._oral_taxon_215 = "#ce1256",
	Streptococcus = "#00441b", 
	Streptococcus_aginosus = "#006d2c",
	Streptococcus_australis = "#238b45",
	Streptococcus_constellatus = "#41ae76",
	Streptococcus_cristatus = "#66c2a4",
	Streptococcus_gordonii = "#99d8c9",
	Streptococcus_infantis = "#ccece6",
	Streptococcus_intermedius = "#004529",
	Streptococcus_mitis = "#006837",
	Streptococcus_oralis = "#238443",
	Streptococcus_parasanguinis = "#41ab5d",
	Streptococcus_pneumoniae = "#78c679",
	Streptococcus_salivarius = "#addd8e",
	Streptococcus_sanguinis = "#d9f0a3",
	Streptococcus_sinensis = "#f7fcb9",
	Streptococcus_sp._oral_taxon_066 = "#016c59",
	Streptococcus_sp._oral_taxon_066 = "#02818a",
	Treponema = "#4d004b",
	Treponema_medium = "#810f7c",
	Treponema_vincentii = "#88419d",
	Other = "#838383")
```

### 6. Create stacked histograms for each sample colored by genus (200 days RNA data)

```R
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

### 6a. Create stacked histograms for each sample colored by species (200 days RNA data)

```R
# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")

# first need to collapse dataframe by Genus_Species
collapsed <- count.twice %>%
  mutate(
    Genus_grouped = case_when(Genus %in% keepgroup ~ Genus, TRUE ~ "Other"),
    Genus_species = if_else(Genus_grouped == "Other", "Other", paste(Genus, Species, sep = "_"))) %>%
  group_by(Genus_species) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# replace HMT with oral_taxon
collapsed <- collapsed %>%
	mutate(Genus_species = gsub("sp._HMT", "sp._oral_taxon", Genus_species))

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus_species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder species
collapsed_long <- collapsed_long %>%
  mutate(Species = factor(Genus_species, levels = c( 
  	"Actinomyces_graevenitzii", 
  	"Actinomyces_johnsonii",
  	"Actinomyces_naeslundii", 
  	"Actinomyces_oris", 
  	"Actinomyces_sp._oral_taxon_169",
  	"Actinomyces_sp._oral_taxon_170",
  	"Actinomyces_sp._oral_taxon_171",
  	"Actinomyces_sp._oral_taxon_175",
  	"Actinomyces_viscosus",
  	"Kingella_oralis", 
  	"Leptotrichia_sp._oral_taxon_212", 
  	"Leptotrichia_sp._oral_taxon_215", 
  	"Streptococcus_aginosus",
  	"Streptococcus_australis",
  	"Streptococcus_constellatus",
  	"Streptococcus_cristatus", 
  	"Streptococcus_gordonii", 
  	"Streptococcus_infantis",
  	"Streptococcus_intermedius",
  	"Streptococcus_mitis", 
  	"Streptococcus_oralis", 
  	"Streptococcus_parasanguinis", 
  	"Streptococcus_pneumoniae",
  	"Streptococcus_salivarius",
  	"Streptococcus_sanguinis", 
  	"Streptococcus_sinensis", 
  	"Streptococcus_sp._HMT_056",
  	"Streptococcus_sp._HMT_066",
  	"Treponema_medium", 
  	"Treponema_vincentii", 
  	"Other")))

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

pdf("RNA_log10_ADS_200days-species.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = spcols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat RNA_log10_ADS_200days-species.pdf")

```

### 7. Create stacked histograms for each sample colored by genus (400 days RNA data)

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
# get list of species to pull from rpoC data later
keepsp.triple <- unique(paste(count.triple$Genus, count.triple$Species, sep="_"))

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

### 7a. Create stacked histograms for each sample colored by species (400 days RNA data)

```R
# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")

# first need to collapse dataframe by Genus_Species
collapsed <- count.triple %>%
  mutate(
    Genus_grouped = case_when(Genus %in% keepgroup ~ Genus, TRUE ~ "Other"),
    Genus_species = if_else(Genus_grouped == "Other", "Other", paste(Genus, Species, sep = "_"))) %>%
  group_by(Genus_species) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# replace HMT with oral_taxon
collapsed <- collapsed %>%
	mutate(Genus_species = gsub("sp._HMT", "sp._oral_taxon", Genus_species))

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Genus_species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder species
collapsed_long <- collapsed_long %>%
  mutate(Species = factor(Genus_species, levels = c( 
  	"Actinomyces_johnsonii",
  	"Actinomyces_naeslundii", 
  	"Actinomyces_oris", 
  	"Actinomyces_sp._oral_taxon_169",
  	"Actinomyces_sp._oral_taxon_170",
  	"Actinomyces_sp._oral_taxon_171",
  	"Actinomyces_sp._oral_taxon_175",
  	"Actinomyces_viscosus",
  	"Kingella_oralis", 
  	"Leptotrichia_sp._oral_taxon_212", 
  	"Leptotrichia_sp._oral_taxon_215", 
  	"Streptococcus_aginosus",
  	"Streptococcus_australis",
  	"Streptococcus_constellatus",
  	"Streptococcus_cristatus", 
  	"Streptococcus_gordonii", 
  	"Streptococcus_infantis",
  	"Streptococcus_intermedius",
  	"Streptococcus_mitis", 
  	"Streptococcus_oralis", 
  	"Streptococcus_parasanguinis", 
  	"Streptococcus_pneumoniae",
  	"Streptococcus_salivarius",
  	"Streptococcus_sanguinis", 
  	"Streptococcus_sinensis", 
  	"Streptococcus_sp._HMT_056",
  	"Streptococcus_sp._HMT_066",
  	"Treponema_medium", 
  	"Treponema_vincentii", 
  	"Other")))

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
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("RNA_log10_ADS_400days-species.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = spcols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat RNA_log10_ADS_400days-species.pdf")
```

### 8. Load in rpoC ASV count data 

```R
# load and clean up rpoC data
map <- read.table("~/domhain_RNAseq/map.txt", header=T, sep="\t")
map$aliquot_type <- sub("-", "", map$aliquot_type)
row.names(map) <- map$sample_id
# sequence table
seqtab <- read.table("~/domhain_RNAseq/11-rpoc_processing/sequence_table.merged.txt", header=T, sep="\t", row.names=1)
tax <- read.table("~/domhain_RNAseq/11-rpoc_processing/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
notinmeta <- setdiff(row.names(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), row.names(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw
# should both come back as character(0) 
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)))
ps.dat

# # generate a quick rarefaction curve plot for supplementary figures
# library(ranacapa)
# rarefaction_plot <- ggrare(
#     physeq = ps.dat,  # Replace with your phyloseq object
#     step = 100,                     # Step size for rarefaction
#     color = "hiv_status",           # Grouping variable (metadata column)
#     se = TRUE                      # Include standard error
# ) +
#     theme_bw() +
#     labs(title = "Rarefaction Curve", x = "Sequencing Depth", y = "Observed ASVs")
# pdf("rarefaction_plot.pdf")
# rarefaction_plot
# dev.off()
# system("/home/allie/.iterm2/imgcat rarefaction_plot.pdf")
```

### 9. Function for filtering out species from ASV table

```R
# filter phyloseq object to only include those species we want to keep (in any column)
filter_ps_by_sp <- function(ps.dat, keepsp, keepgen){
	# convert tax table to dataframe
	tax_df <- as.data.frame(as(tax_table(ps.dat), "matrix"))
	last_col <- names(tax_df)[ncol(tax_df)]
	# identify taxa that match species vector in any column
	keeptax1 <- rownames(tax_df)[apply(tax_df, 1, function(row) any(row %in% keepsp))]
	# keep taxa if they are at the genus level in the last column (e.g., unassigned to species level)
	keeptax2 <- rownames(tax_df)[tax_df[[last_col]] %in% keepgen]
	# combine ASV names 
	keeptax <- unique(c(keeptax1, keeptax2))

	# subset phyloseq object
	ps.dat.filt <- prune_taxa(keeptax, ps.dat)
	return(ps.dat.filt)
}
```

### 10. Run histograms for teeth sampled ~200 days apart first at genus level (rpoC data)

```R
# fix species names so that they match between the RNA and rpoC datasets
keepsp.twice <- gsub("_sp._HMT", "_sp._oral_taxon", keepsp.twice)
# get list of genera to keep if they are unassigned at the species level
keepgen.twice <- unique(sub("_.*", "", keepsp.twice))

# need to clean up taxonomy at L10 to be able to collapse correctly, remove trailing strain identifiers (but leaving sp._oral_taxon)
tax_table <- as.data.frame(tax_table(ps.dat), stringsAsFactors=F)
tax_table$V10 <- ifelse(
  str_detect(tax_table$V10, "_sp\\._oral_taxon_\\d+"),  # if oral taxon, leave alone
  tax_table$V10,  
  str_replace(tax_table$V10, "^(.*?_.*?)(_.*)?$", "\\1")  # otherwise, remove trailing strain identifier
)

# replace V10
tax_table(ps.dat) <- as.matrix(tax_table)

# run filter function
ps.dat.filt <- filter_ps_by_sp(ps.dat, keepsp.twice, keepgen.twice)

# filter samples to only include our twice samples
ps.dat.twice <- prune_samples(twice, ps.dat.filt)
ps.dat.twice

# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.twice, taxrank = "V8")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))

# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V8)
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

### 11. Stacked barchart for 200 day sampled teeth at genus level (rpoC data)

```R
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

### 12. Run histograms for teeth sampled ~200 days apart at the species level (rpoC data)

```R
# aggregate counts by species level
ps.dat.sp <- tax_glom(ps.dat.twice, taxrank = "V10")

# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.sp))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.sp))
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true

# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V10)
asvtab <- as.data.frame(asvtab)

# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by species and sum across rows
merge.asvtab <- asvtab %>% group_by(V25) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))

# rename columns
merge.asvtab <- merge.asvtab %>% rename(Species = V25)
# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")

# collapse and sum rows
collapsed <- merge.asvtab %>%
	mutate(Genus = str_extract(Species, "^[^_]+"),
		Species = ifelse(Genus %in% keepgroup, Species, "Other")) %>%
	select(-Genus) %>%
	group_by(Species) %>%
	summarize(across(where(is.numeric), sum, na.rm = TRUE))

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder species
collapsed_long <- collapsed_long %>%
  mutate(Species = factor(Species, levels = c("Actinomyces", 
  	"Actinomyces_graevenitzii", 
  	"Actinomyces_naeslundii", 
  	"Actinomyces_oris", 
  	"Actinomyces_sp._oral_taxon_169",
  	"Actinomyces_sp._oral_taxon_170",
  	"Kingella", 
  	"Leptotrichia",
  	"Leptotrichia_sp._oral_taxon_212", 
  	"Leptotrichia_sp._oral_taxon_215", 
  	"Streptococcus",
  	"Streptococcus_aginosus",
  	"Streptococcus_australis",
  	"Streptococcus_constellatus",
  	"Streptococcus_cristatus", 
  	"Streptococcus_gordonii", 
  	"Streptococcus_infantis",
  	"Streptococcus_intermedius",
  	"Streptococcus_mitis", 
  	"Streptococcus_oralis", 
  	"Streptococcus_parasanguinis", 
  	"Streptococcus_pneumoniae",
  	"Streptococcus_salivarius",
  	"Streptococcus_sanguinis", 
  	"Streptococcus_sinensis", 
  	"Treponema",
  	"Treponema_medium", 
  	"Treponema_vincentii", 
  	"Other")))

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

pdf("rpoC_log10_ADS_200days-species.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = spcols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat rpoC_log10_ADS_200days-species.pdf")
```

### 12. Run histograms for teeth sampled ~400 days apart at the genus level (rpoC)

```R
# fix species names so that they match between the RNA and rpoC datasets
keepsp.triple <- gsub("_sp._HMT", "_sp._oral_taxon", keepsp.triple)
# get list of genera to keep if they are unassigned at the species level
keepgen.triple <- unique(sub("_.*", "", keepsp.triple))

# need to clean up taxonomy at L10 to be able to collapse correctly, remove trailing strain identifiers (but leaving sp._oral_taxon)
tax_table <- as.data.frame(tax_table(ps.dat), stringsAsFactors=F)
tax_table$V10 <- ifelse(
  str_detect(tax_table$V10, "_sp\\._oral_taxon_\\d+"),  # if oral taxon, leave alone
  tax_table$V10,  
  str_replace(tax_table$V10, "^(.*?_.*?)(_.*)?$", "\\1")  # otherwise, remove trailing strain identifier
)

# replace V10
tax_table(ps.dat) <- as.matrix(tax_table)

# run filter function
ps.dat.filt <- filter_ps_by_sp(ps.dat, keepsp.triple, keepgen.triple)

# filter samples to only include our triple samples
ps.dat.triple <- prune_samples(triple, ps.dat.filt)
ps.dat.triple

# aggregate counts by genus
ps.dat.gen <- tax_glom(ps.dat.triple, taxrank = "V8")
# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.gen))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.gen))

# merge tax table and asvtab by ASVID
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true
# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V8)
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

### 13. Stacked barchart for 400 day sampled teeth (genus level, rpoC)

```R
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

### 12. Run histograms for teeth sampled ~400 days apart at the species level (rpoC data)

```R
# aggregate counts by species level
ps.dat.sp <- tax_glom(ps.dat.triple, taxrank = "V10")

# extract ASV table
asvtab <- as.data.frame(otu_table(ps.dat.sp))
asvtab <- t(asvtab)

# extract taxonomic information to map to ASV ids
tax_table <- as.data.frame(tax_table(ps.dat.sp))
# check that locus tags match between the two dataframes
table(rownames(tax_table)==rownames(asvtab)) # should return all true

# if all are true, merge together
asvtab <- cbind(asvtab, tax_table$V10)
asvtab <- as.data.frame(asvtab)

# convert all but last column to numeric
asvtab <- asvtab %>%
  mutate(across(-ncol(asvtab), as.numeric))
# merge by species and sum across rows
merge.asvtab <- asvtab %>% group_by(V13) %>% summarize(across(where(is.numeric), sum, na.rm=TRUE))

# rename columns
merge.asvtab <- merge.asvtab %>% rename(Species = V13)
# groups to collapse by
keepgroup <- c("Actinomyces", "Kingella", "Leptotrichia", "Streptococcus", "Treponema")

# collapse and sum rows
collapsed <- merge.asvtab %>%
	mutate(Genus = str_extract(Species, "^[^_]+"),
		Species = ifelse(Genus %in% keepgroup, Species, "Other")) %>%
	select(-Genus) %>%
	group_by(Species) %>%
	summarize(across(where(is.numeric), sum, na.rm = TRUE))

# long format
collapsed_long <- collapsed %>%
  pivot_longer(cols = -Species, names_to = "sample", values_to = "count")
# convert to log10 values
collapsed_long <- collapsed_long %>%
  mutate(log10_count = log10(count + 1))  # Add 1 to avoid log10(0)
# reorder species
collapsed_long <- collapsed_long %>%
  mutate(Species = factor(Species, levels = c("Actinomyces", 
  	"Actinomyces_naeslundii", 
  	"Actinomyces_oris", 
  	"Actinomyces_sp._oral_taxon_169",
  	"Actinomyces_sp._oral_taxon_170",
  	"Kingella", 
  	"Leptotrichia",
  	"Leptotrichia_sp._oral_taxon_212", 
  	"Leptotrichia_sp._oral_taxon_215", 
  	"Streptococcus",
  	"Streptococcus_aginosus",
  	"Streptococcus_australis",
  	"Streptococcus_constellatus",
  	"Streptococcus_cristatus", 
  	"Streptococcus_gordonii", 
  	"Streptococcus_infantis",
  	"Streptococcus_intermedius",
  	"Streptococcus_mitis", 
  	"Streptococcus_oralis", 
  	"Streptococcus_parasanguinis", 
  	"Streptococcus_pneumoniae",
  	"Streptococcus_salivarius",
  	"Streptococcus_sanguinis", 
  	"Streptococcus_sinensis", 
  	"Treponema",
  	"Treponema_medium", 
  	"Treponema_vincentii", 
  	"Other")))

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
sample_positions <- seq(2.5, length(unique(collapsed_long$sample)) -0.5, by = 2)

pdf("rpoC_log10_ADS_400days-species.pdf", width = 20, height = 5)
ggplot(collapsed_long, aes(x = sample, y = log10_count, fill = Species)) +
	geom_bar(stat="identity", position="stack") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = spcols) +
	geom_vline(xintercept = sample_positions, linetype = "dotted", color = "black", linewidth = 0.75)
dev.off()
system("/home/allie/.iterm2/imgcat rpoC_log10_ADS_400days-species.pdf")
```



















