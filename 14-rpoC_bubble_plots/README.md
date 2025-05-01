## Bubble plots of average rpoC counts across taxa of interest

### 1. Load environment

```bash
conda activate 2024-HIV_RNASeq
cd ~/domhain_RNAseq/10-rpoC_bubble_plots
```

### 2. Load libraries

```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(phyloseq)
library(vegan)
```

### 3. Load in rpoC ASV count data and filter out species of interest

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
# change sp. HMT to oral taxon to match nomenclature differences
filter_sp <- c("Streptococcus_anginosus", "Streptococcus_australis", "Streptococcus_cristatus", "Streptococcus_gordonii", "Streptococcus_sanguinis", "Streptococcus_sp._oral_taxon_056", "Streptococcus_mitis", "Streptococcus_oralis", "Streptococcus_parasanguinis", "Streptococcus_salivarius", "Streptococcus_sp._oral_taxon_066", "Streptococcus_constellatus", "Streptococcus_intermedius", "Actinomyces_naeslundii", "Actinomyces_oris", "Actinomyces_johnsonii", "Actinomyces_sp._oral_taxon_169", "Actinomyces_sp._oral_taxon_170", "Actinomyces_sp._oral_taxon_175", "Actinomyces_viscosus", "Cryptobacterium_curtum", "Cutibacterium_acnes", "Kingella_oralis", "Oribacterium_asaccharolyticum")
# do these actually exist in my phyloseq object?
species_df <- data.frame(QuerySpecies = filter_sp, stringsAsFactors = FALSE)
# Fuzzy match with taxonomy table
tax_levels <- paste0("V", 2:13)
# Find out if these species exist anywhere in the taxonomy and pull asvids
tax_df <- as.data.frame(tax_table(ps.dat)) %>%
  rownames_to_column("ASVID")

# Find matching ASVs for each query species
result_list <- lapply(species_df$QuerySpecies, function(x) {
  # Find rows where the species appears in any taxonomic level
  matches <- tax_df %>%
    filter(if_any(all_of(tax_levels), ~ . == x))
  
  if(nrow(matches) > 0) {
    data.frame(QuerySpecies = x, 
              ASVID = matches$ASVID,
              TaxLevel = apply(matches[, tax_levels], 1, function(row) {
                names(which(row == x))[1]
              }))
  } else {
    NULL
  }
})

# Combine results
matches_df <- bind_rows(result_list)
# Get all unique matching ASV IDs
matching_asvs <- unique(matches_df$ASVID)
# filter by matching ASVIDs
ps.dat.filt <- prune_taxa(matching_asvs, ps.dat)
ps.dat.filt
# head(tax_table(ps.dat.filt))
```

### 4. Prevalence filtering and conditional mean across HIV status groups

Note: mean abundance is not a good way of looking at microbiome data because of a high prevalence of zeros -- therefore the mean abundance by itself can be misleading as it 1) dilutes true biological signals by averaging with potentially higher numbers of zeros in the data, 2) it doesn't distinguish between true absence versus technical zeros (undetected but present) and 3) it underrepresents low-abundance but consistently present taxa

Instead, we are using a log10 transformed conditional mean --> average abundance of microbial taxa calculated only from samples where the taxon is actually present (i.e. ignoring zeros)

```R
# first need to get a standardized species column (since the species show up in different levels)
matches_df <- matches_df %>%
  mutate(StandardizedSpecies = QuerySpecies)

# First calculate total number of samples in each comparison group
sample_counts <- sample_data(ps.dat.filt) %>% 
  as_tibble() %>% 
  count(hiv_status, tooth_health, name = "total_samples")
# Conditional mean abundance
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status", "tooth_health")) %>%
  # Calculate metrics
  group_by(Species, hiv_status, tooth_health) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
```

### 5. Create bubble chart showing abundance across HIV status groups, split by H vs D

```R
df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)

# I want to add in missing species that do not show up in our rpoC data to make the figures comparable
df <- df %>%
  # Ensure all species are included (even missing ones)
  complete(
    Species = filter_sp,
    hiv_status = unique(df$hiv_status),
    tooth_health = unique(df$tooth_health),
    fill = list(conditional_mean = 0, prevalence = 0)
  ) %>%
  # Reapply factor levels to Species
  mutate(Species = factor(Species, levels = filter_sp))

# order by species level
df <- df[order(df$Species),]
# set levels of x axis
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))
df <- df %>%
  filter(tooth_health != "E")
# and order of grid
df$tooth_health <- factor(df$tooth_health, levels = c("H", "D"))

pdf("HvD_rpoC_bubble_plot.pdf", width = 10)
ggplot(df,
  aes(
    x = hiv_status,
    y = Species,
    size = log10(conditional_mean + 1),
    color = tooth_health
  )
) +
  geom_point(alpha = ifelse(df$conditional_mean == 0, 0.1, 0.7)) +
  scale_size(range = c(0.5, 10), name = "Log10 Conditional Mean") +  
  scale_y_discrete(limits = rev) +  # Reverse y-axis order (top-to-bottom)
  scale_color_manual(values = c("#22A146", "#9B002F")) +
  labs(
    x = "Source",
    y = "Species",
    color = "Tooth Health"
  ) +
  facet_grid(. ~ tooth_health, switch = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Improve x-axis label readability
  )
dev.off()
system("/home/allie/.iterm2/imgcat HvD_rpoC_bubble_plot.pdf")
```

### 6. Create bubble chart showing abundance across HIV status groups, split by HCF vs HCD

```R
filter_sp <- c("Streptococcus_australis", "Streptococcus_cristatus", "Streptococcus_gordonii", "Streptococcus_mitis", "Streptococcus_oralis", "Streptococcus_sanguinis", "Streptococcus_constellatus", "Streptococcus_sinensis", "Streptococcus_anginosus", "Streptococcus_parasanguinis", "Streptococcus_intermedius", "Streptococcus_salivarius", "Streptococcus_sp._oral_taxon_056", "Actinomyces_naeslundii", "Actinomyces_oris", "Actinomyces_sp._oral_taxon_170", "Actinomyces_viscosus", "Actinomyces_sp._oral_taxon_175", "Treponema_medium", "Treponmea_vincentii", "Bulleidia_extructa", "Fusobaterium_sp._oral_taxon_370", "Solobacterium_moorei")

# First calculate total number of samples in each comparison group
sample_counts <- sample_data(ps.dat.filt) %>% 
  as_tibble() %>% 
  count(hiv_status, aliquot_type, name = "total_samples")
# Conditional mean abundance
abundmean <- psmelt(ps.dat.filt) %>%
  # Join with standardized species names
  left_join(
    matches_df %>% select(OTU = ASVID, Species = QuerySpecies),
    by = "OTU"
  ) %>%
  # Filter for target species
  filter(Species %in% filter_sp) %>% 
  # Join with sample counts
  left_join(sample_counts, by = c("hiv_status", "aliquot_type")) %>%
  # Calculate metrics
  group_by(Species, hiv_status, aliquot_type) %>%
  summarize(
    conditional_mean = mean(Abundance[Abundance > 0]),  # Mean of non-zero values
    prevalence = sum(Abundance > 0) / first(total_samples),  # Prevalence calculation
    .groups = "drop"
  ) %>%
  # Replace NaN (from 0/0) with 0
  mutate(conditional_mean = ifelse(is.nan(conditional_mean), 0, conditional_mean))
df <- abundmean
# reorder species column to match RNA seq figures
df$Species <- factor(df$Species, levels = filter_sp)

# I want to add in missing species that do not show up in our rpoC data to make the figures comparable
df <- df %>%
  # Ensure all species are included (even missing ones)
  complete(
    Species = filter_sp,
    hiv_status = unique(df$hiv_status),
    aliquot_type = unique(df$aliquot_type),
    fill = list(conditional_mean = 0, prevalence = 0)
  ) %>%
  # Reapply factor levels to Species
  mutate(Species = factor(Species, levels = filter_sp))

# order by species level
df <- df[order(df$Species),]
# set levels of x axis
df$hiv_status <- factor(df$hiv_status, levels = c("HUU", "HEU", "HI"))
df <- df %>%
  filter(aliquot_type == "HCF" | aliquot_type == "HCD")
# and order of grid
df$aliquot_type <- factor(df$aliquot_type, levels = c("HCF", "HCD"))

pdf("HCFvHCD_rpoC_bubble_plot.pdf", width = 10)
ggplot(df,
  aes(
    x = hiv_status,
    y = Species,
    size = log10(conditional_mean + 1),
    color = aliquot_type
  )
) +
  geom_point(alpha = ifelse(df$conditional_mean == 0, 0.1, 0.7)) +
  scale_size(range = c(0.5, 10), name = "Log10 Conditional Mean") +  
  scale_y_discrete(limits = rev) +  # Reverse y-axis order (top-to-bottom)
  scale_color_manual(values = c("#22A146", "#9B002F")) +
  labs(
    x = "Source",
    y = "Species",
    color = "Aliquot Type"
  ) +
  facet_grid(. ~ aliquot_type, switch = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Improve x-axis label readability
  )
dev.off()
system("/home/allie/.iterm2/imgcat HCFvHCD_rpoC_bubble_plot.pdf")
```

### 7. Finally, I want some ordination plots that show the different categories, highlighting the ones that we've highlighted here 

```R
# clr normalization
ps.clr <- microbiome::transform(ps.dat, "clr") 
# CLR-transformed PCA by tooth health
ord.pca <- ordinate(
  ps.clr,
  method = "RDA",  # Redundancy Analysis
  formula = ~ tooth_health,
  distance = "euclidean"
)

# Plot
pdf("clr_pca.tooth_health.pdf")
plot_ordination(
  ps.clr, 
  ord.pca, 
  color = "tooth_health",  # Map groups to border color (temporarily)
  title = "PCA (CLR-transformed)"
) +
  geom_point(
    size = 3,
    shape = 21,                     # Fillable circle
    stroke = 0.5,                   # Border thickness
    aes(fill = tooth_health),       # Map fill to groups
    show.legend = TRUE              # Ensure legend appears
  ) +
  scale_fill_manual(
    values = c("H" = "#22A146", "D" = "#9B002F", "E" = "#F0F032"),  # Inner colors
    name = "Tooth health"
  ) +
  scale_color_manual(
    values = c("H" = "black", "D" = "black", "E" = "black"),      # Force borders to black
    guide = "none"                                 # Hide color legend
  ) +
  theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat clr_pca.tooth_health.pdf")

# Calculate permanova results
dist_mat <- phyloseq::distance(ps.clr, "euclidean")
sampledf <- data.frame(sample_data(ps.clr))
adonis2(dist_mat ~ tooth_health, data = sampledf, permutations = 999)
# adonis2(formula = dist_mat ~ tooth_health, data = sampledf, permutations = 999)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     2    13273 0.02277 1.0484  0.188
# Residual 90   569709 0.97723
# Total    92   582983 1.00000

# CLR-transformed PCA by aliquot type
ord.pca <- ordinate(
  ps.clr,
  method = "RDA",  # Redundancy Analysis
  formula = ~ aliquot_type,
  distance = "euclidean"
)

# Plot
pdf("clr_pca.aliquot_type.pdf")
plot_ordination(
  ps.clr, 
  ord.pca, 
  color = "aliquot_type",  # Map groups to border color (temporarily)
  title = "PCA (CLR-transformed)"
) +
  geom_point(
    size = 3,
    shape = 21,                     # Fillable circle
    stroke = 0.5,                   # Border thickness
    aes(fill = aliquot_type),       # Map fill to groups
    show.legend = TRUE              # Ensure legend appears
  ) +
  scale_fill_manual(
    values = c("DCD" = "#AA0A3B", 
    		   "ECD" = "#F87850", 
    		   "ECE" = "#F0F032",
    		   	"HCD" = "#2F5AC8",
				"HCE" = "#3FD2DC",
				"HCF" = "#24B45A"),
    name = "Tooth and oral health"
  ) +
  scale_color_manual(
    values = c("DCD" = "black", 
    		   "ECD" = "black", 
    		   "ECE" = "black", 
    		   "HCD" = "black", 
    		   "HCE" = "black", 
    		   "HCF" = "black"),      
    guide = "none"                                 # Hide color legend
  ) +
  theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat clr_pca.aliquot_type.pdf")

# Calculate permanova results
adonis2(dist_mat ~ aliquot_type, data = sampledf, permutations = 999)
# adonis2(formula = dist_mat ~ aliquot_type, data = sampledf, permutations = 999)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     5    37448 0.06423 1.1944  0.001 ***
# Residual 87   545535 0.93577
# Total    92   582983 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

One more at the HIV status level for supplement

```R
# CLR-transformed PCA by hiv status
ord.pca <- ordinate(
  ps.clr,
  method = "RDA",  # Redundancy Analysis
  formula = ~ hiv_status,
  distance = "euclidean"
)

# Plot
pdf("clr_pca.hiv_status.pdf")
plot_ordination(
  ps.clr, 
  ord.pca, 
  color = "hiv_status",  # Map groups to border color (temporarily)
  title = "PCA (CLR-transformed)"
) +
  geom_point(
    size = 3,
    shape = 21,                     # Fillable circle
    stroke = 0.5,                   # Border thickness
    aes(fill = hiv_status),       # Map fill to groups
    show.legend = TRUE              # Ensure legend appears
  ) +
  scale_fill_manual(
    values = c("HI" = "black", "HEU" = "grey", "HUU" = "white"),  # Inner colors
    name = "Tooth health"
  ) +
  scale_color_manual(
    values = c("HI" = "black", "HEU" = "black", "HUU" = "black"),      # Force borders to black
    guide = "none"                                 # Hide color legend
  ) +
  theme_minimal()
dev.off()
system("/home/allie/.iterm2/imgcat clr_pca.hiv_status.pdf")

# Calculate permanova results
adonis2(dist_mat ~ hiv_status, data = sampledf, permutations = 999)
# adonis2(formula = dist_mat ~ hiv_status, data = sampledf, permutations = 999)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     2    18026 0.03092 1.4358  0.001 ***
# Residual 90   564957 0.96908
# Total    92   582983 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```