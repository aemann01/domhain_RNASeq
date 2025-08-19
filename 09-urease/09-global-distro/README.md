# 1. Precent of Global Reads that are A. naeslundii
Get annots file
```sh
cp annotations.merge.txt ../09-urease/09-global-distro
paste annotations.merge.txt <(awk '{print $4, $5}' annotations.merge.txt | sed 's/ /_/') | sed 's/ /\t/g' > temp
mv temp annotations.merge.txt
```
```py
import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/homd_map/read_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.merge.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus_Species


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('species_reads.txt', sep="\t", index=True, header = True)
```
# 2. Get distro of genus
genus_distro.py
```py
import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/homd_map/read_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.merge.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('genus_reads.txt', sep="\t", index=True, header = True)
```
# 3. See gene distro
gene_distro.py
```py
import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/homd_map/read_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.merge.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = gene


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('gene_reads.txt', sep="\t", index=True, header = True)
```
# 4. See rpoC taxonomic distro
```sh
parallel -a <(grep rpoC annotations.merge.txt | awk '{print $1}') -j 7 -k "grep '{}' ../../homd_map/read_counts.txt"> rpoC_counts.txt
cat <(head -n 1 ../../homd_map/read_counts.txt) rpoC_counts.txt > temp
mv temp rpoC_counts.txt
grep rpoC annotations.merge.txt > annotations.rpoc.txt
sed -i 's/rpoC_[12]/rpoC/' annotations.rpoc.txt
cat <(head -n 1 annotations.merge.txt) annotations.rpoc.txt > temp
mv temp annotations.rpoc.txt
```
Per species
```py
import pandas as pd
import numpy as np  # For handling NaN values
import matplotlib.pyplot as plt

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = './rpoC_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.rpoc.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus_Species


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('species_rpoC.txt', sep="\t", index=True, header = True)

pd.set_option('display.max_rows', None)  # None for unlimited rows
pd.set_option('display.max_columns', None)  # None for unlimited columns

# Print the entire DataFrame
print(percentage_reads_by_sample.Actinomyces_naeslundii)

# Reset options to default after use if needed
pd.reset_option('display.max_rows')
pd.reset_option('display.max_columns')
```
Per genus
```py
import pandas as pd
import numpy as np  # For handling NaN values
import matplotlib.pyplot as plt

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = './rpoC_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annotations.rpoc.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index
percentage_reads_by_sample.to_csv('genus_rpoC.txt', sep="\t", index=True, header = True)

pd.set_option('display.max_rows', None)  # None for unlimited rows
pd.set_option('display.max_columns', None)  # None for unlimited columns

# Print the entire DataFrame
print(percentage_reads_by_sample.Actinomyces)

# Reset options to default after use if needed
pd.reset_option('display.max_rows')
pd.reset_option('display.max_columns')
```
# 5. Make plots

rpoC
```R
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
meta <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
counts <- read.csv("species_rpoC.txt", header=T,sep = "\t")
data <- melt(counts)
data <- left_join(data, meta, by = join_by(sample ==  sample_id))

data$genus <- gsub(x = data$variable, pattern = "_.*", replacement = "")

pdf("rna_rpoc.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=sample, y=value,fill=genus),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./rna_rpoc.pdf")

#actinomyces plot
sub_data <- data[data$genus == "Actinomyces",]
pdf("actino_rna_rpoc.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./actino_rna_rpoc.pdf")
#actinomyces naeslundii plot
sub_data <- data[data$variable == "Actinomyces_naeslundii",]
pdf("anaes_rna_rpoc.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./anaes_rna_rpoc.pdf")
pdf("anaes_rna_rpoc.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./anaes_rna_rpoc.hist.pdf")
# Strep sang
sub_data <- data[data$variable == "Streptococcus_sanguinis",]
pdf("ssang_rna_rpoc.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./ssang_rna_rpoc.pdf")
pdf("ssang_rna_rpoc.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./ssang_rna_rpoc.hist.pdf")
mean(sub_data$value)
# Strep sang
sub_data <- data[data$variable == "Streptococcus_oralis",]
pdf("soral_rna_rpoc.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./soral_rna_rpoc.pdf")
pdf("soral_rna_rpoc.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./soral_rna_rpoc.hist.pdf")
mean(sub_data$value)
#Find the top 5 species from rpoC rna

top_species <- data %>%
  group_by(variable) %>%
  summarise(Total_Count = sum(value)/93) %>%
  arrange(desc(Total_Count)) 

print(top_species, n =50)

top_rpoC <- data %>%
  group_by(V8) %>%
  summarise(Total_Count = (sum(Abundance)/93)*100) %>%
  arrange(desc(Total_Count)) 
print(top_rpoC, n =50)

```
Full dataset
```R
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
meta <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
counts <- read.csv("species_reads.txt", header=T,sep = "\t")
data <- melt(counts)
data <- left_join(data, meta, by = join_by(sample ==  sample_id))

data$genus <- gsub(x = data$variable, pattern = "_.*", replacement = "")

pdf("rna_reads.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=sample, y=value,fill=genus),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./rna_reads.pdf")

#actinomyces plot
sub_data <- data[data$genus == "Actinomyces",]
pdf("actino_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./actino_rna_reads.pdf")
#actinomyces naeslundii plot
sub_data <- data[data$variable == "Actinomyces_naeslundii",]
pdf("anaes_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./anaes_rna_reads.pdf")
mean(sub_data$value)

pdf("anaes_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./anaes_rna_reads.hist.pdf")
# Strep sang
sub_data <- data[data$variable == "Streptococcus_sanguinis",]
pdf("ssang_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./ssang_rna_reads.pdf")
pdf("ssang_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./ssang_rna_reads.hist.pdf")
mean(sub_data$value)
# Trep denticola
sub_data <- data[data$variable == "Treponema_denticola",]
pdf("tdent_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./tdent_rna_reads.pdf")
pdf("tdent_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./tdent_rna_reads.hist.pdf")
mean(sub_data$value)
# Porphyromonas_gingivalis
sub_data <- data[data$variable == "Porphyromonas_gingivalis",]
pdf("pging_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./pging_rna_reads.pdf")
pdf("pging_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./pging_rna_reads.hist.pdf")
mean(sub_data$value)
# Trep denticola
sub_data <- data[data$variable == "Tannerella_forsythia",]
pdf("tfors_rna_reads.pdf", width =15, heigh =10)
ggplot(sub_data)+
  geom_bar(aes(x=sample, y=value,fill=variable),stat="identity", position="stack")+
  facet_grid(~ hiv_status, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
      legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./tfors_rna_reads.pdf")
pdf("tfors_rna_reads.hist.pdf", width =15, heigh =10)
ggplot(sub_data, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", aes(y = ..count..)) +
  labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
  theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./tfors_rna_reads.hist.pdf")
mean(sub_data$value)
#Find the top 5 species from rpoC rna

top_species <- data %>%
  group_by(variable) %>%
  summarise(Total_Count = sum(value)/93) %>%
  arrange(desc(Total_Count)) 

print(top_species, n =50)

top_rpoC <- data %>%
  group_by(V8) %>%
  summarise(Total_Count = (sum(Abundance)/93)*100) %>%
  arrange(desc(Total_Count)) 
print(top_rpoC, n =50)
```






```py
import pandas as pd
import numpy as np

# Example DataFrame of read counts with genes as rows
data = {
    'Sample1': [100, 200, 300],
    'Sample2': [150, 250, 350],
    'Sample3': [180, 280, 380]
}
read_counts_df = pd.DataFrame(data, index=['Gene1', 'Gene2', 'Gene3'])

# Example dictionary mapping genes to species
gene_to_species = {
    'Gene1': 'SpeciesA',
    'Gene2': 'SpeciesB',
    'Gene3': 'SpeciesA'
}

# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(gene_to_species.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = gene_to_species.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)


```
```py
import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = 'count_test.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None)  # Assuming sample IDs are in the first column
read_counts_df = read_counts_df.iloc[:, :-1]
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = 'annot_test.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        locus_tag, gene, SEQ_ID, Genus, Species, Genus_Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Genus_Species


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(sample_id)
    
    # Initialize dictionary to store total reads for each species in the current sample
    total_reads_sample_by_species = {species: 0 for species in species_list}
    
    # Iterate over each gene in the read counts DataFrame
    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)  # Get species for the current gene
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]
    
    # Calculate percentage of reads for each species in the current sample
    for species, total_reads in total_reads_sample_by_species.items():
        if total_reads_sample == 0:
            percentage_reads = np.nan  # Handle division by zero or no reads case
        else:
            percentage_reads = (total_reads / total_reads_sample) * 100
        
        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)


```