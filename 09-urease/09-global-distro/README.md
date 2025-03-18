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
        locus_to_taxonomy[locus_tag] = Species


# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
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
        print(percentage_reads_by_sample)

# Print or use percentage_reads_by_sample DataFrame as needed
print("Percentage of reads for each species by sample:")
print(percentage_reads_by_sample)
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
        locus_tag, gene, SEQ_ID, Genus, Species = line.strip().split('\t')
        locus_to_taxonomy[locus_tag] = Species


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