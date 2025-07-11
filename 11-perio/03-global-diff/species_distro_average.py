import pandas as pd
import numpy as np
from collections import defaultdict

# Step 1: Read the read counts file
read_counts_file = '~/rna_dohmain/11-perio/02-pgap/gene_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None, header=0)
read_counts_df = read_counts_df.iloc[:, :-1]  # Remove last column if not needed
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping and genome-to-species relationships
locus_to_taxonomy = {}
genome_to_species = defaultdict(set)  # Maps species to a set of genomes

locus_taxonomy_file = '/home/suzanne/rna_dohmain/11-perio/02-pgap/gene_annots.txt'
with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        columns = line.strip().split('\t')
        if len(columns) >= 7:
            genome, tag, gene, product, genus, species, strain = columns[:7]
            locus_to_taxonomy[tag] = species
            genome_to_species[species].add(genome)

# Step 3: Count the number of genomes per species
genomes_per_species = {species: len(genomes) for species, genomes in genome_to_species.items()}

# Step 4: Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Step 5: Iterate over each sample to compute normalized read percentages
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()
    print(f"Processing sample: {sample_id}")

    total_reads_sample_by_species = {species: 0 for species in species_list}

    for gene in read_counts_df.index:
        species = locus_to_taxonomy.get(gene)
        if species:
            total_reads_sample_by_species[species] += read_counts_df.loc[gene, sample_id]

    # Normalize by number of genomes and calculate percentage
    for species, total_reads in total_reads_sample_by_species.items():
        num_genomes = genomes_per_species.get(species, 1)  # Avoid division by zero
        normalized_reads = total_reads / num_genomes

        if total_reads_sample == 0:
            percentage_reads = np.nan
        else:
            percentage_reads = (normalized_reads / total_reads_sample) * 100

        percentage_reads_by_sample.loc[sample_id, species] = percentage_reads

# Finalize DataFrame
percentage_reads_by_sample.index.name = "sample"

# Save the results to a file
percentage_reads_by_sample.to_csv('species_reads_by_genome.txt', sep="\t", index=True, header=True)

print("Percentage of reads (normalized by number of genomes per species) saved to 'species_reads_by_genome.txt'")
