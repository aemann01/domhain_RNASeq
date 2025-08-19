import pandas as pd
import numpy as np  # For handling NaN values

# Step 1: Read the read counts file (assuming it's a tabular format)
read_counts_file = '~/rna_dohmain/11-perio/02-pgap/gene_counts.txt'
read_counts_df = pd.read_csv(read_counts_file, sep='\t', index_col=None, header = 0)  # Assuming sample IDs are in the first column
# read_counts_df = read_counts_df.iloc[:, :-1]  # Remove last column if it's not needed
read_counts_df.set_index('Geneid', inplace=True)
read_counts_df.index.name = None  # Remove the name of the index

# Step 2: Read the locus tag to taxonomy mapping file into a dictionary
locus_taxonomy_file = '/home/suzanne/rna_dohmain/11-perio/02-pgap/gene_annots.txt'
locus_to_taxonomy = {}

with open(locus_taxonomy_file, 'r') as f:
    for line in f:
        # Split the line into columns and avoid duplicate 'species' assignment
        columns = line.strip().split('\t')
        if len(columns) >= 7:  # Check if the line has enough columns
            genome, tag, gene, product, genus, species, strain = columns[:7]
            locus_to_taxonomy[tag] = species  # Using 'tag' as the locus tag key

# Initialize DataFrame to store percentage of reads for each species by sample
species_list = list(set(locus_to_taxonomy.values()))  # List of unique species
percentage_reads_by_sample = pd.DataFrame(index=read_counts_df.columns, columns=species_list)

# Calculate total reads for each species across all samples
total_reads_by_species = {species: 0 for species in species_list}

# Iterate over each sample
for sample_id in read_counts_df.columns:
    total_reads_sample = read_counts_df[sample_id].sum()  # Total reads for the current sample
    print(f"Processing sample: {sample_id}")

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

# Set index name for the output
percentage_reads_by_sample.index.name = "sample"  # Remove the name of the index

# Save the results to a file
percentage_reads_by_sample.to_csv('species_reads.txt', sep="\t", index=True, header=True)
