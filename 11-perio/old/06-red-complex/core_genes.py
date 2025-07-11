import os
from Bio import SeqIO
import re  
import gzip
from collections import defaultdict

# get all files that end in .fna
directory = "/home/suzanne/rna_dohmain/11-perio/06-red-complex/tree"  
fna_files = [f for f in os.listdir(directory) if f.endswith(".fna.gz")]

# make function to get gene names
def extract_unique_genes(file_path):
    gene_counts = {}  
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(r'\[gene=([^\]]+)\]', record.description)
            if match:
                gene_name = match.group(1)
                if gene_name in gene_counts:
                    gene_counts[gene_name] += 1 
                else:
                    gene_counts[gene_name] = 1  
    return {gene for gene, count in gene_counts.items() if count == 1}

# list to store genes that appear exactly once in each file
unique_genes_per_file = [extract_unique_genes(os.path.join(directory, f)) for f in fna_files]

#get genes that show up 
common_unique_genes = set.intersection(*unique_genes_per_file)

# Output the names of the common genes
for gene in common_unique_genes:
    print(gene)