#!/usr/bin/python3

from Bio import SeqIO
import pandas as pd

# Path to your GenBank file
genbank_file = 'combined_genbank_file.gbk'

# List to store the extracted data
data = []

# Parse the GenBank file
for record in SeqIO.parse(genbank_file, 'genbank'):
    for feature in record.features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
            protein_id = feature.qualifiers.get('protein_id', [''])[0]
            data.append({
                'Locus_tag': locus_tag,
                'Protein_ID': protein_id
            })

# Convert list of dictionaries to a DataFrame
df = pd.DataFrame(data)
df_uniprot = pd.read_csv('gene_refseq_uniprotkb_collab', delimiter="\t")
merged_df = pd.merge(df, df_uniprot, left_on='Protein_ID', right_on='#NCBI_protein_accession', how='inner')
final_df = merged_df[['Locus_tag', 'UniProtKB_protein_accession']]
with open("locus2uniprot.txt", "w") as output:
    final_df.to_csv(output, sep="\t", index=False)

