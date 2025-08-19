import pandas as pd

# Load kegg_ids
kegg_ids = pd.read_csv('kegg_ids', sep='\t', header=None, names=['gene', 'protein_accession', 'unirpot', 'ko'])

# Load kegg2map
kegg2map = pd.read_csv('kegg2map', sep='\t', header=None, names=['ko', 'path'])

# Load pathway
pathway = pd.read_csv('pathway', sep='\t', header=None, names=['path', 'description'])

# read in files
read_counts = pd.read_csv("~/rna_dohmain/11-perio/06-red-complex/red_counts.txt", delimiter="\t")
# clean up empty column
read_counts = read_counts.dropna(axis=1)
# Remove everything before the colon in both columns of kegg2map
kegg2map['ko'] = kegg2map['ko'].str.split(':').str[1]
kegg2map['path'] = kegg2map['path'].str.split(':').str[1]

# Merge kegg_ids with kegg2map on 'ko'
merged = kegg_ids.merge(kegg2map, on='ko', how='left')

# Merge the result with pathway on 'path'
final_merged = merged.merge(pathway, on='path', how='left')

# Replace NaN in 'description' with 'no_pathway'
final_merged['description'].fillna('no_pathway', inplace=True)

read_counts['Geneid_base'] = read_counts['Geneid'].apply(lambda x: '_'.join(x.split('_')[:2]))

merged2 = pd.merge(final_merged, read_counts, right_on="Geneid_base", left_on="gene")
# Save the final merged result to a new file
merged2.to_csv('ko_read_counts.tsv', sep='\t', index=False)

print("Merging completed. Output saved to 'ko_read_counts.tsv'.")
