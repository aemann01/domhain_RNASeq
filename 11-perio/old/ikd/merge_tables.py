import pandas as pd

# read in files
read_counts = pd.read_csv("~/rna_dohmain/11-perio/06-red-complex/red_counts.txt", delimiter="\t")
# clean up empty column
read_counts = read_counts.dropna(axis=1)
unprot = pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/locus2go.txt", delimiter="\t")
# uni2ref = pd.read_csv("gene_refseq_uniprotkb_collab", delimiter="\t")
# gene2acc = pd.read_csv("gene_protein_acc", delimiter="\t")
# kegg = pd.read_csv("ALL_genomes_KO.txt", delimiter="\t")
goterms = pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/gene2go.sub", delimiter="\t")

# merge files
unprot['Locus_tag'] = 'gene-' + unprot['Locus_tag']
read_counts['Geneid_base'] = read_counts['Geneid'].apply(lambda x: '_'.join(x.split('_')[:2]))


merged = pd.merge(unprot, read_counts, right_on="Geneid_base", left_on="Locus_tag")
# Set indexes before merging

merged2 = pd.merge(goterms, merged, right_on="Go_term", left_on="GO_ID")

# write merged data to file
with open("read_counts_goKegg.txt", "w") as outfile:
	merged2.to_csv(outfile, sep="\t", index=False)



import pandas as pd

read_counts = pd.read_csv("~/rna_dohmain/11-perio/06-red-complex/red_counts.txt", delimiter="\t")
unprot = pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/locus2go.txt", delimiter="\t")
goterms = pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/gene2go.sub", delimiter="\t")


# Set indices on merge columns
unprot.set_index('Locus_tag', inplace=True)

read_counts.set_index('Geneid_base', inplace=True)
goterms.set_index('GO_ID', inplace=True)

# Merge using indices
merged = unprot.join(read_counts, how='inner', on='Locus_tag')
merged2 = merged.join(goterms, how='inner', on='Go_term')

# Reset index if needed
merged2.reset_index(drop=True, inplace=True)