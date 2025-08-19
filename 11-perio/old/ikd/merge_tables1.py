import pandas as pd

# read in files
read_counts = 
pd.read_csv("~/rna_dohmain/11-perio/06-red-complex/red_counts.txt", 
delimiter="\t")
# clean up empty column
read_counts = read_counts.dropna(axis=1)
unprot = 
pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/locus2go.txt", 
delimiter="\t")
# uni2ref = 
pd.read_csv("gene_refseq_uniprotkb_collab", 
delimiter="\t")
# gene2acc = pd.read_csv("gene_protein_acc", 
delimiter="\t")
# kegg = pd.read_csv("ALL_genomes_KO.txt", 
delimiter="\t")
goterms = 
pd.read_csv("~/rna_dohmain/11-perio/09-KEGG/gene2go.sub", 
delimiter="\t")

# merge files
unprot['Locus_tag'] = 'gene-' + unprot['Locus_tag']
read_counts['Geneid_base'] = 
read_counts['Geneid'].apply(lambda x: 
'_'.join(x.split('_')[:2]))


merged = pd.merge(unprot, read_counts, 
right_on="Geneid_base", left_on="Locus_tag")
# Set indexes before merging
with open("go_reads.txt", "w") as outfile:
	merged.to_csv(outfile, sep="\t", index=False)
