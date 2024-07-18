import pandas as pd
df1 = pd.read_table("acc_taxids.tsv",delimiter='\t', header= None)
df1.columns =['seq', 'GCA', 'taxid']
print(df1.head())

df2 = pd.read_table("locus_tags",delimiter='\t', header= None)
df2.columns =['seq', 'locus_tag']
print(df2.head())

combined = pd.merge(df1, df2, on ='seq')


combined.to_csv('prep_headers.tsv', sep="\t", index=False, header=None)