import pandas as pd
df1 = pd.read_table("seqf2GCA_seqs",delimiter='\t', header= None)
df1.columns =['seq', 'GCA']
print(df1.head())

df2 = pd.read_table("GCA_2_taxid",delimiter='\t', header= None)
df2.columns =['GCA', 'TaxID']
print(df2.head())

combined = pd.merge(df1, df2, on ='GCA')

print(combined['GCA'].unique())

combined.to_csv('acc_taxids.tsv', sep="\t", index=False, header=None)