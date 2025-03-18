import pandas as pd
df1 = pd.read_table("annotations.denovo.txt",delimiter='\t', header= 0)
print(df1.head())

df2 = pd.read_table("all.ids",delimiter='\t', header= None)
df2.columns =['denovo', 'locus_tag']
print(df2.head())

combined = pd.merge(df1, df2, on ='locus_tag')

df3 = pd.read_table("taxonomy.txt",delimiter='\t', header= None)
df3.columns =['denovo', 'kingdom','phylum','class','order','family','genus','species','strain']
print(df3.head())

combined2 = pd.merge(combined, df3, on ='denovo', how = 'left')

combined2.to_csv('annotations.denovo.txt', sep="\t", index=False, header=True)