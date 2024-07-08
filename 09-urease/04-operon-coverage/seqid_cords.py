import pandas as pd
#H
df1 = pd.read_table("sig.H.depth",delimiter='\t')
df1.columns =['SEQ', 'cord', 'depth']

df2 = pd.read_table("line_cords",delimiter='\t')
df2.columns =['SEQ', 'A1', 'A2', 'B1', 'B2', 'C1', 'C2']

df_all = pd.merge(df1, df2, on ='SEQ', how='left', indicator=False)

df_all.to_csv('all_H.depth', sep="\t", index=False, header = False)

#D
df3 = pd.read_table("sig.D.depth",delimiter='\t')
df3.columns =['SEQ', 'cord', 'depth']

df_all = pd.merge(df3, df2, on ='SEQ', how='left', indicator=False)
df_all.to_csv('all_D.depth', sep="\t", index=False, header = False)

