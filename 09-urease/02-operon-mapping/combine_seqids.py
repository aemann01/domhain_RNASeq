import pandas as pd
import numpy as np
df1 = pd.read_table("../01-operon-identify/ureABC.operon",delimiter='\t', header =0)

df2 = pd.read_table("full.id",delimiter='\t', header=None)
df2.columns =['seq', 'full', 'A1']
print(df2.head())

df_all = pd.merge(df1, df2, on =['seq','A1'], how='left', indicator=False)

df_all.to_csv('ids.cords', sep="\t", index=False, header =True)

#select the cords
df3 = df_all[["full","A1", "A2","B1","B2", "C1", "C2"]]
df3=df3.set_index('full')
df3.values.sort()
df3.to_csv('sorted.cords', sep="\t", index=True, header =True)
