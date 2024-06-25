import pandas as pd
df = pd.read_table("ids.cords",delimiter='\t', header= None)
df.columns =['seq', 'C1', 'C2', 'full']
s = df['C1'] > df['C2']
df.loc[s, ['C1','C2']] = df.loc[s, ['C2','C1']].values
df.to_csv('cord.ids', sep="\t", index=False, header = False)
