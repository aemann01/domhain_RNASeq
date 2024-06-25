import pandas as pd
df = pd.read_table("operon.seqs",delimiter='\t', header= None)
df.columns =['seq', 'C1', 'C2']
s = df['C1'] > df['C2']
df.loc[s, ['C1','C2']] = df.loc[s, ['C2','C1']].values
df.to_csv('operon.seq', sep="\t", index=False, header = False)
