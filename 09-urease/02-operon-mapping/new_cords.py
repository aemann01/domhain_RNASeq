import pandas as pd
import numpy
#read file in
ureABC = pd.read_table("ids.cords",delimiter='\t', header= 0)
#find the smallest coordinate
df3 = ureABC[["full","A1", "A2","B1","B2", "C1", "C2"]]
df3=df3.set_index('full')
df3.values.sort()
subs = df3.A1 -1
subs = pd.to_numeric(subs)
subs=subs.reset_index(drop=True)
subs.columns =['rem']

#merge to the dataframe
df = pd.concat([ureABC, subs], axis=1)
df.columns =['seq','idA','idB','idC','btwnab','btwnbc',	'A1',	'A2',	'B1',	'B2',	'C1',	'C2',	'strand','full','rem']

df['newA1'] = df.A1 - df.rem
df['newA2'] = df.A2 - df.rem
df['newB1'] = df.B1 - df.rem
df['newB2'] = df.B2 - df.rem
df['newC1'] = df.C1 - df.rem
df['newC2'] = df.C2 - df.rem

df.to_csv('new.cords', sep="\t", index=False, header =True)
