import pandas as pd
import numpy as np
from Bio import SeqIO

#get ureB copy 1 and flip cords
ureB = pd.read_table("ureB.cords1",delimiter='\t', header= None)

ureB.columns =['seq', 'idB', 'B1', 'B2', 'strandB']

s = ureB['B1'] > ureB['B2']
ureB.loc[s, ['B1','B2']] = ureB.loc[s, ['B2','B1']].values

#get ureB copy 2 and flip cords

ureb = pd.read_table("ureB.cords2",delimiter='\t', header= None)

ureb.columns =['seq', 'idb', 'b1', 'b2', 'strandb']

s = ureb['b1'] > ureb['b2']
ureb.loc[s, ['b1','b2']] = ureb.loc[s, ['b2','b1']].values

#merge by seq
ureBs = pd.merge(ureB, ureb, on ='seq', how='left', indicator=False)
ureBs
#get ones that are on the same strand
ureBs = ureBs[ureBs[['strandB','strandb']].nunique(axis=1) == 1]
#swtich so first coordinate is smaller
s = ureBs['B1'] > ureBs['b1']
ureBs.loc[s, ['B1','b1']] = ureBs.loc[s, ['b1','B1']].values
ureBs.loc[s, ['B2','b2']] = ureBs.loc[s, ['b2','B2']].values
ureBs.loc[s, ['idB','idb']] = ureBs.loc[s, ['idb','idB']].values

#see differences between the end of the first gene cord and start of second gene cord
ureBs['btwn1'] = ureBs.b1 - ureBs.B2
#is the strarting of the bigger coordinate between the first annotation or if the ending of the first coordinate is less than 300 basepairs from the beginning of the second cord
ureBs['inbtwn'] = ureBs.apply(lambda x: "True" if x['b1'] <=
                     x['B2'] and x['b1']
                     >= x['B1'] or x['btwn1'] <= 300
                     else "False", axis=1)
#flip end cords to make sure the largest is at the end
s = ureBs['b2'] < ureBs['B2']
ureBs.loc[s, ['b2','B2']] = ureBs.loc[s, ['B2','b2']].values

#make file
ureBs.to_csv('ureB.comb', sep="\t", index=False, header = True)