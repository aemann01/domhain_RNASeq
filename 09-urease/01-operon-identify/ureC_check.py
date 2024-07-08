import pandas as pd
import numpy as np
from Bio import SeqIO

#get ureC copy 1 and flip cords
ureC = pd.read_table("ureC.cords1",delimiter='\t', header= None)

ureC.columns =['seq', 'idC', 'C1', 'C2', 'strandC']

s = ureC['C1'] > ureC['C2']
ureC.loc[s, ['C1','C2']] = ureC.loc[s, ['C2','C1']].values

#get ureC copy 2 and flip cords

urec = pd.read_table("ureC.cords2",delimiter='\t', header= None)

urec.columns =['seq', 'idc', 'c1', 'c2', 'strandc']

s = urec['c1'] > urec['c2']
urec.loc[s, ['c1','c2']] = urec.loc[s, ['c2','c1']].values

#merge by seq
ureCs = pd.merge(ureC, urec, on ='seq', how='left', indicator=False)
ureCs
#get ones that are on the same strand
ureCs = ureCs[ureCs[['strandC','strandc']].nunique(axis=1) == 1]
#swtich so first coordinate is smaller
s = ureCs['C1'] > ureCs['c1']
ureCs.loc[s, ['C1','c1']] = ureCs.loc[s, ['c1','C1']].values
ureCs.loc[s, ['C2','c2']] = ureCs.loc[s, ['c2','C2']].values
ureCs.loc[s, ['idC','idc']] = ureCs.loc[s, ['idc','idC']].values

#see differences between the end of the first gene cord and start of second gene cord
ureCs['btwn1'] = ureCs.c1 - ureCs.C2
#is the strarting of the bigger coordinate between the first annotation or if the ending of the first coordinate is less than 300 basepairs from the beginning of the second cord
ureCs['inbtwn'] = ureCs.apply(lambda x: "True" if x['c1'] <=
                     x['C2'] and x['c1']
                     >= x['C1'] or x['btwn1'] <= 300
                     else "False", axis=1)
#flip end cords to make sure the largest is at the end
s = ureCs['c2'] < ureCs['C2']
ureCs.loc[s, ['c2','C2']] = ureCs.loc[s, ['C2','c2']].values

#make file
ureCs.to_csv('ureC.comb', sep="\t", index=False, header = True)