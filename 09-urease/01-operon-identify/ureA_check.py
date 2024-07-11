import pandas as pd
import numpy as np
from Bio import SeqIO

#get ureA copy 1 and flip cords
ureA = pd.read_table("ureA.cords1",delimiter='\t', header= None)

ureA.columns =['seq', 'idA', 'A1', 'A2', 'strandA']

s = ureA['A1'] > ureA['A2']
ureA.loc[s, ['A1','A2']] = ureA.loc[s, ['A2','A1']].values

#get ureA copy 2 and flip cords

urea = pd.read_table("ureA.cords2",delimiter='\t', header= None)

urea.columns =['seq', 'ida', 'a1', 'a2', 'stranda']

s = urea['a1'] > urea['a2']
urea.loc[s, ['a1','a2']] = urea.loc[s, ['a2','a1']].values

#merge by seq
ureAs = pd.merge(ureA, urea, on ='seq', how='left', indicator=False)
ureAs
#get ones that are on the same strand
ureAs = ureAs[ureAs[['strandA','stranda']].nunique(axis=1) == 1]
#swtich so first coordinate is smaller
s = ureAs['A1'] > ureAs['a1']
ureAs.loc[s, ['A1','a1']] = ureAs.loc[s, ['a1','A1']].values
ureAs.loc[s, ['A2','a2']] = ureAs.loc[s, ['a2','A2']].values
ureAs.loc[s, ['idA','ida']] = ureAs.loc[s, ['ida','idA']].values

#see differences between the end of the first gene cord and start of second gene cord
ureAs['btwn1'] = ureAs.a1 - ureAs.A1
#is the strarting of the bigger coordinate between the first annotation or if the ending of the first coordinate is less than 300 basepairs from the beginning of the second cord
ureAs['inbtwn'] = ureAs.apply(lambda x: "True" if x['a1'] <=
                     x['A2'] and x['a1']
                     >= x['A1'] or x['btwn1'] <= 300
                     else "False", axis=1)
#flip end cords to make sure the largest is at the end
s = ureAs['a2'] < ureAs['A2']
ureAs.loc[s, ['a2','A2']] = ureAs.loc[s, ['A2','a2']].values

#make file
ureAs.to_csv('ureA.comb', sep="\t", index=False, header = True)