import pandas as pd
import numpy
import matplotlib.pyplot as plt

#get ureA and flip cords
ureA = pd.read_table("ureA.cords.filt",delimiter='\t', header= None)
ureA.columns =['seq', 'idA', 'A1', 'A2', 'strandA']
s = ureA['A1'] > ureA['A2']

#check if any need to be changed
test = pd.DataFrame(s,columns =['Name'])
test.drop(test[test.Name == False].index)
ureA.loc[s, ['A1','A2']] = ureA.loc[s, ['A2','A1']].values

#get ureB and flip cords
ureB = pd.read_table("ureB.cords.filt",delimiter='\t', header= None)
ureB.columns =['seq', 'idB', 'B1', 'B2', 'strandB']
s = ureB['B1'] > ureB['B2']

#check if any need to be changed
test = pd.DataFrame(s,columns =['Name'])
test.drop(test[test.Name == False].index)
ureB.loc[s, ['B1','B2']] = ureB.loc[s, ['B2','B1']].values

#merge by seq
ureAB = pd.merge(ureA, ureB, on ='seq', how='left', indicator=False)
ureAB
ureAB = ureAB.dropna()
#get ones that are on the same strand
strandAB = ureAB[ureAB[['strandA','strandB']].nunique(axis=1) == 1]

#postive strand
pos =strandAB.drop(strandAB[strandAB.strandB == "-"].index)
#get some stats for postive strand
pos['btwn'] = pos.B1 - pos.A2
pos.sort_values("btwn")
# pos['ureA']=pos.A2 - pos.A1
# pos['ureB']=pos.B2 - pos.B1
# pos.sort_values("ureA")
# pos.sort_values("ureB")
#make distro plot
pos.hist(column="btwn")
plt.savefig('ureAB.pos.png')
pos.median(numeric_only=True)
pos.mean(numeric_only=True)
pos.mode(numeric_only=True)
df =pos.drop(pos[pos.btwn > 2000].index)
df2 =df.drop(df[df.btwn < -50].index)
df2.to_csv('ureAB.pos.seqsid', sep="\t", index=False, header = True)
df2.hist(column="btwn")
plt.savefig('ureAB.pos.png')

#negative strand
#select which ones belong in a range between -2000 to 2000 basepairs of each other
neg =strandAB.drop(strandAB[strandAB.strandB == "+"].index)
#test for if ureA starts before B
neg['btwn'] = neg.B1 - neg.A2
s = neg['B1'] < neg['A1']
test = pd.DataFrame(s,columns =['Name'])
neg['bigger'] = test.drop(test[test.Name == True].index)
#calculate distance from ureB to ureA (in that order)
neg['btwn'] = neg.A1 - neg.B2

# pos['ureA']=pos.A2 - pos.A1
# pos['ureB']=pos.B2 - pos.B1
# pos.sort_values("ureA")
# pos.sort_values("ureB")
#make distro plot
neg.hist(column="btwn")
plt.savefig('ureAB.neg.png')
neg.median(numeric_only=True)
neg.mean(numeric_only=True)
neg.mode(numeric_only=True)
#drop stuff that falls out of this range
df =neg.drop(neg[neg.btwn > 2000].index)
df2 =df.drop(df[df.btwn < -50].index)
#make files
df2.to_csv('ureAB.neg.seqsid', sep="\t", index=False, header = True)
df2.hist(column="btwn")
plt.savefig('ureAB.neg.png')


