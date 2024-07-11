import pandas as pd
import numpy
import matplotlib.pyplot as plt


#get ureAB cords for neg
ureAB = pd.read_table("ureAB.neg.seqsid",delimiter='\t', header= 0)

#get ureB and flip cords
ureC = pd.read_table("ureC.cords.filt",delimiter='\t', header= None)

ureC.columns =['seq', 'idC', 'C1', 'C2', 'strandC']

s = ureC['C1'] > ureC['C2']
#check if any need to be changed
test = pd.DataFrame(s,columns =['Name'])
test.drop(test[test.Name == False].index)
ureC.loc[s, ['C1','C2']] = ureC.loc[s, ['C2','C1']].values

#merge by seq
ureABC = pd.merge(ureAB, ureC, on ='seq', how='left', indicator=False)
ureABC
#test for if ureB starts before ureC
ureABC['btwn'] = ureABC.C1 - ureABC.B2
s = ureABC['C1'] < ureABC['B1']
test = pd.DataFrame(s,columns =['Name'])
ureABC['bigger'] = test.drop(test[test.Name == True].index)
#get ones that are on the same strand
strandABC = ureABC[ureABC[['strandB','strandC']].nunique(axis=1) == 1]

strandABC['btwnbc'] = strandABC.B1 - strandABC.C2
strandABC.sort_values("btwnbc")
strandABC['btwn'] = strandABC.A1 - strandABC.B2
# strandABC['ureC']=strandABC.C2 - strandABC.C1
# strandABC.sort_values("ureC")

#select which ones belong in a range between -2000 to 2000 basepairs of each other
strandABC.median(numeric_only=True)
strandABC.mean(numeric_only=True)
strandABC.mode(numeric_only=True)
#check distro
strandABC.hist(column="btwnbc")
plt.savefig('ureBC.neg.png')

#make file and drop things that are too long
df =strandABC.drop(strandABC[strandABC.btwnbc > 1500].index)
df2 =df.drop(df[df.btwnbc < -50].index)
df2.to_csv('ureABC.neg.seqsid', sep="\t", index=False, header = True)
df2.sort_values("btwnbc")
#make distro plot
df2.hist(column="btwnbc")
plt.savefig('ureBC.neg.png')
