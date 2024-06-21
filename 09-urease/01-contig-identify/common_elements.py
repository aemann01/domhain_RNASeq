import pandas as pd
df1 = pd.read_table("ureA.seqs",delimiter='\t')
df2 = pd.read_table("ureB.seqs",delimiter='\t')
df3 = pd.read_table("ureC.seqs",delimiter='\t')

df1.columns =['seq1']
df2.columns =['seq2']
df3.columns =['seq3']

list1 = df1['seq1'].tolist()
list2 = df2['seq2'].tolist()
list3 = df3['seq3'].tolist()

def find_common(list1, list2, list3):
    common = set()
    for elem in list1:
        if elem in list2 and elem in list3:
            common.add(elem)
    return common
 
common = find_common(list1, list2, list3)

seqs = pd.DataFrame(common)

seqs.to_csv('ureABC.seqs', sep="\t", index=False, header = False)
