import pandas as pd
# in and out files
df = pd.read_table("gene_cluster.tab",delimiter='\t', header= 0, index_col =0)
result = (df == 1).all(axis=1)
reslts = pd.DataFrame(data=result)
reslts.columns =['core']
result2 = reslts[reslts["core"] == True]
rows = result2.index.tolist()

pd.DataFrame(rows).to_csv('single-copy-core-tags', header=False, index=False)
