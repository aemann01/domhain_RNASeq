import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#get ureB copy 1 and flip cords
ureABC = pd.read_table("ureABC.operon",delimiter='\t', header= 0)
#get opeon length
ureABC['length']=abs(ureABC['A1'] - ureABC['C2'])
#sort by length
ureABC.sort_values("length")
ureABC.median(numeric_only=True)
ureABC.mode(numeric_only=True)
#dist
ureABC.hist(column="length")
plt.savefig('ureABC.operon.png')