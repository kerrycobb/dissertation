import pandas as pd
from functools import reduce

df1 = pd.read_csv("qmat-clust-90-indel-16-samples-139-include-amer-terr-2-K-2.csv")
df2 = pd.read_csv("qmat-clust-90-indel-16-samples-149-include-amer-terr-1-K-2.csv")
df3 = pd.read_csv("qmat-clust-90-indel-16-samples-179-include-americanus-group-2-K-3.csv")

df1.columns = ["id", "at1-clust1", "at1-clust2"] 
df2.columns = ["id", "at2-clust1", "at2-clust2"] 
df3.columns = ["id", "ag-clust1", "ag-clust2", "ag-clust3"] 

merged = reduce(lambda left, right: pd.merge(left, right, on="id", how="outer"), 
    [df1, df2, df3]) 

merged.to_csv("merged-qmat.csv", index=False)
