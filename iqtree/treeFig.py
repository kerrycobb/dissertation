from ete3 import Tree, TreeStyle, NodeStyle
import ete3 as et
import pandas as pd

df = pd.read_csv("../structure/qmat-clust-90-indel-16-samples-198-include-americanus-group-1-K-3.csv")
# df = pd.read_csv("../structure/qmat-clust-90-indel-16-samples-179-include-americanus-group-2-K-3.csv")
df.rename(columns={df.columns[0]: "id"}, inplace=True)

t = et.Tree("clust-80-indel-16-snps-0p3-samples-140-exclude-2.phy.contree")
t.set_outgroup("hera10484_unknown")
t.ladderize(direction=0)

def layout(node):
    if node.is_leaf():
        node.set_style(NodeStyle(dict(size=0)))
        search = df[df["id"] == node.name]
        if len(search) == 1:
            row = search.iloc[0]
            if row["Cluster1"] > 0.9:
                node.set_style(NodeStyle(dict(fgcolor="blue")))
            elif row["Cluster2"] > 0.9:
                node.set_style(NodeStyle(dict(fgcolor="red")))
            elif row["Cluster3"] > 0.9: 
                node.set_style(NodeStyle(dict(fgcolor="green")))
            else:
                node.set_style(NodeStyle(dict(fgcolor="pink")))
    else:
        node.set_style(NodeStyle(dict(size=0)))


ts = et.TreeStyle()
# ts.branch_vertical_margin = 1
ts.layout_fn = layout
ts.scale = 10000 


# ns = et.NodeStyle()
# ns["size"] = 0
# t.set_style(ns)
# for n in t.traverse():
    # n.set_style(ns)

# for n in t.iter_leaves():
#     search = df[df["id"] == n.name]
#     if len(search) == 1:
#         row = search.iloc[0]
#         if row["Cluster1"] > 0.95:
#             n.set_style(NodeStyle(dict(fgcolor="blue")))
#         if row["Cluster2"] > 0.95:
#             n.set_style(NodeStyle(dict(fgcolor="")))
#         if row["Cluster3"] > 0.95: 
#             n.set_style(NodeStyle(dict(fgcolor="green")))



#     n.set_style(ns)

t.show(tree_style=ts)