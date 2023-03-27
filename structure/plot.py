#!/usr/bin/env python

import fire
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import math

colors = ["rgb(215,25,28)", "rgb(253,174,97)", "rgb(171,221,164)", "rgb(43,131,186)", "rgb(140, 86, 75)"]

class Structure:
    def __init__(self, structFile, dataFile, K, labels=None, extent=None):
        self.K = K
        self.labels = labels
        self.extent = extent
        self.pops = [i for i in range(1, self.K+1)]

        fh = open(structFile, "r")
        while True:
            line = fh.readline()
            if not line:
                quit("Error: no \"Inferred ancestry of individuals\" block")
            elif line.strip() == "Inferred ancestry of individuals:": 
                fh.readline()
                break
        rows = []
        while True:
            line = fh.readline().strip()
            if len(line) == 0:
                break
            s = line.split()
            row = [s[1]]
            for i in range(K):
                row.append(float(s[4+i]))
            rows.append(row)
        fh.close()
        structDf =  pd.DataFrame(rows, columns=["id"] + [i for i in range(1, K+1)]) 
        sampleDf = pd.read_csv(dataFile)        
        self.df = structDf.merge(sampleDf, how="left", left_on="id", right_on="sample_id2") 
        # self.df = structDf

    def adjustProps(self):
        ## Adjust proportions so they sum to 1.0
        self.df[self.K] = self.df[self.K] + 1.0 - self.df[self.pops].sum(axis=1)

    def groupSimiliar(self):
        ## Group populations by cluster similarity
        self.df["pop"] = self.df[self.pops].idxmax(axis=1)
        ascend = True
        groups = []
        for i in self.pops:
            subset = list(self.pops) 
            subset.remove(i)
            newDf = self.df[self.df["pop"] == i].copy()
            newDf.sort_values(by=subset, ascending=ascend, inplace=True)
            groups.append(newDf)
            ascend = not ascend
        self.df = pd.concat(groups) 
        self.df.reset_index(inplace=True)

    def barPlot(self):
        if self.labels:
            labels = self.labels 
        else:
            labels = self.pops
        data = []
        for i in self.pops:
            data.append(
                go.Bar(
                    name=labels[i-1], 
                    x=self.df.index, 
                    y=self.df[i],
                    text=self.df["id"]
                )
            )
        fig = go.Figure(data=data)
        fig.update_layout(
            barmode='stack',
            bargap=0.05,
            title="Structure",
            title_x=0.5,
            xaxis_title="Samples",
            yaxis_title="Ancestry Proportions",
        )
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(range=[0, 1]) 
        fig.show()

    def piePlot(self, size=1, extent=None):
        ## Extent: [left lower long, left lower lat, right upper long, right upper lat]
        ## Size: Radius of pie graph in degrees 
        assert size > 0.0

        if self.labels:
            labels = self.labels 
        else:
            labels = self.pops

        data = []
        x = [[] for _ in self.pops]
        y = [[] for _ in self.pops]
        text = []
        for _, row in self.df.iterrows():
            cum = 0
            text.append(row["rename"])
            for j in self.pops:
                rads = np.linspace(2 * np.pi, cum * 2 * np.pi, 100) 
                x0 = row["longitude"]
                y0 = row["latitude"]
                data.append(go.Scattergeo(
                    lon=[x0] + (size * np.cos(rads) + x0).tolist() + [x0],
                    lat=[y0] + (size * np.sin(rads) + y0).tolist() + [y0],
                    name=labels[j-1], 
                    text=row["rename"],
                    fillcolor=colors[j-1],
                    line_width=0,
                    fill="toself",
                    mode="lines",
                ))
                cum += row[j] 

        data.append(go.Scattergeo())
        fig = go.Figure(data=data)

        fig.update_geos(
            projection_type="mercator",
            visible=True, 
            resolution=50, 
            # scope="world",
            scope="north america",
            showcountries=True, 
        #     countrycolor="Black",
            showsubunits=True, 
        #     subunitcolor="grey",
        )
        if extent == None:
            fig.update_geos(fitbounds='locations')
        else:
            fig.update_geos(
                lonaxis_range=[extent[0], extent[2]],
                lataxis_range=[extent[1], extent[3]]
            )
        fig.show()

dataPath = "../toad-data.csv"
ext1 = [-120, 28, -65, 50]


## Americanus group 1 75%
# structPath = "clust-90-indel-16-samples-165-include-americanus-group-1-K-3/clust-90-indel-16-samples-165-include-americanus-group-1.out_f"
# S = Structure(structPath, dataPath, 3)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

# structPath = "clust-90-indel-16-samples-165-include-americanus-group-1-K-4/clust-90-indel-16-samples-165-include-americanus-group-1.out_f"
# S = Structure(structPath, dataPath, 4)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

# structPath = "clust-90-indel-16-samples-165-include-americanus-group-1-K-5/clust-90-indel-16-samples-165-include-americanus-group-1.out_f"
# S = Structure(structPath, dataPath, 5)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])


## Americanus group 1 90%
# structPath = "clust-90-indel-16-samples-198-include-americanus-group-1-K-3/clust-90-indel-16-samples-198-include-americanus-group-1-K-3.out_f"
# S = Structure(structPath, dataPath, 3)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

# structPath = "clust-90-indel-16-samples-198-include-americanus-group-1-K-4/clust-90-indel-16-samples-198-include-americanus-group-1-K-4.out_f"
# S = Structure(structPath, dataPath, 4)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

# structPath = "clust-90-indel-16-samples-198-include-americanus-group-1-K-5/clust-90-indel-16-samples-198-include-americanus-group-1-K-5.out_f"
# S = Structure(structPath, dataPath, 5)
# S.adjustProps()
# S.groupSimiliar()
# S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])


## Americanus group 2 - 90%
structPath = "clust-90-indel-16-samples-179-include-americanus-group-2-K-3/clust-90-indel-16-samples-179-include-americanus-group-2-K-3.out_f"
S = Structure(structPath, dataPath, 3)
S.adjustProps()
S.groupSimiliar()
S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

structPath = "clust-90-indel-16-samples-179-include-americanus-group-2-K-4/clust-90-indel-16-samples-179-include-americanus-group-2-K-4.out_f"
S = Structure(structPath, dataPath, 4)
S.adjustProps()
S.groupSimiliar()
S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])

structPath = "clust-90-indel-16-samples-179-include-americanus-group-2-K-5/clust-90-indel-16-samples-179-include-americanus-group-2-K-5.out_f"
S = Structure(structPath, dataPath, 5)
S.adjustProps()
S.groupSimiliar()
S.barPlot()
# # S.piePlot(size=0.5, extent=ext1)
# # S.piePlot(size=0.01, extent=[-88.5, 31.5, -85, 34])


## Amer-terr 75%