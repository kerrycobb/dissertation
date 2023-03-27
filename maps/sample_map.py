
import pandas as pd
import matplotlib.pyplot as plt
from cartopy import crs
# from cartopy.feature import ShapelyFeature, NaturalEarthFeature
import cartopy.feature as cfeature
import math
import seaborn as sns
import numpy as np

def rgb(r, g, b):
    return (r/255, g/255, b/255) 


toadData = pd.read_csv("~/Desktop/toad-project/toad-data.csv")
toadData = toadData[toadData["species"] != "?"]
toadData = toadData[toadData["species"] != "fowleri / woodhousii intergrade"]
toadData["species"] = toadData["species"].astype("category")
colorMap = dict(zip(toadData["species"].cat.categories.tolist(), plt.cm.tab20.colors))

# Read in data
# weatherby_samples = pd.read_csv("weatherby-sites.csv")
# my_samples = pd.read_csv("cobb-sites.csv")
# fall_line = pd.read_csv("fall-line.csv")

# Filter data
# terr = my_samples[my_samples["species"] == "terrestris"]
# amer = my_samples[my_samples["species"] == "americanus"]


# # Getting Mid point
# # Work in progress, net yet working
# def degToEquirectangular(lon, lat):
#     # From wikipedia https://en.wikipedia.org/wiki/Equirectangular_projection
#     print(lon)
#     print(lat)
#     return 6371 * lon * np.cos(lat)
# def equirectangularToDeg(x, lat):
#     return x / (6371 * np.cos(lat)) 
# # Get midpoint
# def mid(mn, mx):
#     return mx - (mx - mn) / 2
# midLat = mid(toadData["latitude"].min(), toadData["latitude"].max())
# print(midLat)
# x = degToEquirectangular(toadData["longitude"], toadData["latitude"])
# mid_x = mid(x.min(), x.max())
# midLon = equirectangularToDeg(mid_x, midLat)
# print(midLon)
# print(midLat)

res = "50m"
lineWidth = .5
waterFaceColor = rgb(197, 213, 227)
waterEdgeColor = rgb(88, 135, 155)
borderColor = (0.6, 0.6, 0.6)
extent = [-120.38896371457172, -67.05643017812572, 23.33769204708392, 49.75492690666086]

proj = crs.Orthographic(central_longitude=-98.5795, central_latitude=36)
transform = crs.Geodetic()


# All species map
plt.figure()
ax = plt.axes(projection=proj)
ax.add_feature(cfeature.STATES.with_scale(res), edgecolor=borderColor, linewidth=lineWidth, zorder=1)
ax.add_feature(cfeature.OCEAN.with_scale(res), facecolor=waterFaceColor, edgecolor=waterEdgeColor, linewidth=lineWidth, zorder=2)
ax.add_feature(cfeature.LAKES.with_scale(res), facecolor=waterFaceColor, edgecolor=waterEdgeColor, linewidth=lineWidth/2, zorder=2)
for name, group in toadData.groupby("species"):
    ax.scatter(group["longitude"], group["latitude"], label=name, transform=transform, color=colorMap[name], edgecolor="black", linewidth=0.3, s=15, zorder=100)
# ax.scatter(midLon, midLat, s=100, color="black")

# ax.gridlines(linewidth=0.25, color="gray", linestyle=(0, (10, 10)), alpha=1.0)
ax.set_extent(extent)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.title("Sampling Distribution")
plt.savefig('sample-map.svg', transparent=True, bbox_inches='tight')


# Single species maps
for name, group in toadData.groupby("species"):
    plt.figure()
    ax = plt.axes(projection=proj)
    ax.add_feature(cfeature.STATES.with_scale(res), edgecolor=borderColor, linewidth=lineWidth)
    ax.add_feature(cfeature.OCEAN.with_scale(res), facecolor=waterFaceColor, edgecolor=waterEdgeColor, linewidth=lineWidth)
    ax.add_feature(cfeature.LAKES.with_scale(res), facecolor=waterFaceColor, edgecolor=waterEdgeColor, linewidth=lineWidth)
    ax.scatter(group["longitude"], group["latitude"], transform=transform, color=colorMap[name], edgecolor="black", linewidth=0.3, s=20, zorder=100)
    ax.gridlines(linewidth=0.5, color="gray", linestyle='--', alpha=0.5)
    ax.set_extent(extent)
    plt.title("${}$ Sampling Distribution".format(name))
    plt.savefig('species-sampling/{}.svg'.format(name), transparent=True, bbox_inches='tight')

