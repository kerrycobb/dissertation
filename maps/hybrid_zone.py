
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

# Read in data
weatherby = pd.read_csv("weatherby-sites.csv")
weatherby = weatherby.drop_duplicates(['lon', 'lat'])
fallLine = pd.read_csv("fall-line.csv")

# Filter data
terr = toadData[toadData["species"] == "terrestris"]
amer = toadData[toadData["species"] == "americanus"]


res = "10m"
blue = rgb(27, 119, 178)
orange = rgb(255, 127, 30)
waterFaceColor = rgb(197, 213, 227)
waterEdgeColor = rgb(88, 135, 155)
landColor = (0.99, 0.99, 0.99)
extent1 = [-88.6, -84.8, 30, 35.1]
extent2 = [-87.4, -85.0, 32.1, 33.9]

proj = crs.Mercator()
transform = crs.Geodetic()


# Setup Map
ax = plt.axes(projection=proj)
ax.set_extent(extent1)
# Names for NaturalEarthFeature found here: https://github.com/nvkelso/natural-earth-vector
ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_2_counties", "10m"), facecolor=landColor, linewidth=0.25, edgecolor=(0.6, 0.6, 0.6))
ax.add_feature(cfeature.STATES.with_scale(res), facecolor="none", linewidth=0.5)
ax.add_feature(cfeature.OCEAN.with_scale(res), facecolor=waterFaceColor)
ax.add_feature(cfeature.COASTLINE.with_scale(res), edgecolor=waterEdgeColor, linewidth=0.3)


# Add Cities
lonOffset = .05 
latOffset = 0 
textSize = 6 
citySize = 3 
ax_aub_t = ax.text(-85.45 + lonOffset, 32.57 + latOffset, "Auburn", ha='left', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_aub = ax.scatter(-85.45, 32.59, transform=transform, color="black", marker="*", s=citySize, zorder=1000)

ax_birm_t = ax.text(-86.8104 - lonOffset, 33.5186 + latOffset, "Birmingham", ha='right', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_birm = ax.scatter(-86.8104, 33.5186, transform=transform, color="black", marker="*", s=citySize, zorder=1000)

ax_mont_t = ax.text(-86.3077 - lonOffset, 32.3792 + latOffset, "Montgomery", ha='right', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_mont = ax.scatter(-86.3077, 32.3792, transform=transform, color="black", marker="*", s=citySize, zorder=1000)


# Weatherby
ax_weatherby = ax.scatter(weatherby["lon"], weatherby["lat"], transform=transform, marker="o", edgecolor="black", linewidth=0.3, s=7, zorder=100)
plt.title("Weatherby 1984, Morphological Intermediates")
plt.savefig('hybrid-zone-weatherby.svg', transparent=True, bbox_inches='tight')
ax_weatherby.remove()


# Sampled Data
ax_amer = ax.scatter(amer["longitude"], amer["latitude"], label="americanus", transform=transform, color=blue, edgecolor="black", linewidth=0.3, s=7, zorder=101)
ax_terr = ax.scatter(terr["longitude"], terr["latitude"], label="terrestris", transform=transform, color=orange, edgecolor="black", linewidth=0.3, s=7, zorder=100)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=6)
plt.title("Sampling Distribution")
plt.savefig("hybrid-zone-sample.svg", transparent=True, bbox_inches='tight')



##########
# Zoom 

# Setup Map
plt.figure()
ax = plt.axes(projection=proj)
ax.set_extent(extent2)
# Names for NaturalEarthFeature found here: https://github.com/nvkelso/natural-earth-vector
ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_2_counties", "10m"), facecolor="none", linewidth=0.25, edgecolor=(0.6, 0.6, 0.6))
ax.add_feature(cfeature.STATES.with_scale(res), facecolor="none", linewidth=0.5)


# Add Cities
lonOffset = .025 
latOffset = 0 
textSize = 9 
citySize = 10 
markerSize = 20 
ax_aub_t = ax.text(-85.4808 + lonOffset, 32.59 + latOffset, "Auburn", ha='left', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_aub = ax.scatter(-85.4808, 32.59, transform=transform, color="black", marker="*", s=citySize, zorder=1000)

ax_birm_t = ax.text(-86.8104 - lonOffset, 33.5186 + latOffset, "Birmingham", ha='right', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_birm = ax.scatter(-86.8104, 33.5186, transform=transform, color="black", marker="*", s=citySize, zorder=1000)

ax_mont_t = ax.text(-86.25 - lonOffset, 32.34 + latOffset, "Montgomery", ha='right', 
    va="center", transform=transform, fontsize=textSize, zorder=1000)
ax_mont = ax.scatter(-86.3077, 32.3792, transform=transform, color="black", marker="*", s=citySize, zorder=1000)


# Weatherby Zoom
ax_weatherby = ax.scatter(weatherby["lon"], weatherby["lat"], transform=transform, marker="o", edgecolor="black", linewidth=0.3, s=30, zorder=100)
plt.title("Weatherby 1984, Morphological Intermediates")
plt.savefig('hybrid-zone-weatherby-zoom.svg', transparent=True, bbox_inches='tight')
ax_weatherby.remove()

# Sampled Data Zoom
ax.add_feature(cfeature.LAKES.with_scale(res), facecolor=waterEdgeColor, edgecolor=waterEdgeColor)
ax.add_feature(cfeature.RIVERS.with_scale(res), edgecolor=waterEdgeColor)
ax_amer = ax.scatter(amer["longitude"], amer["latitude"], label="americanus", transform=transform, color=blue, edgecolor="black", linewidth=0.3, s=markerSize, zorder=101)
ax_terr = ax.scatter(terr["longitude"], terr["latitude"], label="terrestris", transform=transform, color=orange, edgecolor="black", linewidth=0.3, s=markerSize, zorder=100)
plt.legend(loc='upper left', borderaxespad=0)
plt.title("Sampling Distribution")
plt.savefig("hybrid-zone-sample-zoom.svg", transparent=True, bbox_inches='tight')

