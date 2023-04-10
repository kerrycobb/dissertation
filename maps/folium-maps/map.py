#!/usr/bin/env python

import pandas as pd
import folium
from folium.plugins import MarkerCluster
import fire

colors = {
    "americanus": "red",
    "terrestris": "blue",
    "quercicus": "gray",
    "fowleri": "orange",
    "woodhousii": "gray",
    "speciosus": "beige",
    "punctatus": "green",
    "microscaphus": "darkblue",
    "cognatus": "black",
    "hemiophrys": "lightblue",
    "baxteri": "darkpurple",
    "debilis": "pink",
    "retiformis": "lightgreen",
    "fowleri_woodhousii_intergrade": "black"}

def filter(df, path, oper):
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("#") and not line.isspace():
                ids.append(line)
                if not line in df["sample_id2"].values:
                    quit("Error!: {} not a valid sampled id".format(line))
    series = pd.Series(ids)
    if oper == "in":
        return df[df["sample_id2"].isin(series)]
    elif oper == "not":
        return df[~df["sample_id2"].isin(series)]


def makeMap(output, data, include=None, exclude=None):
    df = pd.read_csv(data)
    df = df.dropna(subset=["latitude", "longitude"])
    if include and exclude:
        quit("Can't have both include and exclude")
    if include:
        df = filter(df, include, "in") 
    if exclude:
        df = filter(df, exclude, "not") 

    m = folium.Map(tiles="Stamen Terrain")
    tot = 0
    for name, group in df.groupby(by="species"):
        tot += len(group)
        print("{}: {}".format(name, len(group)))
        # feature_group = folium.FeatureGroup(name=name).add_to(m)
        feature_group = MarkerCluster(name=name).add_to(m)
        for ix, row  in group.iterrows():
            feature_group.add_child(folium.Marker(
                location=[row["latitude"], row["longitude"]],
                popup=folium.Popup(
                    "{}<br>{}".format(row["species"], row["sample_id"]),
                    max_width=200),
                icon=folium.Icon(color=colors[name])))

    m.fit_bounds(m.get_bounds())
    folium.LayerControl().add_to(m)
    m.save(output)
    print("Total: {}".format(tot))


if "__main__" == __name__:
    fire.Fire(makeMap)