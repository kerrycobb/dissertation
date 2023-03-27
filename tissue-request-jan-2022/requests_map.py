
import pandas as pd
import folium
from folium.plugins import MarkerCluster

colors = ['purple', 'red', 'lightred', 'black', 'orange', 'lightgray', 
    'darkblue', 'lightgreen', 'cadetblue', 'darkred', 'blue', 'pink', 
    'gray', 'darkgreen', 'white', 'green', 'beige', 'lightblue', 'darkpurple']

data_path = "samples-for-request.csv"
df = pd.read_csv(data_path)
# df = df[(df.latitude != "?") & (df.longitude != "?")]
# df = df.dropna(subset=["latitude", "longitude"])
# df['latitude'] = df['latitude'].astype('float').round(6)
# df['longitude'] = df['longitude'].astype('float').round(6)


m = folium.Map(tiles="Stamen Terrain")

cix = 0
for name, group in df.groupby(by="species"):
    # feature_group = folium.FeatureGroup(name=name).add_to(m)
    feature_group = MarkerCluster(name=name).add_to(m)
    for ix, row  in group.iterrows():
        feature_group.add_child(folium.Marker(
            location=[row["latitude"], row["longitude"]],
            popup=folium.Popup(
                "{}<br>{}<br>{}".format(row["species"], row["year"], row["id"]),
                max_width=200),
            icon=folium.Icon(color=colors[cix]) 
            ))
    cix += 1

m.fit_bounds(m.get_bounds())
folium.LayerControl().add_to(m)
m.save("requests-map.html")