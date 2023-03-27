# import plotly.plotly as py
# from plotly.graph_objs import *
import plotly.graph_objects as go
import numpy as np
# from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
# init_notebook_mode(connected=True)

lon_center = 45
lat_center = 0
radius = 5

thetas = np.linspace(0, 2*np.pi, 50)
lon = radius * np.cos(thetas) + lon_center
lat = radius * np.sin(thetas) + lat_center

trace1 = go.Scattergeo(
    lat=lat,
    line=go.Line(
        color='rgb(213,62,79)',
        width=3,
    ),
    fill="toself",
    lon=lon,
    mode='lines')

data = go.Data([trace1])
layout = go.Layout(
    geo=dict(
        lataxis=dict(
            gridcolor='rgb(102, 102, 102)',
            gridwidth=0.5,
            showgrid=True
        ),
        lonaxis=dict(
            gridcolor='rgb(102, 102, 102)',
            gridwidth=0.5,
            showgrid=True
        ),
        projection=dict(
            rotation=dict(
                lat=-10,
                lon=0,
                roll=-103
            ),
            type='orthographic'
        ),
        coastlinewidth=0
    ),
    showlegend=False,
)
fig = go.Figure(data=data, layout=layout)
fig.show()
# iplot(fig)