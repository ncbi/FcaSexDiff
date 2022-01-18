#!/usr/bin/env python3
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd

infile = snakemake.input[0]
outfile = snakemake.output[0]
fcaver = snakemake.wildcards.fcaver
tissue = snakemake.wildcards.tissue

nodes = pd.read_hdf(infile, key='node')
links = pd.read_hdf(infile, key='link')
legend = pd.read_hdf(infile, key='legend')

def factorize(s):
    a = pd.factorize(s, sort=True)[0]
    return (a + 0.01) / (max(a) + 0.1)

nodes = (
    nodes
    .assign(x = lambda df: df.i/df.i.max())
    .assign(y = lambda df: df.j/df.j.max())
    .assign(x = lambda df: df.x.apply(lambda t: min(max(t, 0.001), 0.999)))
    .assign(y = lambda df: df.y.apply(lambda t: min(max(t, 0.001), 0.999)))
    #.assign(x = lambda df: factorize(df.i)) #/df.i.max())
    #.assign(y = lambda df: factorize(df.j)) #/df.i.max())
)

print(nodes)

names = nodes.name.tolist()

fig = make_subplots(rows=1, cols=1, row_heights=[1])

# create dummy scatter plots at orgin in all color levels
# just to make sure that the legend appears in the diagram
for idx, row in legend.iterrows():
    color = row['color']
    lab = row['label']
    leg = go.Scatter(
        x = [0], y = [0], fill = "toself",
        fillcolor = color,
        line_color = color,
        name = lab,
        #visible = 'legendonly',
        marker = dict(
            size = 0,
            line = dict(width=0, color='DarkSlateGrey')
        ),
    )
    fig.add_trace(leg, row=1, col=1)


fig.add_trace(go.Sankey(
    arrangement = 'snap',
    valueformat = "d",
    domain = {
        'x': [0, 1],
        'y': [0, 1],
    },
    node = dict(
        pad = 10,
        thickness = 20,
        line = dict(color = "gray", width = 0.5),
        label = nodes['label'],
        color = nodes['color'],
        x = nodes['x'],
        y = nodes['y'],
        customdata = nodes['details'],
        hovertemplate = 'Node: %{label}<br>'+
            '#cells: %{value}<br>'+
            '%{customdata}<br>'+
            '<extra></extra>'
    ),
    link = dict(
        source = [names.index(x) for x in links['source']],
        target = [names.index(x) for x in links['target']],
        value = links['weight'],
        color = links['color'],
        hovertemplate = 'Source: %{source.label}<br>'+
            'Target: %{target.label}<br>'+
            '#cells: %{value}<br>'+
            '<extra></extra>',
    )
))

fig.update_yaxes(visible = False, showticklabels = False)
fig.update_xaxes(visible = False, showticklabels = False)

fig.update_layout(
    title_text = f'Sankey Diagram | FCA {fcaver} | {tissue}',
    title_x = 0.5,
    font_size = 15,
    autosize = True,
    width = 2000,
    height = 2000,
    legend_bgcolor = 'rgba(0,0,0,0)',
    template = 'simple_white',
    legend = dict(
        title = dict(text = 'Bias', side = 'left'),
        orientation = "h",
        x = 0,
        y = 1,
        xanchor = "left",
        yanchor = "bottom",
        valign = 'middle',
    )
)
#fig.show()
fig.write_html(outfile)
