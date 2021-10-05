import pandas as pd
import anndata as ad
import numpy as np
import plotly.graph_objects as go
from matplotlib import colors
from pathlib import Path

ColorsG = [
    (u'indigo', u'#4B0082'), (u'gold', u'#FFD700'), (u'firebrick', u'#B22222'),
    (u'indianred', u'#CD5C5C'), (u'yellow', u'#FFFF00'),
    (u'darkolivegreen', u'#556B2F'), (u'darkseagreen', u'#8FBC8F'),
    (u'slategrey', u'#708090'), (u'darkorange', u'#FF8C00'),
    (u'darkslategrey', u'#2F4F4F'), (u'mediumvioletred', u'#C71585'),
    (u'mediumorchid', u'#BA55D3'), (u'papayawhip', u'#FFEFD5'),
    (u'mediumslateblue', u'#7B68EE'), (u'black', u'#000000'),
    (u'springgreen', u'#00FF7F'), (u'orange', u'#FFA500'),
    (u'lightsalmon', u'#FFA07A'), (u'brown', u'#A52A2A'),
    (u'turquoise', u'#40E0D0'), (u'olivedrab', u'#6B8E23'),
    (u'lightcyan', u'#E0FFFF'), (u'cyan', u'#00FFFF'), (u'silver', u'#C0C0C0'),
    (u'skyblue', u'#87CEEB'), (u'gray', u'#808080'), (u'darkturquoise', u'#00CED1'),
    (u'goldenrod', u'#DAA520'), (u'darkgreen', u'#006400'), (u'darkviolet', u'#9400D3'),
    (u'darkgray', u'#A9A9A9'), (u'lightpink', u'#FFB6C1'), (u'teal', u'#008080'),
    (u'darkmagenta', u'#8B008B'), (u'lightgoldenrodyellow', u'#FAFAD2'),
    (u'lavender', u'#E6E6FA'), (u'yellowgreen', u'#9ACD32'), (u'thistle', u'#D8BFD8'),
    (u'violet', u'#EE82EE'), (u'navy', u'#000080'), (u'dimgrey', u'#696969'),
    (u'orchid', u'#DA70D6'), (u'blue', u'#0000FF'), (u'ghostwhite', u'#F8F8FF'),
    (u'honeydew', u'#F0FFF0'), (u'cornflowerblue', u'#6495ED'), (u'darkblue', u'#00008B'),
    (u'coral', u'#FF7F50'), (u'darkkhaki', u'#BDB76B'), (u'mediumpurple', u'#9370DB'),
    (u'cornsilk', u'#FFF8DC'), (u'red', u'#FF0000'), (u'bisque', u'#FFE4C4'),
    (u'slategray', u'#708090'), (u'darkcyan', u'#008B8B'), (u'khaki', u'#F0E68C'),
    (u'wheat', u'#F5DEB3'), (u'deepskyblue', u'#00BFFF'), (u'rebeccapurple', u'#663399'),
    (u'darkred', u'#8B0000'), (u'steelblue', u'#4682B4'), (u'aliceblue', u'#F0F8FF'),
    (u'lightslategrey', u'#778899'), (u'gainsboro', u'#DCDCDC'), (u'c', (0, 0.75, 0.75)),
    (u'mediumturquoise', u'#48D1CC'), (u'g', (0, 0.5, 0)), (u'k', (0, 0, 0)),
    (u'purple', u'#800080'), (u'lightgrey', u'#D3D3D3'), (u'burlywood', u'#DEB887'),
    (u'darksalmon', u'#E9967A'), (u'beige', u'#F5F5DC'), (u'w', (1, 1, 1)),
    (u'azure', u'#F0FFFF'), (u'lightsteelblue', u'#B0C4DE'), (u'oldlace', u'#FDF5E6'),
    (u'greenyellow', u'#ADFF2F'), (u'royalblue', u'#4169E1'), (u'lightseagreen', u'#20B2AA'),
    (u'sienna', u'#A0522D'), (u'lightcoral', u'#F08080'), (u'orangered', u'#FF4500'),
    (u'navajowhite', u'#FFDEAD'), (u'lime', u'#00FF00'), (u'palegreen', u'#98FB98'),
    (u'mistyrose', u'#FFE4E1'), (u'seashell', u'#FFF5EE'), (u'mediumspringgreen', u'#00FA9A'),
    (u'fuchsia', u'#FF00FF'), (u'chartreuse', u'#7FFF00'), (u'blanchedalmond', u'#FFEBCD'),
    (u'peru', u'#CD853F'), (u'aquamarine', u'#7FFFD4'), (u'white', u'#FFFFFF'),
    (u'darkslategray', u'#2F4F4F'), (u'lightgray', u'#D3D3D3'), (u'ivory', u'#FFFFF0'),
    (u'darkgoldenrod', u'#B8860B'), (u'lawngreen', u'#7CFC00'), (u'chocolate', u'#D2691E'),
    (u'crimson', u'#DC143C'), (u'forestgreen', u'#228B22'), (u'darkgrey', u'#A9A9A9'),
    (u'olive', u'#808000'), (u'mintcream', u'#F5FFFA'), (u'antiquewhite', u'#FAEBD7'),
    (u'b', (0, 0, 1)), (u'floralwhite', u'#FFFAF0'), (u'moccasin', u'#FFE4B5'),
    (u'limegreen', u'#32CD32'), (u'saddlebrown', u'#8B4513'), (u'grey', u'#808080'),
    (u'darkslateblue', u'#483D8B'), (u'lightskyblue', u'#87CEFA'), (u'r', (1, 0, 0)),
    (u'deeppink', u'#FF1493'), (u'plum', u'#DDA0DD'), (u'cadetblue', u'#5F9EA0'),
    (u'dodgerblue', u'#1E90FF'), (u'maroon', u'#800000'), (u'sandybrown', u'#F4A460'),
    (u'aqua', u'#00FFFF'), (u'magenta', u'#FF00FF'), (u'tan', u'#D2B48C'),
    (u'rosybrown', u'#BC8F8F'), (u'pink', u'#FFC0CB'), (u'lightblue', u'#ADD8E6'),
    (u'palevioletred', u'#DB7093'), (u'mediumseagreen', u'#3CB371'),
    (u'slateblue', u'#6A5ACD'), (u'dimgray', u'#696969'), (u'powderblue', u'#B0E0E6'),
    (u'seagreen', u'#2E8B57'), (u'snow', u'#FFFAFA'), (u'mediumblue', u'#0000CD'),
    (u'midnightblue', u'#191970'), (u'paleturquoise', u'#AFEEEE'),
    (u'palegoldenrod', u'#EEE8AA'), (u'whitesmoke', u'#F5F5F5'),
    (u'darkorchid', u'#9932CC'), (u'salmon', u'#FA8072'), (u'lightslategray', u'#778899'),
    (u'lemonchiffon', u'#FFFACD'), (u'lightgreen', u'#90EE90'), (u'tomato', u'#FF6347'),
    (u'hotpink', u'#FF69B4'), (u'lightyellow', u'#FFFFE0'), (u'lavenderblush', u'#FFF0F5'),
    (u'm', (0.75, 0, 0.75)), (u'linen', u'#FAF0E6'), (u'mediumaquamarine', u'#66CDAA'),
    (u'green', u'#008000'), (u'blueviolet', u'#8A2BE2'), (u'y', (0.75, 0.75, 0)),
    (u'peachpuff', u'#FFDAB9')
]


inpath = Path(snakemake.input[0])
outpath = Path(snakemake.output[0])


data = (
    ad.read_h5ad(snakemake.input[0])
    .obs[['sex', 'annotation', 'annotation_broad', 'scpopcorn_cluster']]
    .assign(sex = lambda df: df.sex.apply(lambda x: x[0].upper()))
    .assign(annotation = lambda df: df.apply(
        lambda x: x['sex'] + '_' + x['annotation'],
        axis = 1
    ))
    .assign(annotation_broad = lambda df: df.apply(
        lambda x: x['sex'] + '_' + x['annotation_broad'],
        axis = 1
    ))
    .reset_index()
)

print(data)


def prepare_sankey_data(df, label_key, cluster_key):
    edges = (
        df[['CellID', 'sex', label_key, cluster_key]]
        .rename(columns={label_key:'label', cluster_key:'cluster'})
        .assign(cluster = lambda df: df["cluster"].apply(
            lambda x: f'C{int(x)}'
        ))
        .groupby(["cluster", "label"])
        .agg({'CellID':'count', 'sex':'first'})
        .reset_index()
        .assign(ncell_cluster = lambda df: df.groupby("cluster").CellID.transform("sum"))
        .assign(ncell_annot = lambda df: df.groupby("label").CellID.transform("sum"))
        .assign(percent_annot = lambda df: df["CellID"]*100/df["ncell_annot"])
    )

    print(edges)

    return edges

def export_plot_data(edges, writer, sheet_name="AllEdges"):

    res = (
        edges[['cluster', 'ncell_cluster', 'label', 'ncell_annot', 'CellID', 'percent_annot']]
        .sort_values(by=['ncell_cluster', 'percent_annot', 'label'],
                     ascending=[False, False, True])
        .rename(columns={'cluster': 'scPopCornCluster',
                         'ncell_cluster': 'nCellInCluster',
                         'label': 'Annotation',
                         'CellID': 'nCellWithAnnotationInCluster',
                         'ncell_annot': 'nCellWithAnnotation',
                         'percent_annot': 'percentCellsWithAnnotationInCluster'})
    )
    print(res)

    res.to_excel(writer, sheet_name=sheet_name, index=False)




def plot_sankey(edges, title):

  females = (
          edges.query('sex == "F"')
          [["label"]]
          .drop_duplicates()
          .assign(i=0, color = "red")
  )

  males = (
          edges.query('sex == "M"')
          [["label"]]
          .drop_duplicates()
          .assign(i=2, color = "blue")
  )

  clusters = (
          edges[["cluster"]]
          .rename(columns={"cluster":'label'})
          .drop_duplicates()
          .assign(i=1, color = lambda df: df.label.apply(
              lambda x: f'rgba{colors.to_rgba(ColorsG[int(x[1:])+1][1])}'
           ))
  )

  nodes = (
      pd.concat([females, clusters, males], ignore_index = True, axis = 0)
      .assign(j = lambda df: df.groupby('i').label.transform(lambda s:range(s.shape[0])))
      .assign(x = lambda df: 0.001 + df.i/4)
      .assign(y = lambda df: (df.j + 1)/(1 + df.groupby('i').count().values.max()))
  )

  edges = (
      edges
      .assign(source = lambda tdf: [nodes.label.tolist().index(x) for x in tdf.cluster])
      .assign(target = lambda tdf: [nodes.label.tolist().index(x) for x in tdf.label])
  )

  fig = go.Figure(data=[go.Sankey(
      arrangement='snap',
      valueformat="d",
      node = dict(
          pad = 15,
          thickness = 20,
          line = dict(color = "black", width = 0.5),
          label = nodes.label,
          color = nodes.color,
          x = nodes.x,
          y = nodes.y
      ),
      link = dict(
        source = edges.source,
        target = edges.target,
        value = edges.CellID
    ))])

  fig.update_layout(
      title_text=title,
      font_size=20,
      autosize=False,
      width=2500,
      height=1500
  )
  #fig.show()
  fig.write_image(str(outpath.with_suffix('.png')))
  fig.write_html(str(outpath.with_suffix('.html')))


edges = prepare_sankey_data(data, 'annotation', 'scpopcorn_cluster')

writer = pd.ExcelWriter(str(outpath.with_suffix('.xlsx')), engine='xlsxwriter')

export_plot_data(edges, writer, 'AllEdges')

edges = edges.query("percent_annot > 5.0")

export_plot_data(edges, writer, "ThickEdges")

writer.save()

plot_sankey(edges, f"FCA {snakemake.wildcards.tissue} male/female integration using scPopCorn")


