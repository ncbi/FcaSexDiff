#!/usr/bin/env python
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import anndata as ad


adata_file = snakemake.input[0]
mphate_file = snakemake.input[1]
out_file = snakemake.output[0]
alluvial_file = snakemake.output[1]
fcaver = snakemake.wildcards['fcaver']
tissue = snakemake.wildcards['tissue']

print(adata_file)
print(mphate_file)
print(out_file)
print(tissue)
print(fcaver)

adata = ad.read_h5ad(adata_file)

nlevel_use = 6

mp_clusters = (
    pd.read_table(mphate_file, index_col = 0)
    [[f"MP{i}" for i in range(nlevel_use)]]
)
print(mp_clusters)

meta = (
    adata.obs[['sex', 'annotation']]
    .assign(sex = lambda df: df.sex.astype(str))
    .assign(annotation = lambda df: df.annotation.astype(str))
    .merge(mp_clusters, how = "right", left_index = True, right_index = True)
    .query('sex in ["male", "female"]')
)
print(meta)

resolutions = (
  ['sex']
  + mp_clusters.columns.tolist()
  + ['annotation']
)

(
    meta[resolutions]
    .assign(freq = 1)
    .groupby(resolutions)
    .count()
    .to_csv(alluvial_file, sep='\t')
)


#thresholds = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
#thresholds = [-1, -0.5, -0.25, 0, 0.25, 0.5, 1]
#thresholds = [-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1]
#thresholds = [-1, -0.5, 0.5, 1]
#
#
#def get_level(x):
#    for i, th in enumerate(thresholds):
#        if x < th:
#            return i
#    return len(thresholds)
#
#
#def convert_to_rgb(val, minval, maxval):
#    # "colors" is a series of RGB colors delineating a series of
#    # adjacent linear color gradients between each pair.
#    colors = [(255, 0, 0), (250, 250, 250), (0, 0, 255)]  # [RED, WHITE, BLUE]
#    # Determine where the given value falls proportionality within
#    # the range from minval->maxval and scale that fractional value
#    # by the total number in the "colors" pallette.
#    i_f = float(val-minval) / float(maxval-minval) * (len(colors)-1)
#    # Determine the lower index of the pair of color indices this
#    # value corresponds and its fractional distance between the lower
#    # and the upper colors.
#    i, f = int(i_f // 1), i_f % 1  # Split into whole & fractional parts.
#    # Does it fall exactly on one of the color points?
#    if f < sys.float_info.epsilon:  # Smallest possible difference.
#        return colors[i]
#    else:  # Otherwise return a color within the range between them.
#        (r1, g1, b1), (r2, g2, b2) = colors[i], colors[i+1]
#        return int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))
#
#def get_color(x):
#    rgb = convert_to_rgb(x, minval=0, maxval=len(thresholds))
#    return f'rgb{rgb}'
#
#
#legend = []
#labs = ['-inf'] + [str(th) for th in thresholds] + ['inf']
#for i in range(len(thresholds)+1):
#    label = '({},{}]'.format(labs[i], labs[i+1])
#    if (i==0):
#        val = thresholds[i]-1
#    elif (i==len(thresholds)):
#        val = thresholds[i-1]+1
#    else:
#        val = (thresholds[i]+thresholds[i-1])/2
#    color = get_color(get_level(val))
#    legend.append((label, color))
#
#legend = pd.DataFrame(legend, columns=['label', 'color'])
#print(legend)
count_bias_colors = {
    "FemaleOnly" : "#A60000", #"red",
    "FemaleSignificant" : "#FF5233", #"red",
    "FemaleNonsignificant" : "#FFbbbc", # "#BC544B", #"pink",
    "Unbiased" : "gray",
    "MaleNonsignificant" : "#B2D8FF", #"#73C2FB", #maya
    "MaleSignificant" : "#023672", #"blue",
    "MaleOnly" : "#021732" #"blue"
}
legend = pd.DataFrame(dict(
    label = count_bias_colors.keys(),
    color = count_bias_colors.values(),
))


def order_resolutions(df):
    return df.merge(pd.DataFrame(
        enumerate(resolutions),
        columns=['i', 'resolution']
    ))

node = pd.concat([
    pd.DataFrame(dict(
        resolution = res,
        cluster = meta[res].unique()
    ))
    for res in resolutions
])
print(node)

def rearrange(x):
    o2n = {v:i for i,v in enumerate(sorted(x))}
    return [o2n[v] for v in x]


node = (
#    pd.read_table(countsfile)
#    .rename(columns={
#        'log2_male_to_female':'ratio',
#        'padj_binom':'padj'
#    })
    node
    .assign(color = lambda df: df.cluster.apply(
        lambda x: count_bias_colors["FemaleOnly"] if x == "female" else
                  count_bias_colors["MaleOnly"] if x == "male" else
                  "gray"
    ))
#    .assign(level = lambda df: df['ratio'].apply(get_level))
#    .assign(color = lambda df: df['level'].apply(get_color))
#    .assign(sig = lambda df: df['padj'].apply(pvalstr))
#    .assign(label = lambda df: df.apply(
#        lambda x: f'{x["cluster"]}{x["sig"]}' if x['resolution'] == 'annotation' else
#                  f'{x["cluster"]}|{x["major_annotation"]}{x["sig"]}',
#        axis=1
#    ))
    .assign(label = lambda df: df.cluster)
#    #.query('resolution in ["L0.4"]')
#    .pipe(print_pass)
#    .append(pd.DataFrame({
#        'resolution' : ['root', 'root'],
#        'cluster': ['female', 'male'],
#        'label': ['female', 'male'],
#        'color': ['red', 'blue'],
#        'ratio': ['0', 'inf'],
#        'padj': [0,0]
#    }), sort=False, ignore_index=True)
#    .sort_values(['resolution', 'padj', 'ratio'])
    .assign(j = lambda df: (
        df.groupby('resolution')['cluster'].transform(rearrange)
    ))
#    .pipe(print_pass)
#    # now assign x as the order of the resolution first sex, then annotation
#    # and finally all other Leiden levels
    .pipe(order_resolutions)
#    .pipe(print_pass)
    .assign(details = "")
#    .pipe(print_pass)
#    .assign(annotation = lambda df: df['annotation'].fillna('{}'))
#    .assign(details = lambda df: df.apply(lambda x: '<br>'.join([
#        "log2(M/F): {}".format(round(float(x['ratio']), 5)),
#        "annotations:" ] + [
#            f"&nbsp;&nbsp;{k}:{v}" for k,v in
#            literal_eval(x['annotation']).items()
#    ]), axis=1))
#    .pipe(print_pass)
)
#print(node[['resolution', 'cluster', 'padj', 'i', 'x', 'j', 'y']])
#print(node.cluster.unique())
#print(node.cluster.unique().shape)


link = pd.DataFrame(columns='source target color weight'.split())

print(meta.columns.tolist())

for col in range(len(resolutions)-1):
    select = resolutions[col:col+2]
    if not "sex" in select:
        select.append('sex')
    tdf = (
        meta[select]
        .assign(weight=1)
        #.pipe(print_pass)
        .groupby(select)
        .count()
        .reset_index()
        .assign(source = lambda df: df[select[0]])
        .assign(target = lambda df: df[select[1]])
    )
    print(tdf)
    col_dict = {'male': 'rgba(0,0,255,0.09)', 'female': 'rgba(255,0,0,0.09)'}
    tdf['color']  = tdf['sex'].apply(lambda x: col_dict[x])
    link = link.append(tdf['source target color weight'.split()])

print(link)

print(node)
node.to_excel("haha.xlsx")
node = (
    node["cluster label color i j details".split()]
    .rename({"cluster": "name"}, axis=1)
)
node.to_hdf(out_file, key='node')
link.to_hdf(out_file, key='link')
legend.to_hdf(out_file, key='legend')

