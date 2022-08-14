import pandas as pd


df = pd.read_csv("../data/ScienceDirect_files_05May2022_18-41-51.867 2/TableS1.csv")

print(df)

print(df.columns)


df["Female_Fatbody_DamIDseq1"] = df["Female_Fatbody_DamIDseq_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)
df["Male_Fatbody_DamIDseq1"] = df["Male_Fatbody_DamIDseq_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)
df["Female_Fatbody_DamIDchip1"] = df["Female_Fatbody_DamID-chip_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)
df["Female_Ovary_DamIDseq1"] = df["Female_Ovary_DamIDseq_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)
df["DSXM_S2_ChIPseq1"] = df["DSXM_S2_ChIPseq_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)
df["DSXF_S2_ChIPseq1"] = df["DSXF_S2_ChIPseq_Gene_Level_Occupancy_Gene+1kb_Peaksum_Ranked_Score"].rank(method ="max", pct=True)

print(df)


df = df[(df["Female_Fatbody_DamIDseq1"] > 0.9) |
        (df["Male_Fatbody_DamIDseq1"] > 0.9) |
        (df["Female_Fatbody_DamIDchip1"] > 0.9) |
        (df["Female_Ovary_DamIDseq1"] > 0.9) |
        (df["DSXM_S2_ChIPseq1"] > 0.9) |
        (df["DSXF_S2_ChIPseq1"] > 0.9)]

print(df)

df = df[df["Gene_Flybase_ID"].str.startswith("FBgn")]
print(df)

df[["Gene_Flybase_ID", "Gene_Name"]].to_csv("DSX_targets.tsv", sep="\t", index=False)
