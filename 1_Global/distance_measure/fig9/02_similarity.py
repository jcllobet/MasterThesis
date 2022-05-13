import pandas as pd
import backboning as bb

df = pd.read_csv("dist_table.csv", sep = "\t")
df = df.corr(method = "spearman").unstack().reset_index()
df.columns = ("src", "trg", "nij")
df = df[(df["src"] != df["trg"]) & (df["nij"] > 0)]

df_nc = bb.noise_corrected(df, undirected = True)
bb.test_densities(df_nc, 0.1, 0.15, 0.01) # Should find the highest theshold ensuring a single connected component, to properly evaluate all similarities
df_nc_bb = bb.thresholding(df_nc, 0.145)
df_nc_bb.to_csv("fig9.csv", sep = "\t", index = False)


