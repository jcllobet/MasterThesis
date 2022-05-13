import pandas as pd

df = pd.read_csv("year_country_product.csv", sep = "\t")
ps = pd.read_csv("PS_SITC_edges_4digit", sep = "\t", header = None, names = ("src", "trg"))

# Countries need to be present from first to last year
_ = df[["year", "exporter"]].drop_duplicates().groupby(by = "exporter").size()
countries = set(_[_ == _.max()].index)

# Products need to be present in the product space
products = set(ps["src"]) | set(ps["trg"])
df = df[df["exporter"].isin(countries) & df["sitc"].isin(products)]

# Let's make decades averages
df["decade"] = pd.cut(df["year"], 5, labels = ["1960s", "1970s", "1980s", "1990s", "2000s"])
df4 = df.groupby(by = ["decade", "exporter", "sitc"])["value"].mean().reset_index()
df4["value"] = df4["value"].fillna(0)
df4 = pd.pivot_table(df4, index = ["decade", "exporter"], columns = "sitc", values = "value")
df4 = df4.div(df4.sum(axis = 1), axis = 0).reset_index()
df4.to_csv("country_vectors_PS4.csv", sep = "\t", index = False)

# 3-digit
df3 = df.copy()
df3["sitc"] = df3["sitc"].astype(str).str.rjust(width = 4, fillchar = '0').map(lambda x : x[:3])
df3 = df3.groupby(by = ["decade", "exporter", "sitc"])["value"].mean().reset_index()
df3["value"] = df3["value"].fillna(0)
df3 = pd.pivot_table(df3, index = ["decade", "exporter"], columns = "sitc", values = "value")
df3 = df3.div(df3.sum(axis = 1), axis = 0).reset_index()
df3.to_csv("country_vectors_PS3.csv", sep = "\t", index = False)

# 2-digit
df2 = df.copy()
df2["sitc"] = df2["sitc"].astype(str).str.rjust(width = 4, fillchar = '0').map(lambda x : x[:2])
df2 = df2.groupby(by = ["decade", "exporter", "sitc"])["value"].mean().reset_index()
df2["value"] = df2["value"].fillna(0)
df2 = pd.pivot_table(df2, index = ["decade", "exporter"], columns = "sitc", values = "value")
df2 = df2.div(df2.sum(axis = 1), axis = 0).reset_index()
df2.to_csv("country_vectors_PS2.csv", sep = "\t", index = False)

# 1-digit
df1 = df.copy()
df1["sitc"] = df1["sitc"].astype(str).str.rjust(width = 4, fillchar = '0').map(lambda x : x[:1])
df1 = df1.groupby(by = ["decade", "exporter", "sitc"])["value"].mean().reset_index()
df1["value"] = df1["value"].fillna(0)
df1 = pd.pivot_table(df1, index = ["decade", "exporter"], columns = "sitc", values = "value")
df1 = df1.div(df1.sum(axis = 1), axis = 0).reset_index()
df1.to_csv("country_vectors_PS1.csv", sep = "\t", index = False)

