import sys, copy
import numpy as np
import pandas as pd
import networkx as nx
sys.path.append("..")
from implementation import network_distance as nd

def eci(df):
   products = df.drop(["exporter", "decade"], 1).columns
   product_share = df.set_index("exporter").groupby(by = "decade").mean().reset_index()
   result = df.merge(product_share, on = "decade", suffixes = ("_own", "_world"))
   rca = result[["exporter", "decade"]]
   for product in products:
      rca[product] = result["%s_own" % product] / result["%s_world" % product]
   rca = (rca.fillna(0).set_index(["exporter", "decade"]) > 1).reset_index()
   result = rca[["exporter",]].drop_duplicates()
   for decade in ("1970s", "1980s", "1990s", "2000s"):
      eci = rca[rca["decade"] == decade].drop("decade", 1).set_index("exporter")
      ubiquity = eci.sum(axis = 0)
      diversity = eci.sum(axis = 1)
      ubiquity = ubiquity[ubiquity > 0]
      eci = eci[ubiquity.index]
      Q = eci / ubiquity
      R = (eci.T / diversity).T
      eci = Q.dot(R.T)
      values, vectors = np.linalg.eig(eci.T)
      values = np.real(values)
      vectors = np.real(vectors)
      sorted_index = values.argsort()[::-1]
      values = values[sorted_index]
      result[decade] = vectors[:,1]
   result = result.set_index("exporter").unstack().reset_index()
   result.columns = ("decade", "country", "eci")
   result["decade"] = result["decade"].map({"1970s": 1, "1980s": 2, "1990s": 3, "2000s": 4})
   return result.set_index(["country", "decade"]).to_dict()["eci"]

def diversity(df):
   products = df.drop(["exporter", "decade"], 1).columns
   product_share = df.set_index("exporter").groupby(by = "decade").mean().reset_index()
   df2 = df.merge(product_share, on = "decade", suffixes = ("_own", "_world"))
   rca = df2[["exporter", "decade"]]
   for product in products:
      rca[product] = df2["%s_own" % product] / df2["%s_world" % product]
   rca = (rca.set_index(["exporter", "decade"]) > 1).sum(axis = 1).reset_index()
   rca = pd.pivot_table(data = rca, index = "exporter", columns = "decade", values = 0).reset_index()
   rca[1] = rca["1970s"] - rca["1960s"]
   rca[2] = rca["1980s"] - rca["1970s"]
   rca[3] = rca["1990s"] - rca["1980s"]
   rca[4] = rca["2000s"] - rca["1990s"]
   rca = rca.set_index("exporter").drop(["1960s", "1970s", "1980s", "1990s", "2000s"], 1).unstack().reset_index()
   rca.columns = ("decade", "country", "diversity")
   return rca.set_index(["country", "decade"]).to_dict()["diversity"]

for digit in (1, 2, 3, 4):
   df = pd.read_csv("country_vectors_PS%s.csv" % digit, sep = "\t")
   if digit == 4:
      G = nx.read_edgelist("PS_SITC_edges_%sdigit" % digit, nodetype = int)
   else:
      G = nx.read_edgelist("PS_SITC_edges_%sdigit" % digit)
   nodes = sorted(G.nodes)
   _ = nx.OrderedGraph()
   _.add_nodes_from(nodes)
   _.add_edges_from(G.edges)
   G = _
   ge_Q = nd._ge_Q(G)
   gft_v = nd._gft_v(G)
   spl = nd.calculate_spl(G, nodes, 10, return_as_dict = True)
   diameter = nx.diameter(G)
   mmc_P = nd._mmc_P(G, diameter)
   ann_Q = nd._annihilation_Q(nx.adjacency_matrix(G, nodelist = nodes).todense().astype(float))
   ecidf = eci(df)
   diversitydf = diversity(df)
   with open("country_distances_%s.csv" % digit, 'w') as f:
      f.write("country\tdecade\tmeasure\tdistance\n")
      for country in sorted(list(set(df["exporter"]))):
         sys.stderr.write("%s %s\r" % (country, digit))
         _ = df[df["exporter"] == country].set_index(["exporter"])
         decades = [_[_["decade"] == d].set_index("decade").values[0] for d in sorted(list(set(df["decade"])))]
         trg = {nodes[j]: decades[0][j] for j in range(len(nodes))}
         for i in range(1, len(decades)):
            src = copy.deepcopy(trg) 
            trg = {nodes[j]: decades[i][j] for j in range(len(nodes))}
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "lapl", nd.ge(src, trg, G, Q = ge_Q)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "mmc", nd.mmc(src, trg, G, P = mmc_P, time_steps = diameter)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "annihil", nd.annihilation(src, trg, G, Q = ann_Q)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "otp", nd.emd(src, trg, G, shortest_path_lengths = spl)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "spls", nd.spl(src, trg, G, linkage = "single", shortest_path_lengths = spl)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "spla", nd.spl(src, trg, G, linkage = "avg", shortest_path_lengths = spl)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "splc", nd.spl(src, trg, G, linkage = "complete", shortest_path_lengths = spl)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "gft", nd.gft(src, trg, G, v = gft_v)))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "eci", ecidf[country, i]))
            f.write("%s\t%s\t%s\t%s\n" % (country, i, "count", diversitydf[country, i]))
            f.flush()

sys.stderr.write("Done!\n")
