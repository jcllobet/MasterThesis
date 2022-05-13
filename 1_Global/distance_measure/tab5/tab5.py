import sys, copy
import pandas as pd
import networkx as nx
from collections import defaultdict
sys.path.append("..")
from implementation import network_distance as nd

outbreaks = pd.read_csv("outbreaks.csv", sep = "\t", dtype = {"month": str})
flights = nx.read_edgelist("flights_undirected.csv")

nodes = list(flights.nodes)
ge_Q = nd._ge_Q(flights)
gft_v = nd._gft_v(flights)
spl = nd.calculate_spl(flights, nodes, 10, return_as_dict = True)
diameter = nx.diameter(flights)
mmc_P = nd._mmc_P(flights, diameter)
ann_Q = nd._annihilation_Q(nx.adjacency_matrix(flights, nodelist = nodes).todense().astype(float))

nodemap = {list(flights.nodes)[i]: i for i in range(len(flights.nodes))}
flights_mapp = nx.relabel_nodes(flights, nodemap)
spl_mapp = nd.calculate_spl(flights_mapp, list(flights_mapp.nodes), 10, return_as_dict = True)

disease_measure_distances = defaultdict(lambda : defaultdict(list))
for disease in set(outbreaks["disease"]):
   df = outbreaks[outbreaks["disease"] == disease]
   months = sorted(list(df["month"].drop_duplicates()))
   countries_touched = set()
   for i in range(1, len(months)):
      sys.stderr.write("%s %s of %s\n" % (disease, i , len(months)))
      yeardelta = int(months[i][:4]) - int(months[i - 1][:4])
      monthdelta = int(months[i][4:]) - int(months[i - 1][4:])
      monthlag = (yeardelta * 12) + monthdelta
      src = {country: 1 for country in countries_touched if country in flights.nodes}
      src.update({country: 1 for country in df[df["month"] == months[i - 1]]["country"] if country in flights.nodes})
      trg = copy.deepcopy(src)
      trg.update({country: 1 for country in df[df["month"] == months[i]]["country"] if country in flights.nodes})
      countries_touched |= set(trg.keys())
      src_mapp = {nodemap[x]: src[x] for x in src}
      trg_mapp = {nodemap[x]: trg[x] for x in trg}
      disease_measure_distances[disease]["lapl"].append(nd.ge(src, trg, flights, Q = ge_Q) / monthlag)
      disease_measure_distances[disease]["mmc"].append(nd.mmc(src, trg, flights, P = mmc_P, time_steps = diameter) / monthlag)
      disease_measure_distances[disease]["annihil"].append(nd.annihilation(src, trg, flights, Q = ann_Q) / monthlag)
      disease_measure_distances[disease]["otp"].append(nd.emd(src, trg, flights, shortest_path_lengths = spl) / monthlag)
      disease_measure_distances[disease]["spls"].append(nd.spl(src, trg, flights, linkage = "single", shortest_path_lengths = spl) / monthlag)
      disease_measure_distances[disease]["spla"].append(nd.spl(src, trg, flights, linkage = "avg", shortest_path_lengths = spl) / monthlag)
      disease_measure_distances[disease]["splc"].append(nd.spl(src, trg, flights, linkage = "complete", shortest_path_lengths = spl) / monthlag)
      disease_measure_distances[disease]["gft"].append(nd.gft(src, trg, flights, v = gft_v) / monthlag)
      disease_measure_distances[disease]["mapp"].append(nd.mapp(src_mapp, trg_mapp, flights_mapp, shortest_path_lengths = spl_mapp) / monthlag)

print("\\begin{table}")
print("\\centering")
print("\\begin{tabular}{l|llll}")
print("Measure & Ebola & Dengue & Avian Flu & Zika\\\\")
print("\\hline")
for distance in ("lapl", "mmc", "annihil", "otp", "spls", "spla", "splc", "gft", "mapp"):
   line = []
   line.append(distance)
   for disease in ("ebola", "dengue", "avianflu", "zika"):
      line.append("%1.4f" % (sum(disease_measure_distances[disease][distance]) / len(disease_measure_distances[disease][distance])))
   print("%s\\\\" % " & ".join(line))

print("\\end{tabular}")
print("\\caption{}")
print("\\label{tab:si}")
print("\\end{table}")
