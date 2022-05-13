import pandas as pd
import numpy as np
import networkx as nx
import network_distance as nd

G = nx.read_edgelist("/path/to/anobii/social/net") # Space-separated user-user network, 4core

nodes = list(G.nodes)
ge_Q = nd._ge_Q(G, hermitian = True)
gft_v = nd._gft_v(G)
spl = nd.calculate_spl(G, nodes, 21, return_as_dict = True)
diameter = nx.diameter(G)
mmc_P = nd._mmc_P(G, diameter)

ann_Q = nd._annihilation_Q(nx.adjacency_matrix(G, nodelist = nodes).todense().astype(float))

df = pd.read_csv("/path/to/user/data", sep = "\t") # Tab-separated user-book-week table
measure_book_week_distance = []
book_n = len(set(df["book"]))
i = 0
for book in set(df["book"]):
   print("%1.2f%%" % (100 * i / book_n))
   df_b = df[df["book"] == book]
   w = []
   for _ in range(1, 7):
      w.append({node: 1 for node in G.nodes if node in set(df_b[df_b["week"] == _]["user"])})
   for _ in range(len(w) - 1):
      measure_book_week_distance.append(("lapl", book, _, nd.ge(w[_], w[_ + 1], G, Q = ge_Q)))
      measure_book_week_distance.append(("mmc", book, _, nd.mmc(w[_], w[_ + 1], G, P = mmc_P, time_steps = diameter)))
      measure_book_week_distance.append(("annihil", book, _, nd.annihilation(w[_], w[_ + 1], G, Q = ann_Q)))
      measure_book_week_distance.append(("otp", book, _, nd.emd(w[_], w[_ + 1], G, shortest_path_lengths = spl)))
      measure_book_week_distance.append(("spls", book, _, nd.spl(w[_], w[_ + 1], G, linkage = "single", shortest_path_lengths = spl)))
      measure_book_week_distance.append(("spla", book, _, nd.spl(w[_], w[_ + 1], G, linkage = "avg", shortest_path_lengths = spl)))
      measure_book_week_distance.append(("splc", book, _, nd.spl(w[_], w[_ + 1], G, linkage = "complete", shortest_path_lengths = spl)))
      measure_book_week_distance.append(("gft", book, _, nd.gft(w[_], w[_ + 1], G, v = gft_v)))
   i += 1

df = pd.DataFrame(data = measure_book_week_distance, columns = ("measure", "book", "week", "distance"))
df.to_csv("measure_book_week_distance", sep = "\t", index = False)
