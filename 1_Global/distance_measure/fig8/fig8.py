import sys
import networkx as nx
sys.path.append("..")
from implementation import network_distance as nd

with open("fig8.csv", 'w') as f:
   f.write("chain length\tge\tmmc\tannihil\temd\tspl single\tspl avg\tspl complete\tgft\tmapp\n")
   for n in range(2, 101):
      print(n)
      G = nx.path_graph(n)
      src = {0: 1}
      trg = {n - 1: 1}
      ge = nd.ge(src, trg, G)
      mmc = nd.mmc(src, trg, G)
      annihil = nd.annihilation(src, trg, G)
      emd = nd.emd(src, trg, G)
      spl_s = nd.spl(src, trg, G, linkage = "single")
      spl_a = nd.spl(src, trg, G, linkage = "avg")
      spl_c = nd.spl(src, trg, G, linkage = "complete")
      gft_e = nd.gft(src, trg, G, linkage = "euclidean")
      mapp = nd.mapp(src, trg, G)
      f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (n, ge, mmc, annihil, emd, spl_s, spl_a, spl_c, gft_e, mapp))
      f.flush()

