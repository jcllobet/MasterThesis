import igraph, sys, os
import numpy as np
import network_sampling as ns

policies = {
   "sphl": ns.APIPolicy(edges_page_size = 10, seconds_per_call = 10),      # 1 e/s
   "spll": ns.APIPolicy(edges_page_size = 5, seconds_per_call = 2),        # 2.5 e/s
   "lphl": ns.APIPolicy(edges_page_size = 100, seconds_per_call = 20),     # 5 e/s
   "lpll": ns.APIPolicy(edges_page_size = 40, seconds_per_call = 4)        # 10 e/s
}

algorithms = {
   "nc": ns.nc,
   "ff": ns.ff,
   "mhrw": ns.mhrw,
   "rw": ns.rw,
   "dfs": ns.dfs,
   "bfs": ns.bfs,
   "snowball": ns.snowball
}

budget = int(sys.argv[1])
osn = sys.argv[2]
alg = sys.argv[3]
net_label = sys.argv[4]
weighted = sys.argv[4][0] == "w"

G_orig = igraph.Graph.Read_Ncol("../data/%s.edgelist" % net_label, directed = False, weights = weighted)

if not os.path.isdir("samples_%s/%s" % (net_label, alg)):
   os.makedirs("samples_%s/%s" % (net_label, alg), exist_ok = True)
if not os.path.isdir("samples_%s/%s/%s" % (net_label, alg, osn)):
   os.makedirs("samples_%s/%s/%s" % (net_label, alg, osn), exist_ok = True)

seeds = np.random.randint(0, G_orig.vcount(), size = 100)
api_engine = ns.APIEngine(G_orig, policies[osn])

for iteration in range(100):
   sys.stderr.write("%s %s %s %s\n" % (alg, budget, osn, iteration))
   G_smpl = algorithms[alg](api_engine, seed = seeds[iteration], budget = budget, weighted = weighted)
   G_smpl.write_edgelist("samples_%s/%s/%s/%s_%s" % (net_label, alg, osn, budget, iteration))

