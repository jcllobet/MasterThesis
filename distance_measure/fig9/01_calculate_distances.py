import random, sys, copy, subprocess, signal
import numpy as np
import networkx as nx
sys.path.append("..")
from implementation import network_distance as nd

class TimeoutException(Exception):
   pass

def timeout_handler(signum, frame):
   raise TimeoutException

def generate_network(network_type):
   if network_type == 1:
      G = nx.erdos_renyi_graph(100, .095)
      while nx.number_connected_components(G) > 1:
         G = nx.erdos_renyi_graph(100, .095)
      return G
   elif network_type == 2:
      return nx.barabasi_albert_graph(100, 5)
   elif network_type == 3:
      return nx.powerlaw_cluster_graph(100, 5, .01)
   elif network_type == 4:
      bash_command = "./benchmark -N 100 -k 9 -mu 0.1 -t1 3 -t2 1.1 -maxk 50"
      process = subprocess.Popen(bash_command.split(), stdout = subprocess.PIPE)
      output, error = process.communicate()
      return nx.convert_node_labels_to_integers(nx.read_edgelist("network.dat"))

def generate_change(G, infection):
   src = {_: 1 for _ in random.sample(range(100), 10)}
   trg = copy.deepcopy(src)
   for _ in range(5):
      trg = si_step(G, trg, infection)
   return src, trg

signal.signal(signal.SIGALRM, timeout_handler)

with open("dist_table.csv", 'w') as f:
   f.write("lapl\tmmc\tannihil\totp\tspls\tspla\tsplc\tgfte\tmapp\n")
   i = 0
   fail = 0
   while True:
      G = generate_network(random.randint(1, 4))
      if nx.number_connected_components(G) == 1:
         src = {_: random.random() for _ in random.sample(set(G.nodes), random.randint(2, 10))}
         trg = {_: random.random() for _ in random.sample(set(G.nodes), random.randint(2, 10))}
         signal.alarm(2)
         try:
            mapp = float(nd.mapp(src, trg, G))
            signal.alarm(0)
         except:
            fail += 1
            sys.stderr.write("MAPP failed (%1.2f%% failure rate)\n" % (100 * fail / (fail + i)))
            continue
         i += 1
         sys.stderr.write("%s\n" % i)
         line = []
         line.append("%s" % float(nd.ge(src, trg, G)))   
         line.append("%s" % float(nd.mmc(src, trg, G)))
         line.append("%s" % float(nd.annihilation(src, trg, G)))
         line.append("%s" % float(nd.emd(src, trg, G)))   
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "single")))
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "avg")))
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "complete")))   
         line.append("%s" % float(nd.gft(src, trg, G, linkage = "euclidean")))
         line.append("%s" % mapp)
         f.write("%s\n" % '\t'.join(line))
         f.flush()
         if i == 350:
            break
