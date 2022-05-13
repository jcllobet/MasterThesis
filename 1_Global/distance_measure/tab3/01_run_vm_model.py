import random, sys, copy, subprocess, signal
import numpy as np
import networkx as nx
sys.path.append("..")
from implementation import network_distance as nd

class TimeoutException(Exception):
   pass

def timeout_handler(signum, frame):
   raise TimeoutException

def si_step(G, seeds, contagion_parameter):
   new_seeds = set()
   infected = set(seeds.keys())
   for n in G.nodes:
      neighbors = set(G.neighbors(n))
      if (len(neighbors & infected) / len(neighbors)) > contagion_parameter:
         new_seeds.add(n)
   return {_: 1 for _ in new_seeds}

def generate_network(network_type):
   if network_type == "er":
      G = nx.erdos_renyi_graph(100, .095)
      while nx.number_connected_components(G) > 1:
         G = nx.erdos_renyi_graph(100, .095)
      return G
   elif network_type == "ba":
      return nx.barabasi_albert_graph(100, 5)
   elif network_type == "pc":
      return nx.powerlaw_cluster_graph(100, 5, .01)
   elif network_type == "lfr":
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

network_type = sys.argv[1]
signal.signal(signal.SIGALRM, timeout_handler)

with open("vm_table_%s.csv" % network_type, 'w') as f:
   f.write("infection\tge\tmmc\tannihil\temd\tspl single\tspl avg\tspl complete\tgft\tmapp\tcount\n")
   i = 0
   fail = 0
   while True:
      G = generate_network(network_type)
      infection = random.random()
      line = []
      line.append("%s" % infection)
      src, trg = generate_change(G, infection)
      if len(trg) > 0:
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
         line.append("%s" % float(nd.ge(src, trg, G)))   
         line.append("%s" % float(nd.mmc(src, trg, G)))
         line.append("%s" % float(nd.annihilation(src, trg, G)))
         line.append("%s" % float(nd.emd(src, trg, G)))   
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "single")))
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "avg")))
         line.append("%s" % float(nd.spl(src, trg, G, linkage = "complete")))   
         line.append("%s" % float(nd.gft(src, trg, G)))
         line.append("%s" % float(mapp))
         line.append("%s" % abs(len(src) - len(trg)))
         f.write("%s\n" % '\t'.join(line))
         f.flush()
         if i == 350:
            break
