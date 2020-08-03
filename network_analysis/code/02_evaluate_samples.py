import sys
import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn.metrics import normalized_mutual_info_score

def generate_rwrw_dd(dd_smpl, nnodes):
   dd_smpl = pd.DataFrame(dd_smpl, columns = ("d",))
   dd_smpl["d_inv"] = 1.0 / dd_smpl["d"]
   denominator = dd_smpl["d_inv"].sum()
   dd_smpl = dd_smpl.groupby(by = "d").sum().reset_index()
   dd_smpl["p"] = dd_smpl["d_inv"] / denominator
   dd_smpl["nnodes"] = (dd_smpl["p"] * nnodes).map(np.ceil).astype(int)
   result = []
   for index, row in dd_smpl.iterrows():
      result.extend([row["d"]] * row["nnodes"].astype(int))
   return np.array(result)

def dd_test(G_smpl, dd_orig):
   dd_smpl = np.array(list(dict(G_smpl.degree()).values()))
   dd_smpl = dd_smpl[dd_smpl > 0]
   if mthd == "rwrw":
      dd_smpl = generate_rwrw_dd(dd_smpl, len(dd_orig))
   return stats.ks_2samp(dd_orig, dd_smpl)[0]

def centr_test(G_smpl, centr_orig):
   centr_smpl = pd.DataFrame().from_dict(dict(G_smpl.degree()), orient = "index").reset_index().rename(columns = {"index": "node", 0: "centr_smpl"})
   centr_smpl = centr_smpl[centr_smpl["centr_smpl"] > 0].merge(centr_orig, on = "node")
   return centr_smpl.corr(method = "spearman").loc["centr_orig", "centr_smpl"]

def print_results(f, budget, osn, degdi_results, centr_results, assort, disassort, comm_quality, comm_nmi):
   line = []
   line.append(str(budget))
   line.append(str(osn))
   line.append(str(np.mean(degdi_results)))
   line.append(str(stats.sem(degdi_results)))
   line.append(str(np.mean(centr_results)))
   line.append(str(stats.sem(centr_results)))
   line.append(str(np.mean(assort)))
   line.append(str(stats.sem(assort)))
   line.append(str(np.mean(disassort)))
   line.append(str(stats.sem(disassort)))
   line.append(str(np.mean(comm_quality)))
   line.append(str(stats.sem(comm_quality)))
   line.append(str(np.mean(comm_nmi)))
   line.append(str(stats.sem(comm_nmi)))
   f.write("%s\n" % '\t'.join(line))
   f.flush()

def run_test(G, f, osn, mthd, budget):
   dd_orig = np.array(list(dict(G.degree()).values()))
   centr_orig = pd.DataFrame().from_dict(dict(G.degree()), orient = "index").reset_index().rename(columns = {"index": "node", 0: "centr_orig"})
   assort_orig = nx.attribute_assortativity_coefficient(G, "assortative")
   disassort_orig = nx.attribute_assortativity_coefficient(G, "disassortative")
   comms_orig = list(nx.algorithms.community.label_propagation.label_propagation_communities(G))
   commquality_orig = nx.algorithms.community.quality.coverage(G, comms_orig)
   comms_orig = {n: i for i in range(len(comms_orig)) for n in comms_orig[i]}
   try:
      degdi_results = []
      centr_results = []
      assort = []
      disassort = []
      comm_quality = []
      comm_nmi = []
      for run in range(100):
         G_smpl = nx.read_edgelist("samples_%s/%s/%s/%s_%s" % (net_label, mthd, osn, budget, run), nodetype = int)
         nx.set_node_attributes(G_smpl, assortative_attr, "assortative")
         nx.set_node_attributes(G_smpl, disassortative_attr, "disassortative")
         comms_smpl = list(nx.algorithms.community.label_propagation.label_propagation_communities(G_smpl))
         commquality_smpl = nx.algorithms.community.quality.coverage(G_smpl, comms_smpl)
         comms_smpl = {n: i for i in range(len(comms_smpl)) for n in comms_smpl[i]}
         centr_results.append(centr_test(G_smpl, centr_orig))
         assort.append(abs(nx.attribute_assortativity_coefficient(G_smpl, "assortative") - assort_orig))
         disassort.append(abs(nx.attribute_assortativity_coefficient(G_smpl, "disassortative") - disassort_orig))
         degdi_results.append(dd_test(G_smpl, dd_orig))
         comm_quality.append(abs(commquality_smpl - commquality_orig))
         comm_nmi.append(normalized_mutual_info_score([comms_orig[n] for n in G_smpl.nodes], [comms_smpl[n] for n in G_smpl.nodes]))
      print_results(f, budget, osn, degdi_results, centr_results, assort, disassort, comm_quality, comm_nmi)
   except IOError:
      print("!!!!!!!!!!!!!!IO ERROR!!!!!!!!!!!!!! %s %s %s %s" % (osn, mthd, budget, run))

mthd = sys.argv[1]
net_label = sys.argv[2]

try:
   G = nx.read_edgelist("../data/%s.edgelist" % net_label, nodetype = int)
except TypeError:
   G = nx.read_edgelist("../data/%s.edgelist" % net_label, nodetype = int, data = [("weight", int),])
   
node_attrs = pd.read_csv("../data/%s.nodelist" % net_label, sep = "\t")
assortative_attr = {row["node"]: row["assortative"] for index, row in node_attrs.iterrows()}
disassortative_attr = {row["node"]: row["disassortative"] for index, row in node_attrs.iterrows()}
nx.set_node_attributes(G, assortative_attr, "assortative")
nx.set_node_attributes(G, disassortative_attr, "disassortative")

with open("%s/%s.csv" % (net_label, mthd), 'w') as f:
   f.write("budget\tosn\tdd_ks_mean\tdd_ks_sem\tcentr_spear_mean\tcentr_spear_sem\tassort_mae\tassort_mae_sem")
   f.write("\tdisassort_mae\tdisassort_mae_sem\tcommstrength_mae\tcommstrength_mae_sem\tcomm_nmi\tcomm_nmi_sem\n")
   for budget in (60, 120, 300, 600, 900, 1800):
      sys.stderr.write("Loading %s-%s\n" % (mthd, net_label))
      for osn in ("sphl", "lphl", "lpll", "spll"):
         sys.stderr.write("Testing %s-%s-%s-%s\n" % (mthd, net_label, osn, budget))
         run_test(G, f, osn, mthd, budget)
