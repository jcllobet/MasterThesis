import math, random, copy, igraph
import pandas as pd
import backboning as bb
from itertools import chain
from collections import deque

class APIPolicy():
   def __init__(self, edges_page_size = 30, seconds_per_call = 1, latency = 0):
      self.edges_page_size = edges_page_size # What's the maximum number of edges we return per query?
      self.seconds_per_call = seconds_per_call # How much does a call costs in seconds
      self.latency = latency

class APIEngine():
   def __init__(self, network, policy):
      self.network = network
      self.policy = policy

   def request_neighbors(self, node, budget_left, pages = None, weight = False):
      if pages is None:
         num_calls_available = budget_left // (self.policy.latency + self.policy.seconds_per_call) # This tells me how many calls can the crawler still make
      else:
         num_calls_available = min(pages, budget_left // (self.policy.latency + self.policy.seconds_per_call)) # Some crawlers only want a certain number of pages per node
      neighbors = list(self.network.neighbors(node))[:int(num_calls_available) * self.policy.edges_page_size] # Get the node's neighbors, but only the ones the crawler can afford
      num_calls = int(math.ceil(len(neighbors) / float(self.policy.edges_page_size))) # The number of calls needed is the number of neighbors divided by the edge page size
      cost = (self.policy.latency + self.policy.seconds_per_call) * num_calls # The query cost is how much each call costs in seconds time the number of calls required
      if weight:
         return neighbors, cost, self.network.es[neighbors]["weight"]
      else:
         return neighbors, cost

def bfs(api_engine, seed = 0, budget = 100, directed = False):
   node_fifo = deque([seed,])
   explored_nodes = set()
   edges = {}
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      try:
         node = node_fifo.pop()
         while node in explored_nodes:
            node = node_fifo.pop()
         explored_nodes.add(node)
         neighbors, cost = api_engine.request_neighbors(node, budget)
         edges[node] = neighbors
         node_fifo.extendleft(neighbors) # Fills newcomers to the left, but only the ones that have not been already explored
         budget -= cost
      except IndexError:
         break
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   return G

def dfs(api_engine, seed = 0, budget = 100, directed = False):
   node_lifo = deque([seed,])
   explored_nodes = set()
   edges = {}
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      try:
         node = node_lifo.pop()
         while node in explored_nodes:
            node = node_lifo.pop() # Pops from the right
         explored_nodes.add(node)
         neighbors, cost = api_engine.request_neighbors(node, budget)
         edges[node] = neighbors
         node_lifo.extend(neighbors) # Fills newcomers to the right, but only the ones that have not been already explored
         budget -= cost
      except IndexError:
         break
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   return G

def snowball(api_engine, seed = 0, budget = 100, k = 1, directed = False):
   node_fifo = deque([seed,])
   explored_nodes = set()
   edges = {}
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      try:
         node = node_fifo.pop()
         while node in explored_nodes:
            node = node_fifo.pop() # Pops from the right
         explored_nodes.add(node)
         neighbors, cost = api_engine.request_neighbors(node, budget, pages = k)
         edges[node] = neighbors
         node_fifo.extendleft(neighbors) # Fills newcomers to the left, but only the ones that have not been already explored
         budget -= cost
      except IndexError:
         break
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   return G

def rw(api_engine, seed = 0, budget = 100, restart_p = .1, return_explored_node_set = False, directed = False):
   next_nodes = [seed,]
   explored_nodes = set()
   edges = {}
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      attempts = 0
      if random.random() < restart_p: # We restard from the seed with a 15% random chance
         node = random.randint(0, api_engine.network.vcount() - 1)
      else:
         try:
            node = random.choice(next_nodes) # Selects a neighbor at random
         except IndexError:
            break
      while node in explored_nodes: # If we never saw this node before, we need to get its info...
         attempts += 1
         if attempts > 1000:
            break
         if random.random() < restart_p:
            node = random.randint(0, api_engine.network.vcount() - 1)
         else:
            next_nodes = edges[node]
            try:
               node = random.choice(next_nodes) # Selects a neighbor at random
            except IndexError:
               break
      explored_nodes.add(node)
      neighbors, cost = api_engine.request_neighbors(node, budget)
      edges[node] = neighbors
      budget -= cost
      next_nodes = edges[node] # The new nodes to explore are the neighbors, even if we already explored them (true random walk)
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   if return_explored_node_set:
      return G, explored_nodes, budget
   else:
      return G

def mhrw(api_engine, seed = 0, budget = 100, restart_p = .1, directed = False):
   next_nodes = [seed,]
   explored_nodes = set()
   degree_x = float("inf")
   edges = {}
   attempts = 0
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      if random.random() < restart_p: # We restard from the seed with a 15% random chance
         node = random.randint(0, api_engine.network.vcount() - 1)
      else:
         try:
            node = random.choice(next_nodes) # Selects a neighbor at random
         except IndexError:
            break
      while node in explored_nodes: # If we never saw this node before, we need to get its info...
         attempts += 1
         if attempts > 1000:
            break
         if random.random() < restart_p:
            node = random.randint(0, api_engine.network.vcount() - 1)
         else:
            next_nodes = edges[node]
            try:
               node = random.choice(next_nodes) # Selects a neighbor at random
            except IndexError:
               break
      neighbors, cost = api_engine.request_neighbors(node, budget)
      budget -= cost
      degree_y = float(len(neighbors))
      try:
         if random.random() < (degree_x / degree_y): # We follow the link only if this test succeeeds. This corrects the degree bias of a vanilla RW.
            explored_nodes.add(node)
            edges[node] = neighbors
            next_nodes = neighbors # The new nodes to explore are the neighbors, even if we already explored them (true random walk)
            degree_x = degree_y
      except ZeroDivisionError:
         continue
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   return G

def ff(api_engine, seed = 0, budget = 100, burn_p = .2, directed = False):
   node_fifo = deque([seed,])
   explored_nodes = set()
   edges = {}
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      try:
         node = node_fifo.pop()
         while node in explored_nodes:
            node = node_fifo.pop()
         if random.random() < burn_p: # True if we want to burn
            explored_nodes.add(node)
            neighbors, cost = api_engine.request_neighbors(node, budget)
            edges[node] = neighbors
            node_fifo.extendleft(neighbors) # Fills newcomers to the left, but only the ones that have not been already explored
            budget -= cost
         else: # If we didn't burn the node, we put it back for future attempts
            node_fifo.appendleft(node)
      except IndexError:
         break
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(s, t) for s in edges for t in edges[s]])
   else:
      G.add_edges(set([(min(s, t), max(s, t)) for s in edges for t in edges[s]]))
   return G

def nc(api_engine, seed = 0, budget = 100, directed = False, weighted = False):
   edges = pd.DataFrame(columns = ("src", "trg", "nij"))
   explored_nodes = set()
   node = seed
   while budget > api_engine.policy.seconds_per_call and len(explored_nodes) < api_engine.network.vcount():
      if weighted:
         neighbors, cost, weights = api_engine.request_neighbors(node, budget, weight = weighted)
         new_edges = pd.DataFrame(data = chain.from_iterable(((node, neighbors[i], weights[i]), (neighbors[i], node, weights[i])) for i in range(len(neighbors))), columns = ("src", "trg", "nij"))
      else:
         neighbors, cost = api_engine.request_neighbors(node, budget, weight = weighted)
         new_edges = pd.DataFrame(data = chain.from_iterable(((node, neighbors[i], 1), (neighbors[i], node, 1)) for i in range(len(neighbors))), columns = ("src", "trg", "nij"))
      budget -= cost
      explored_nodes.add(node)
      edges = pd.concat([bb.make_symmetric(edges[["src", "trg", "nij"]]), new_edges], sort = False).drop_duplicates()
      edges["nij"] = edges["nij"].astype(int)
      edges = bb.noise_corrected(edges, undirected = True)
      edges["z"] = edges["score"] / edges["sdev_cij"]
      edge_follow = edges[~(edges["src"].isin(explored_nodes) & edges["trg"].isin(explored_nodes))]
      edge_follow = edge_follow[edge_follow["z"] == edge_follow["z"].max()].sample(1)
      node = edge_follow["src"].values[0]
      if node in explored_nodes:
         node = edge_follow["trg"].values[0]
   G = igraph.Graph(n = api_engine.network.vcount(), directed = directed)
   if directed:
      G.add_edges([(row["src"], row["trg"]) for index, row in edges.iterrows()])
   else:
      G.add_edges(set([(row[["src", "trg"]].min(), row[["src", "trg"]].max()) for index, row in edges.iterrows()]))
   return G

