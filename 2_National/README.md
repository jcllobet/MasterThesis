# Vector Distances in the Propagation of COVID-19 on a Network over Time
This project investigates how the metric of Generalized Euclidean Distance, GED, between node vectors can add to our understanding of how a
spreading event unfolds over time in a network. This metric measures the distance covered by a spreading phenomenon in the nodes of a network between
two steps in time, in a way that takes the topology of the network into account.

The investigation is done by first constructing a mobility network with nodes representing municipalities and weighted edges representing the relative flow of traffic between municipalities. The edges are constructed using passenger count data from the train service provider DSB, and by generating links to isolated nodes using the Gravity Model. The network is then backboned using a threshold, based on the High Salience Skeleton scores of the edges. This results in a network G = (V, E) with |V| = 98, |E| = 1303, and an average degree k = 26.52.
The topology of this network is then used for calculating the GED between node vectors with a weekly frequency for the period of May 2020 - March 2021. These vectors are based on COVID-19 surveillance data provided by SSI and contain the amount of positive tests for each municipality normalized by individual population size.

It is found that the GED metric is informative. This is because the information it captures is different from, i.e., the total activation state in the network or even the standard deviation of the activation states. The GED metric is also found to be meaningful; the pattern we see, when the GED doesn’t follow the standard deviation, intuitively makes sense when inspecting the positions of the nodes with outlying activation state values in the network. 
This metric could be used in future work for evaluating the effectiveness of different network-based non-pharmaceutical interventions, i.e., varying levels of strictness in nationwide stay-at-home orders, in slowing the spread of disease.

# Notebooks

The notebooks GED.ipynb and backboning.ipynb utilize libraries for backboning and calculating network distance. These libraries were generously provided by Michele Coscia on his website. The libraries in question are [Network Backboning](https://www.michelecoscia.com/?page_id=287) and [Node Vector Distance](https://www.michelecoscia.com/?page_id=1733). The associated scripts need to be placed in a folder called modules in the same directory as the notebooks in order for the notebooks to run.

1. **Aggregating the Mobility Data**  
notebook: AggregatingOD.ipynb  
description: Loads in the OD matrix from an excel file, maps all the stations  
to their respective municipalities and aggregates the data by municipality.  
Connections to Germany and Sweden are removed from the data.  

2. **Scraping Inter-Municipal Distances Online**  
notebook: Scraper.ipynb  
description: The web scraper used to get all the inter-municipal distances  
from https://www.afstande.com/. Outputs a file containing per per row a  
string with the travel time in minutes and kilometers. The file contains one  
line per existing edge in the network plus one for each possible combination  
of the isolated nodes with the already connected nodes.

3. **Estimating Edge Weights**  
notebook: GravityModel.ipynb  
description: Loads in the inter-municipal distance data, which is then used  
to derive the best value of the constant k to use in the gravity model. Based  
on the chosen value of k edge weights are estimated for all the new edges  
which are then added to the network connecting the previously isolated  
nodes to every other node in the network. Self-loops are also removed from  
the graph here.  

4. **Backboning the Network**  
notebook: Backboning.ipynb  
description: Prepares the correct input format for the Network Backboning  
module. Then calculates the HSS scores for all edges in the network, and  
performs the thresholding with the best found parameter. Also extracts basic  
statistics for the final network topology and plots degree and edge weight  
distributions.  

5. **Calculating the GEDs**  
notebook: GED.ipynb  
description: The COVID-19 data is aggregated by week, and then turned  
into node 44 vectors using three different methods of normalisation. These  
three sets of 44 vectors are then used with the network topology to calculate  
the GEDs. The results from the three different vector sets are plotted.  

6. **Analysing the Results**  
notebook: Analysis.ipynb  
description: The diff vectors are calculated for each of the three identified  
peaks, and base statistics for each of the six time intervals are calculated  
along with the increases in GED and standard deviations. The distributions  
of diff are plotted, and the node activation states at the middle peak are  
investigated.  
