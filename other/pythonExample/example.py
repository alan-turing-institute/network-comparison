import numpy as np
import networkx as nx
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
orca=importr("orca")
igraph = importr('igraph')
netdist =importr('netdist')
graphConstructor=rpy2.robjects.r['graph.adjacency']

def makeRgdd(G):
    B=nx.to_numpy_matrix(G)
    B=np.array(B)
    nr,nc = B.shape
    Br = robjects.r.matrix(B, nrow=nr, ncol=nc)
    Gr=graphConstructor(Br,mode ="undirected")
    Ggdd=netdist.gdd(Gr)
    return Ggdd

G1=nx.Graph()
G1.add_edges_from([[0,1],[1,2],[2,3]])
G1gdd=makeRgdd(G1)


G2=nx.Graph()
G2.add_edges_from([[0,1],[1,2],[2,3],[3,4]])
G2gdd=makeRgdd(G2)

print netdist.net_emd(G1gdd,G1gdd)
print netdist.net_emd(G2gdd,G2gdd)
print netdist.net_emd(G1gdd,G2gdd)
