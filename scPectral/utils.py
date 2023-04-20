import pickle
import hypernetx as hnx
import numpy as np
from sklearn.cluster import KMeans
import csv

"""
Load in a pickled object. 
"""
def load_pickle(save_path: str):
    with open(save_path, 'rb') as instream:
        obj = pickle.load(instream) 
    
    return obj


def visualize_hyperedge_set(model):
    for i in range(len(model.hedges)):
        hnx.draw(hnx.classes.Hypergraph.from_numpy_array(model.hedges[i]), with_node_labels=False, with_edge_labels=False, node_radius=0.5)


"""
Constructs unweighted hypergraph from a vertex x edge incidence matrix H, collapses common edges, and returns the compressed result.
"""
def get_unweighted_graph(H):
    HG = hnx.classes.Hypergraph.from_numpy_array(M = H)
    HG = hnx.classes.Hypergraph.collapse_edges(HG)
    return HG


def get_clusters(k: int, sorted_eigvecs, gene_labels, min_cluster_size):
    skip = 0 
    kmeans = KMeans(n_clusters=k)
    Y = sorted_eigvecs[:, skip:(k + skip)]
    kmeans.fit(Y)
    labels = kmeans.labels_
    clusters = []
    for i in range(k):
        genes = np.argwhere(labels == i).flatten()
        if genes.shape[0] >= min_cluster_size:
            clusters.append(gene_labels[genes])
    
    return clusters

def write_results(data, outfile:str):
    try:
        with open(outfile, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for key, val in data.items():
                writer.writerow(val)
    except IOError:
        print("I/O error")
