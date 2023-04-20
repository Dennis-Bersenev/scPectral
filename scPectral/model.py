import numpy as np
from copy import deepcopy
import pandas as pd
import glob
import utils
import hypernetx as hnx

class HyperGraph:
    def __init__(self, n_cells, n_genes, infile_dir, gene_labels):
        self.n_cells = n_cells
        self.n_genes = n_genes
        self.HG = hnx.Hypergraph()
        self.edge_hash = 0

        PWgraph = np.zeros((n_cells, n_genes, n_genes))
        self.gene_labels = gene_labels
        
        for i in range(n_cells):
            file = glob.glob(infile_dir + "*gene{i}.csv".format(i=(i+1)))[0]
            df_temp = pd.read_csv(file)
            PWgraph[i, :, :] = np.copy(df_temp.values)

        self.PWgraph = PWgraph    

    """
    Merges tails with common heads to identify all the multi-arity relations to represent as hedges and updates the supplied incidence matrix with these final hedges.
    _edges: a list of tuples (tail_set, head_set) of all the initial single-tailed edges (note: they are sets, not lists!)
    B: a list of n_gene-sized integer numpy arrays to represent the incidence matrix to update.
    """
    def __merge_edges(self, _edges, B, W, cell_nw):
        
        # init
        edges = deepcopy(_edges)
        m = len(edges)

        # Stop when there's only a single edge left (not possible to merge anything else)
        while (m >= 2):
            
            new_edges = [] 
            for i in range(0, m, 1):
                e1 = edges[i]
                for j in range(i + 1, m, 1):
                    e2 = edges[j]
                    common_heads = e1[1].intersection(e2[1])
                    if (len(common_heads) > 0):
                        
                        # The new edge
                        new_tails = e1[0].union(e2[0]) 
                        new_edge = ([new_tails, common_heads])
                        new_edges.append(new_edge)

                        # Updating old edges (their head sets)
                        e1[1] = set(e1[1] - common_heads)
                        e2[1] = set(e2[1] - common_heads)

            
            # Add all edges with non-null head sets to the incidence matrix
            for edge in edges:
                if len(edge[1]) > 0:
                    col = np.zeros(self.n_genes)
                    w = 0
                    for tail in edge[0]:
                        col[tail] = 1 
                        for head in edge[1]:
                            w += cell_nw[tail, head]
                            col[head] = 1
                    
                    B.append(np.copy(col))
                    # add edge weight to edge matrix
                    W.append(w)
                    
                    elt_indices = np.argwhere(col > 0).flatten()
                    elts = self.gene_labels[elt_indices].tolist()
                    hnx_edge = hnx.Entity('e{i}'.format(i = self.edge_hash), elements=elts, weight=w)
                    self.edge_hash += 1
                    self.HG.add_edge(hnx_edge)
            
            edges = deepcopy(new_edges)
            m = len(edges)


    """ 
    Build an incidence matrix to represent the hyperedges found within the given cell's pairwise network. 
    cell_index: the index into the PWscores array for the desired cell.
    p: percentage parameter to use in tolerance calculation; decides which relationships relative to the maximum PW TE. 
    returns: B the incidence matrix of the constructed hypergraph. 
    """
    def __extract_hyper_edges(self, cell_index, p):
        # init
        unique_genes = set()
        cell_nw = self.PWgraph[cell_index, :, :]

        B = []

        # The weight matrix
        W = []

        # all scores within whatever % of the max
        tolerance = np.max(cell_nw) * p

        tails = set(np.arange(start=0, stop=self.n_genes, step=1))

        while (len(tails) > 0):
            edges = []
            new_tails = set()
            heads_to_remember = set()

            for tail in tails:
                heads = set(np.flatnonzero(cell_nw[tail, :] > tolerance))
                
                # Termination condition: once there's no new heads added this never gets hit, no new tails get added and the tail set becomes null at the end of while itr
                if (len(heads) > 0) and (heads.isdisjoint(unique_genes)):
                    for head in heads:
                        # error handling
                        if (cell_nw[tail, head] <= tolerance):
                            raise ValueError("Improper pair discovered.")

                        # update the next-level tails
                        new_tails.add(head)
                        heads_to_remember.add(head)
                        
                    # update the edge set for this hgraph level
                    tail_as_set = set()
                    tail_as_set.add(tail)
                    edges.append([tail_as_set, heads])
            
            for h in heads_to_remember:
                unique_genes.add(h)

            self.__merge_edges(edges, B, W, cell_nw)
            tails = set(deepcopy(new_tails))
        
        return B, W
    
    
    
    """
    Builds the hyper graph using the given tolerance as the metric by which to prune Transfer Entropy edges.
    
    Sets the full incidence matrix H and associated hyper edge weights W for the HyperGraph object.
    """
    def construct_graph(self, tolerance: float):
        hyper_edges = []
        weights = []
        for i in range(self.n_cells):
            G, W = self.__extract_hyper_edges(i, tolerance)
            if (len(G) > 0):
                hyper_edges.append(np.transpose(np.array(G)))
                weights.append(np.array(W)) 

        H = np.concatenate(hyper_edges, axis=1)
        W = np.concatenate(weights)
        if (np.sum(np.where(W < 0))):
            raise ValueError("Negative edge weights are not allowed.")

        self.hedges = hyper_edges
        self.W = W
        self.H = H

    """
    Updates H for this instance by removing vertices from the incidence matrix that are not present in any hyper edges.
    Assumes construct_graph routine has already been successfully run on this instance.
    """
    def prune_graph(self):
        count = 0
        to_prune = []

        for v in range(self.n_genes):
            indices = np.where(self.H[v,:] != 0)
            deg_v = np.sum(self.W[indices])
            if deg_v == 0:
                # remove the vertex, since it isn't involved in any hedges, therefore won't be in any pathways
                count += 1
                to_prune.append(v)

        self.n_genes = self.n_genes - count

        pruned_labels = np.delete(self.gene_labels, to_prune)
        pruned_H = np.delete(self.H, to_prune, axis=0) # delete rows with indices in to_prune

        # Checking it did its job
        if np.intersect1d(pruned_labels, self.gene_labels[to_prune]).shape[0] != 0:
            raise ValueError("Pruning of background genes failed.")

        if (pruned_H.shape[1] != self.H.shape[1]):
            raise ValueError("Bug in pruning, edges were deleted.")
        
        self.H = pruned_H
        self.gene_labels = pruned_labels

    """
    Sets the (normalized) graph Laplacian and associated degree matrices for this instance.
    Assumes construct_graph routine has already been successfully run on this instance.
    """
    def set_graph_Laplacian(self):
        self.Dv = np.zeros((self.n_genes, self.n_genes))
        self.De = np.zeros((self.H.shape[1], self.H.shape[1]))


        for v in range(self.n_genes):
            indices = np.where(self.H[v,:] != 0)
            deg_v = np.sum(self.W[indices])
            if deg_v == 0:
                raise ValueError("Pruning of background genes failed.")
            else:
                self.Dv[v, v] = deg_v 

        for e in range(self.H.shape[1]):
            deg_e = np.sum(self.H[:, e])  
            
            if deg_e == 0: 
                raise ValueError("Graph construction created an empty hyperedge.")
            else:
                self.De[e, e] = deg_e

        Dei = np.linalg.inv(self.De)
        Dvis = np.linalg.inv(np.sqrt(self.Dv)) # replace each diagonal entry with the reciprocal of its square root, i.e. di with 1/sqrt(di)
        W_diag = np.diag(self.W)

        self.L = (Dvis @ self.H @ W_diag @ Dei @ (self.H.T) @ Dvis) 

    def get_incidence_matrix(self):
        return self.H
    
    def get_hedge_weights(self):
        return self.W
    
    def get_Laplacian(self):
        return self.L
    
    """
    Returns the sorted eigenvalues and eigenvectors of the incidence matrix representing this Hypergraph
    """
    def get_sorted_eigens(self):
        eigvals, eigvecs = np.linalg.eig(self.L)
        # sort eigenvals, then sort eigenvecs by the eigenvals
        sorted_indices = np.flip(np.argsort(eigvals))
        sorted_eigvals = eigvals[sorted_indices] # sort the vals
        sorted_eigvecs = eigvecs[:,sorted_indices] # sorts the columns

        # Sanity check to ensure eigenvecs got sorted
        eps = 0.0001
        min_index = np.argmin(eigvals)
        if (np.sum(np.abs(np.abs(eigvecs[:,min_index]) - np.abs(sorted_eigvecs[:,self.n_genes - 1])) < eps) != self.n_genes):
            raise ValueError("Incorrectly sorted eigenvalues")
        
        return sorted_eigvals, sorted_eigvecs