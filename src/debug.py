import anndata as ad
import numpy as np
import pandas as pd
from copy import deepcopy

n_cells = 456
n_genes = 100
PWscores = np.zeros((n_cells, n_genes, n_genes))

for i in range(n_cells):
    path = "./data/TEsmESC/geneXgene{i}.csv".format(i=(i+1))
    df_temp = pd.read_csv(path)
    PWscores[i, :, :] = np.copy(df_temp.values)
    


"""
Merges tails with common heads to identify all the multi-arity relations to represent as hedges and updates the supplied incidence matrix with these final hedges.
_edges: a list of tuples (tail_set, head_set) of all the initial single-tailed edges (note: they are sets, not lists!)
B: a list of n_gene-sized integer numpy arrays to represent the incidence matrix to update.
TODO: Test and be especially make sure Python's memory system isn't doing any weird reference copying and producing the wrong stuff!
Note, all set operations return copies, so this SHOULD be okay.
"""
def merge_edges(_edges, B):
    
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
                    new_tails = e1[0].intersection(e2[0]) # should be the same as doing a union of the tails
                    new_edge = ([new_tails, common_heads])
                    new_edges.append(new_edge)

                    # Updating old edges (their head sets)
                    e1[1] = set(e1[1] - common_heads)
                    e2[1] = set(e2[1] - common_heads)

        
        # Add all edges with non-null head sets to the incidence matrix!
        for edge in edges:
            col = np.zeros(n_genes)
            if len(edge[1]) > 0:
                # TODO: handle weights... how? Each hedge gets assigned a single weight,and those get stored in a different matrix/use a function. 
                for tail in edge[0]:
                    col[tail] = -1
                for head in edge[1]:
                    col[head] = 1

                B.append(np.copy(col))
        
        edges = deepcopy(new_edges)
        m = len(edges)


""" 
Creates a hypergraph, represented using an incidence matrix, from the given cell's PW gene network. 
cell_index: the index into the PWscores array for the desired cell.
p: percentage parameter to use in tolerance calculation; decides which relationships relative to the maximum PW TE. 
returns: B the incidence matrix of the constructed hypergraph. 
"""
def construct_hyper_graph(cell_index, p):
    # init
    unique_genes = set()
    cell_nw = PWscores[cell_index, :, :]

    # The incidence matrix (might use different data structure)
    B = []

    # all scores within whatever % of the max
    tolerance = np.max(cell_nw) * p

    tails = set(np.arange(start=0, stop=n_genes, step=1))

    # debug var
    count = 0

    # TODO: TEST IT ALL 
    while (len(tails) > 0):
        edges = []
        new_tails = set()

        for tail in tails:
            heads = set(np.flatnonzero(cell_nw[tail, :] > tolerance))

            # TODO: This conditional is something to test, in general 
            # Termination condition: once there's no new heads added this never gets hit, no new tails get added and the tail set becomes null at the end of while itr
            if (len(heads) > 0) and (heads.isdisjoint(unique_genes)):
                
                count += 1
                
                for head in heads:
                    # error handling
                    if (cell_nw[tail, head] <= tolerance):
                        print("BUG") # TODO: keep this error cond here and handle it better!

                    # update the unique heads set and next-level tails
                    unique_genes.add(head)
                    new_tails.add(head)
                
                # update the edge set for this hgraph level
                tail_as_set = set()
                tail_as_set.add(tail)
                edges.append([tail_as_set, heads])

        merge_edges(edges, B)
        tails = set(deepcopy(new_tails))
    
    return B


G = construct_hyper_graph(0, 0.5)