import pickle

# TODO: add clustering routines and visuals

"""
Load in a pickled object. 
"""
def load_pickle(save_path: str):
    with open(save_path, 'rb') as instream:
        obj = pickle.load(instream) 
    
    return obj