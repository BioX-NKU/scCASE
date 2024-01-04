#!/usr/bin/env python

import random
import numpy as np
import scanpy as sc


def setup_seed(seed):
    """
    Set random seed.
    
    Parameters
    ----------
    seed
        Number to be set as random seed for reproducibility.
        
    """
    np.random.seed(seed)
    random.seed(seed)

def getNClusters(adata, n_cluster, range_min=0, range_max=3, max_steps=20, method='louvain', key_added=None):
    """
    Tune the resolution parameter in Louvain or Leiden clustering to make the number of clusters and the specified number as close as possible.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    n_cluster
        Specified number of clusters.
    range_min
        Minimum clustering resolution for the binary search. By default, `range_min=0`.
    range_max
        Maximum clustering resolution for the binary search. By default, `range_max=3`.
    max_steps
        Maximum number of steps for the binary search. By default, `max_steps=20`.
    method
        Method (`louvain` or `leiden`) used for cell clustering. By default, `method='louvain'`.
    key_added
        The obs variable name of clustering results. By default, `key_added` is the same as the method name used.

    Returns
    -------
    adata
        AnnData object with clustering assignments in `adata.obs`:

        - `adata.obs['louvain']` - Louvain clustering assignments if `method='louvain'` and `key_added=None`.
        - `adata.obs['leiden']` - Leiden clustering assignments if `method='leiden'` and `key_added=None`.
    """
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    pre_min_cluster = 0
    pre_max_cluster = None
    min_update = False
    max_update = False
    while this_step < max_steps:
        this_step += 1
        if this_step==1: 
            this_resolution = this_min + ((this_max-this_min)/2)
        else:
            if max_update and not min_update:
                this_resolution = this_min + (this_max-this_min) * (n_cluster-pre_min_cluster)/(pre_max_cluster-pre_min_cluster)
            if min_update and not max_update:
                if pre_max_cluster is not None:
                    this_resolution = this_min + (this_max-this_min) * (n_cluster-pre_min_cluster)/(pre_max_cluster-pre_min_cluster)
                else:
                    this_resolution = this_min + ((this_max-this_min)/2)
                    
        if (method == 'louvain') and (key_added==None):
            sc.tl.louvain(adata, resolution=this_resolution)
        elif (method == 'louvain') and (type(key_added)==str):
            sc.tl.louvain(adata, resolution=this_resolution, key_added=key_added)
        elif (method == 'leiden') and (key_added==None):
            sc.tl.leiden(adata,resolution=this_resolution)
        elif (method == 'leiden') and (type(key_added)==str):
            sc.tl.leiden(adata,resolution=this_resolution, key_added=key_added)
        else:
            print('Error settings of method and key_added.')
        
        if key_added==None:
            this_clusters = adata.obs[method].nunique()
        else:
            this_clusters = adata.obs[key_added].nunique()
                
        if this_clusters > n_cluster:
            this_max = this_resolution
            pre_max_cluster = this_clusters
            min_update = False
            max_update = True
        elif this_clusters < n_cluster:
            this_min = this_resolution
            pre_min_cluster = this_clusters
            min_update = True
            max_update = False
        elif this_clusters == n_cluster:
            break