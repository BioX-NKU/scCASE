#!/usr/bin/env python

import gc
import numpy as np
import episcanpy.api as epi
from .estimators import *


def estimate_k(adata, search_list, binary=True, fpeak=0.01, seed=2022):
    """
    Estimate the number of cell types in scCAS data by ASTER.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    search_list
        List of optional numbers of cell types for the estimation.
    binary
        Whether to convert the count matrix into a binary matrix. By default, `binary=True`.
    fpeak
        Select peaks/regions that have at least one read count in at least `fpeak` of the cells in the count matrix. By default, `fpeak=0.01`.
    seed
        Random seed for reproducibility. By default, `seed=None`.
        
    Returns
    -------
    estimated_k
        Estimated number of cell types.

    """
    print('Raw dataset shape: ', adata.shape)
    if binary: epi.pp.binarize(adata)
    epi.pp.filter_features(adata, min_cells=np.ceil(fpeak*adata.shape[0]))
    if np.sum(np.sum(adata.X, axis=1)==0) > 0: 
        print("There are empty cells after filtering features. These cells will be removed. Alternatively, use a smaller fpeak to avoid empty cells.")
        epi.pp.filter_cells(adata, min_features=1)
    print('Dataset shape after preprocessing: ', adata.shape)
    
    adata_sk = adata.copy()
    sk_est = ssd_knee_est(adata_sk, search_list, seed=seed)
    del adata_sk
    gc.collect()
    
    adata_db = adata.copy()
    db_est = davies_bouldin_est(adata_db, search_list, seed=seed)
    del adata_db
    gc.collect()
    
    adata_sil = adata.copy()
    sil_est = silhouette_est(adata_sil, search_list, seed=seed)
    del adata_sil
    gc.collect()
    
    est_sum = 0; cnt = 0
    if sk_est  is not None: est_sum += sk_est;  cnt += 1
    if db_est  is not None: est_sum += db_est;  cnt += 1
    if sil_est is not None: est_sum += sil_est; cnt += 1
    estimated_k = int(np.ceil(1.0*est_sum/cnt))
    
    return estimated_k