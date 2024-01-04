#!/usr/bin/env python

import scipy
import numpy as np
import scanpy as sc
sc.settings.verbosity = 0

from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import silhouette_score
from sklearn.feature_extraction.text import TfidfTransformer
from .utils import *


def ssd_knee_est(adata, search_list, seed=None):
    """
    Estimate the number of cell types by sum of squared distances.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    search_list
        List of optional numbers of cell types for the estimation.
    seed
        Random seed for reproducibility. By default, `seed=None`.
        
    Returns
    -------
    sk_est
        Estimated number of cell types by sum of squared distances.

    """ 
    
    print('Estimating by sum of squared distances...')
    if seed is not None: setup_seed(seed)

    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(adata.X)
    model.idf_ -= 1
    tf_idf = scipy.sparse.csr_matrix(model.transform(adata.X))
    adata.X = tf_idf.copy()
        
    sc.pp.pca(adata, n_comps=50, svd_solver='arpack', use_highly_variable=False)
    
    distances = []
    for k in search_list:
        kmeanModel = KMeans(n_clusters=k)
        kmeanModel.fit(adata.obsm['X_pca'])
        distances.append(kmeanModel.inertia_)
    kl = KneeLocator(search_list, distances, S=1.0, curve="convex", direction="decreasing")
    sk_est = kl.knee
    
    return sk_est


def davies_bouldin_est(adata, search_list, seed=None):
    """
    Estimate the number of cell types by Davies-Bouldin score.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    search_list
        List of optional numbers of cell types for the estimation.
    seed
        Random seed for reproducibility. By default, `seed=None`.
        
    Returns
    -------
    db_est
        Estimated number of cell types by Davies-Bouldin score.

    """ 
    
    print('Estimating by Davies-Bouldin score...')
    if seed is not None: setup_seed(seed)

    count_mat = adata.X.T.copy()
    nfreqs = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    tfidf_mat = np.multiply(nfreqs, np.tile(np.log(1 + 1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1])))
    adata.X = scipy.sparse.csr_matrix(tfidf_mat).T
        
    sc.pp.pca(adata, n_comps=50, svd_solver='arpack', use_highly_variable=False)
    
    scores = []
    for k in search_list:
        kmeans = KMeans(n_clusters=k)
        model = kmeans.fit_predict(adata.obsm['X_pca'])
        score = davies_bouldin_score(adata.obsm['X_pca'], model)
        scores.append(score)
    db_est = search_list[np.argmin(scores)]
    
    return db_est


def silhouette_est(adata, search_list, seed=None):
    """
    Estimate the number of cell types by silhouette coefficient.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    search_list
        List of optional numbers of cell types for the estimation.
    seed
        Random seed for reproducibility. By default, `seed=None`.
        
    Returns
    -------
    sil_est
        Estimated number of cell types by silhouette coefficient.

    """ 
    
    print('Estimating by silhouette coefficient...')
    if seed is not None: setup_seed(seed)

    count_mat = adata.X.T.copy()
    nfreqs = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    tfidf_mat = np.multiply(nfreqs, np.tile(np.log(1 + 1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1])))
    adata.X = scipy.sparse.csr_matrix(tfidf_mat).T
        
    sc.pp.pca(adata, n_comps=50, svd_solver='arpack', use_highly_variable=False)
    sc.pp.neighbors(adata,  n_neighbors=15, n_pcs=50, method='umap', metric='euclidean')
    
    sil_louvain = []
    sil_leiden  = []
    for k in search_list:
        getNClusters(adata, n_cluster=k, method='louvain');
        sil_louvain.append(silhouette_score(adata.obsm['X_pca'], adata.obs['louvain'], metric='correlation'))
        getNClusters(adata, n_cluster=k, method='leiden');
        sil_leiden.append(silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'], metric='correlation'))  
    sil_est = search_list[np.argmax(np.array(sil_louvain) + np.array(sil_leiden))]
    
    return sil_est