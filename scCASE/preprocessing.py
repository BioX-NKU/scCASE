from decimal import InvalidOperation
import numpy as np
import scanpy as sc
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import non_negative_factorization
from numpy import matrix
import episcanpy.api as epi
from sklearn.preprocessing import Binarizer
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.metrics import  silhouette_score
import epiaster as aster
from sklearn.preprocessing import MinMaxScaler

def jaccard_dis(data):
    data = np.array(data)
    data = binarizer(data).T
    res = pdist(data, 'jaccard')
    n = data.shape[0]
    jaccard_sim = (1 - squareform(res))
    jaccard_sim = matrix(jaccard_sim)
    lib = np.sum(jaccard_sim, axis=1)
    jaccard_sim = jaccard_sim/np.tile(lib, (1,n))
    jaccard_sim = np.asarray(jaccard_sim)
    return jaccard_sim

def evaluation_louvain(matrix, cluster):
    '''
    Calculate the louvain clustering silhouette
    '''
    adata = sc.AnnData(matrix)
    epi.pp.lazy(adata)
    epi.tl.louvain(adata, cluster)
    epi.tl.getNClusters(adata, n_cluster=cluster)
    # epi.pl.umap(adata, color=['louvain', 'cell_type'])
    silhouette_s = silhouette_score(matrix, adata.obs['louvain'], metric='correlation')
    return silhouette_s

def evaluation_leiden(matrix, cluster):
    '''
    Calculate the leiden clustering silhouette
    '''
    adata = sc.AnnData(matrix)
    epi.pp.lazy(adata)
    epi.tl.leiden(adata, cluster)
    epi.tl.getNClusters(adata, n_cluster=cluster, method='leiden')
    # epi.pl.umap(adata, color=['leiden', 'cell_type'])
    silhouette_s = silhouette_score(matrix, adata.obs['leiden'], metric='correlation')
    return silhouette_s


def Lambda_calculate(data_tfidf, method):
    """
    Determine Lambda empirically by cell number in data.
    """
    cell_num = data_tfidf.shape[1]
    if method == 'p':
#        Lambda = cell_num*1000
        if cell_num >= 800:
            Lambda = 1000000
        else:
            Lambda = 100000
    else:
        if cell_num >= 800:
            Lambda = 100000
        else:
            Lambda = 10000
    return Lambda

def Estimate_k(data_tfidf,method,search_range ):
    adata = sc.AnnData(data_tfidf.T,dtype="float32")
    search_list = list(search_range)
    estimated_k = aster.ensemble.estimate_k(adata, search_list)
    print(estimated_k)
    if method == "scCASE":
        estimated_k = estimated_k+1
    elif method =="scCASER":
         estimated_k = round(estimated_k*1.3) + estimated_k//2500

    return estimated_k

def select_peak(data, ref_mat, threshold = 0.01):
    select_peak = np.array((np.sum(data, axis=1) >= (data.shape[1] * threshold)).reshape(-1,1))[:,0]
    data_n = data[select_peak,:]
    print('Data shape after feature selection:')
    print(data_n.shape)

    if ref_mat is not None:
        if ref_mat.shape[0] != data.shape[0]:
            ref_mat = ref_mat.T
        ref_mat_n = ref_mat[select_peak,:]
        print('Reference shape after feature selection:')
        print(ref_mat_n.shape)
    else:
        ref_mat_n = ref_mat
    return data_n, ref_mat_n, select_peak

def tf_idf_transform(data):
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(data))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(data)))
    tf_idf = tf_idf.todense()
    tf_idf = np.array(tf_idf)
    return tf_idf

def tfidf(data,bulk = None):
    data_tfidf = tf_idf_transform(data)
    if bulk is not None:
        ref_tfidf =tf_idf_transform(bulk)
        #ref_tfidf = bulk
    else:
        ref_tfidf = None
    return data_tfidf,ref_tfidf


def binarizer(count_matrix,threshold = 0):
    bn = Binarizer(threshold = threshold)
    count_matrix_bn = bn.fit_transform(count_matrix)
    return count_matrix_bn

def myscaler(count_matrix,method = MinMaxScaler()):
    #MinMaxScaler() Normalizer() StandardScaler()
    method.fit(count_matrix)
    count_matrix_sc = method.transform(count_matrix)
    return count_matrix_sc

def nmf_W(X,fixed_W):
    fixed_H = fixed_W.T
    W,H,n_iter = non_negative_factorization(X.T, n_components=fixed_H.shape[0], init='custom', random_state=0, update_H=False, H=fixed_H)
    return W.T