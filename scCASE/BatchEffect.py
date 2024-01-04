from decimal import InvalidOperation
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csc_matrix,csr_matrix
import anndata as ad
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import non_negative_factorization
from numpy import linalg as LA
from sklearn.decomposition import NMF
from numpy import matrix
import os
import time
import episcanpy.api as epi
from sklearn.cluster import KMeans
from sklearn.preprocessing import Binarizer
from scipy.spatial.distance import pdist, jaccard
from scipy.spatial.distance import squareform
from kneed import KneeLocator
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.sparse import csr_matrix
import epiaster as aster
import numpy.matlib
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Normalizer
from scCASE.preprocessing import *

def correct_batch_effect(X2, K, S=0.95, Alpha=1, Lambda=1000000, Gamma=1, Gamma2=1, Theta=0, Stop_rule=1, Seed=100, W2=None, H=None, Z=None, R=None,saveZ = False):
    def J_compute(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta):
        norm1 = pow(LA.norm(np.dot(X2, Z*R) - np.dot(W2, H), ord = 'fro'), 2)
        norm2 = Lambda * pow(LA.norm(Z - np.dot(H.T, H), ord = 'fro'), 2)
        norm3 = Gamma * pow(LA.norm(W2, ord = 'fro'), 2)
        norm4 = Gamma2 *pow(LA.norm(H, ord = 'fro'), 2)
        obj = norm1 + norm2 + norm3 +  norm4
        return obj

    def W2_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj):
        D = -2 * np.einsum("ij,jk,jk,lk->il",X2,Z,R,H,optimize=True) + 2 * np.einsum("ij,jk,lk->il",W2,H,H,optimize=True)+ 2 * Gamma * W2
        EQ1 =  1 * np.einsum("ji,kj,kl,li,li",H,D,X2,Z,R,optimize=True)
        EQ2 =EQ1
        EQ3 = 1 * np.einsum("ji,kj,kl,li",H,W2,D,H,optimize=True)
        EQ4 =EQ3
        EQ5 = 1 * np.einsum("ji,kj,kl,li",H,D,D,H,optimize=True)
        delta = (0 - EQ1 - EQ2 + EQ3 + EQ4) / (2 * EQ5)
        W2 = W2 - delta * D
        W2 = np.maximum(W2, 0)
        return np.array(W2,dtype='float32'), obj

    
    def H_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj):
        obj = J_compute(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta)
        GradH = -2 * np.einsum("ji,jk,kl,kl->il",W2,X2,Z,R,optimize=True) + 2 * np.einsum("ji,jk,kl->il",W2,W2,H,optimize=True) - 2 * Lambda * np.dot(H, (Z+Z.T)) + 4 * Lambda * np.einsum("ij,kj,kl->il",H,H,H,optimize=True) + 2 * Gamma2 * H
        Hn = H
        for i in range(10):
            Hx = H - (0.2**(i)) * GradH
            Hn = (np.abs(Hx) + Hx) / 2.0
            objn = J_compute(X2, Z, Z0, R, W2, Hn, K, Lambda, Gamma, Gamma2, Theta)
            if objn <= obj:
                return np.array(Hn,dtype='float32'), objn
        return np.array(H,dtype='float32'), obj
    
    def Z_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj):
        D = 2 * np.einsum("ji,jk,kl,kl,il->il",X2,X2,Z,R,R, optimize=True)-2 * np.einsum("ji,jk,kl,il->il",X2,W2,H,R, optimize=True) + 2 * Lambda * Z - 2 * Lambda * np.dot(H.T, H)
        EQ1 =  2 * np.einsum("hi,ij,ij,hl,lj,lj",X2,D,R,X2,Z,R,optimize=True)
        EQ2 =  1 * np.einsum("hi,ij,ij,hl,lj,lj",X2,D,R,X2,D,R,optimize=True)
        EQ3 =  2 * np.einsum("ij,jk,im,mk,mk",W2,H,X2,D,R,optimize=True)
        EQ4 =  2 * np.einsum("ji,ji",Z,D,optimize=True)
        EQ5 =  1 * np.einsum("ji,ji",D,D,optimize=True)
        EQ6 =  2 * np.einsum("ij,ik,kj",H,H,D,optimize=True)
        delta = (EQ1 - EQ3 + Lambda * EQ4 - Lambda * EQ6) / (2 * EQ2 + 2 * Lambda * EQ5)
        Z = Z - delta * D
        Z = np.maximum(Z, 0)
        return np.array(Z,dtype='float32'), obj

 
    print("Initializing...")
    d2 = X2.shape
    n = d2[1]
    q = d2[0]
    
    if type(R) == type(None):
        np.random.seed(Seed)
        R = np.array(np.random.binomial(1, S, size=(n, n)))
    
    
    Maxiter_1 = 3
    nmf = NMF(n_components=K, max_iter=10000, random_state=Seed)
    
    W2 = csc_matrix(nmf.fit_transform(X2),dtype='float32')

    H = nmf.components_
    H = matrix(H)
    lib = np.sum(H, axis=1)
    H = H /(np.tile(lib, (1,n))+ np.spacing(1))
    H = csc_matrix(H,dtype='float32')
    
    Z0 = jaccard_dis(X2)
    Z = Z0
    
    if os.path.exists('Z0.npy'):
        print('Loading exist similarity matrix...')
        Z0 = np.load('Z0.npy')
        if Z0.shape[0] != R.shape[0]:
            print('wrong similarity matrix...')
            os.remove('Z0.npy')
            print('Generating similarity matrix...')
            Z0 = jaccard_dis(X2)
            if saveZ == True:
                np.save('Z0.npy', Z0)
    else:
        print('Generating similarity matrix...')
        Z0 = jaccard_dis(X2)
        if saveZ == True:
            np.save('Z0.npy', Z0)
            
    J_record = np.zeros(500)
    J_min = float('inf')
    print("Updating...")
    H= H.A
    X2 = X2.A
    W2 = W2.A
    obj = J_compute(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta)
    
    for iter in range(1, Maxiter_1+1):
        # update R
        '''
        Seed = Seed + 1
        np.random.seed(Seed)
        R = np.random.binomial(1, S, size=(n, n))
        R = np.array(R,dtype='float32')
        '''

        # normalize H
        H = matrix(H)
        lib = np.sum(H, axis=1)
        H = H /(np.tile(lib, (1,n))+ np.spacing(1))
        H = np.asarray(H,dtype='float32')
        HHt = np.dot(H, H.T)

        # update Z
        Z, obj = Z_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj)

        '''
        # normalize Z
        Z = matrix(Z)
        lib = np.sum(Z, axis=1)
        Z = Z /(np.tile(lib, (1,n))+ np.spacing(1))
        Z = np.asarray(Z)
        '''

        # update W2
        W2, obj = W2_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj)

        # update H
        H, obj = H_update(X2, Z, Z0, R, W2, H, K, Lambda, Gamma, Gamma2, Theta, obj)

        HHt = np.dot(H, H.T)

        norm1 = pow(LA.norm(np.dot(X2, Z*R) - np.dot(W2, H), ord = 'fro'), 2)
        norm2 = pow(LA.norm(Z - np.dot(H.T, H), ord = 'fro'), 2)
        norm3 = pow(LA.norm(W2, ord = 'fro'), 2)
        norm4 = pow(LA.norm(H, ord = 'fro'), 2)
        obj = norm1 + Lambda * norm2 + Gamma * norm3 + Gamma2 * norm4 
        J_record[iter-1] = obj
        #if(iter%1 == 0):
            #print('Iteration: '+str(iter)+'\t'+'J value: '+str(obj)+'\t'+'norm1: '+str(norm1)+'\t'+'norm2: '+str(norm2*Lambda)+'\t'+'norm3: '+str(norm3*Gamma)+'\t'+'norm4: '+str(norm4*Gamma2))
        iter = iter + 1
        if iter > 20:
            if (J_record[iter-10] - J_record[iter-9]) / (J_record[iter-9] + np.spacing(1)) < 1e-6:
                break
    
    print('Finished iteration!')

    

    
    return W2, H, Z