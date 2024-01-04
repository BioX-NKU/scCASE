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

def run_scCASER(X2, X1, K1 , K2 , S=0.95, Alpha=0, Lambda=3, Gamma=1, Gamma2=1, Theta=0, Stop_rule=1, Seed=100, W1=None, W2=None, H=None, Z=None, R=None,saveZ = False):
    def J_compute(X2, P1, Z, Z0, R, W1, W2, H, Alpha, Lambda, Gamma, Gamma2, Theta):
        Wm = W1 * M + W2 * N
        norm1 = pow(LA.norm(np.dot(X2, Z*R) - np.dot(Wm, H), ord = 'fro'), 2)
        norm2 = Lambda * pow(LA.norm(Z - np.dot(H.T, H), ord = 'fro'), 2)
        norm3 = Gamma *pow(LA.norm(W2*N, ord = 'fro'), 2)
        norm4 = Gamma2 *pow(LA.norm(H, ord = 'fro'), 2)
        norm5 = Theta *pow(LA.norm(Z - Z0, ord = 'fro'), 2)
        norm6 = Alpha * pow(LA.norm(P1-W1*M, ord = 'fro'), 2)
        obj = norm1+norm2+norm3 +norm4 + norm5 +norm6
        return obj

    def W1_update(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta, obj):
        Wm = W1 * M + W2 * N
        ZR = Z*R
        ZRHt = np.dot(ZR, H.T)
        X2ZRHt = np.dot(X2, ZRHt)
        WmHHt = np.dot(Wm, HHt)
        W1n = W1
        D1 = M * (-2 * X2ZRHt + 2 * WmHHt - 2 * Alpha * P1 + 2 * Alpha * M * W1);
        # eq1
        eq1 = np.dot(H.T, (D1*M).T)
        eq1 = np.dot(eq1, X2)
        eq1 = np.trace(np.dot(eq1, ZR))
        # eq2
        eq2 = np.dot(ZR.T, X2.T)
        eq2 = np.dot(eq2, D1*M)
        eq2 = np.trace(np.dot(eq2, H))
        # eq3
        W2D = np.dot((W2*N).T, D1*M)
        eq3 = np.dot(H.T, W2D+W2D.T)
        eq3 = np.trace(np.dot(eq3, H))
        # eq4
        DW1 = np.dot((D1*M).T, W1*M)
        eq4 = np.dot(H.T, DW1+DW1.T)
        eq4 = np.trace(np.dot(eq4, H))
        # eq5
        eq5 = np.dot((D1*M).T, D1*M)
        eq5 = np.dot(H.T, eq5)
        eq5 = np.trace(np.dot(eq5, H))
        # eq6
        PMD = np.dot(P1.T, D1*M)
        eq6 = np.trace(PMD-PMD.T)
        # eq7
        eq7 = np.trace(DW1+DW1.T)
        # eq8
        eq8 = np.trace(np.dot((D1*M).T, D1*M))

        #step_lr = (0-eq1-eq2+eq3+eq4-Alpha*eq6+Alpha*eq7) / (2*eq5+2*Alpha*eq8)
        step_lr = (0-eq1-eq2+eq3+eq4-(Alpha+Gamma)*eq6+(Alpha+Gamma)*eq7) / (2*eq5+2*(Alpha+Gamma)*eq8)
        W1 = W1 - step_lr * D1
        W1 = np.maximum(W1, 0)
        return W1, obj
    
    
    def W2_update(X2, P1, Z, Z0, R, W1, W2, H, Alpha, Lambda, Gamma, Gamma2, Theta, obj):
        Wm = W1 * M + W2 * N
        ZR = Z*R
        ZRHt = np.dot(ZR, H.T)
        X2ZRHt = np.dot(X2, ZRHt)
        WmHHt = np.dot(Wm, HHt)
        W2n = W2
        D2 = (-2 * X2ZRHt + 2 * WmHHt + 2 * Gamma * W2) * N;

        # eq1
        eq1 = np.dot(H.T, (D2*N).T)
        eq1 = np.dot(eq1, X2)
        eq1 = np.trace(np.dot(eq1, ZR))

        # eq2
        eq2 = np.dot(ZR.T, X2.T)
        eq2 = np.dot(eq2, D2*N)
        eq2 = np.trace(np.dot(eq2, H))

        # eq3
        W1D2 = np.dot((W1*M).T, D2*N)
        eq3 = np.dot(H.T, W1D2+W1D2.T)
        eq3 = np.trace(np.dot(eq3, H))

        # eq4
        D2W2 = np.dot((D2*N).T, W2*N)
        eq4 = np.dot(H.T, D2W2+D2W2.T)
        eq4 = np.trace(np.dot(eq4, H))

        # eq5
        eq5 = np.dot((D2*N).T, D2*N)
        eq5 = np.dot(H.T, eq5)
        eq5 = np.trace(np.dot(eq5, H))

         # eq6
        eq6 = np.dot((D2*N).T, D2*N)
        eq6 = np.trace(eq6)

        #eq7
        eq7 =  np.trace(D2W2)

        step_lr = (0-eq1-eq2+eq3+eq4+2*Gamma*eq7) / (2*eq5 + 2*Gamma*eq6)
        #step_lr = (0-eq1-eq2+eq3+eq4) / 2.0
        #step_lr = (0-eq1-eq2+eq3+eq4) / (2.0*eq5)
        W2 = W2 - step_lr * D2
        W2 = np.maximum(W2, 0)
        return W2, obj
    
    def H_update(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta, obj):
        obj = J_compute(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta)
        Wm = W1 * M + W2 * N
        ZR = Z*R
        WmtX2 = np.dot(Wm.T, X2)
        WmtX2ZR = np.dot(WmtX2, ZR)
        HZZt = np.dot(H, (Z+Z.T))
        HHT = np.dot(H,H.T)
        HHTH = np.dot(HHT, H)
        WmtWm = np.dot(Wm.T, Wm)
        WTWH = np.dot(WmtWm, H)
        GradH = -2 * WmtX2ZR + 2 * WTWH - 2 * Lambda * HZZt + 4 * Lambda * HHTH + 2 * Gamma2 * H
        Hn = H
        for i in range(10):
            Hx = H - (0.2**(i)) * GradH
            Hn = (np.abs(Hx) + Hx) / 2.0
            objn = J_compute(X2, P1, Z, Z0, R, W1, W2, Hn, Alpha, Lambda, Gamma, Gamma2, Theta)
            if objn <= obj:
                return Hn, objn
        return H, obj
    
    def Z_update(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta, obj):
        Wm = W1 * M + W2 * N
        ZR = Z*R
        WmtX2 = np.dot(Wm.T, X2)
        HtH = np.dot(H.T, H)
        X2tWmH = np.dot(WmtX2.T, H)
        RX2tWmH = R * X2tWmH
        X2tX2ZR = np.dot(X2tX2, ZR)
        X2tX2ZRR = X2tX2ZR * R
        D4 = 2 * X2tX2ZRR - 2 * RX2tWmH + 2 * Lambda * Z - 2 * Lambda * HtH + 2 * Theta * Z - 2 * Theta * Z0
        # eq1
        eq1 = np.dot((D4*R).T, X2.T)
        eq1 = np.dot(eq1, X2)
        eq1 = np.dot(eq1, ZR)
        eq1 = np.trace(eq1 + eq1.T)

        # eq2
        D4RX = np.dot((D4*R).T, X2.T)
        eq2 = np.trace(np.dot(D4RX, D4RX.T))

        # eq3
        HWXDR = np.dot(X2tWmH.T, D4*R)
        eq3 = np.trace(HWXDR + HWXDR.T)

        # eq4
        ZD4 = np.dot(Z.T, D4)
        eq4 = np.trace(ZD4 + ZD4.T)

        # eq5
        eq5 = np.trace(np.dot(D4.T, D4))

        # eq6
        HtHD4 = np.dot(HtH, D4)
        eq6 = np.trace(HtHD4 + HtHD4.T)

        step_lr = (eq1 - eq3 + Lambda * eq4 - Lambda * eq6) / (2 * eq2 + 2 * Lambda * eq5)

        Z = Z - step_lr * D4
        Z = np.maximum(Z, 0)
        return Z, obj
 
   
    print("Initializing...")
    d2 = X2.shape
    d3 = X1.shape
    n = d2[1]
    q = d2[0]
    m = d3[1]
    print('m = '+str(m))
    print('n = '+str(n))
    print('q = '+str(q))
    if type(W1) == type(None):
        np.random.seed(Seed)
        W1 = np.matlib.rand(q, K1+K2)

    if type(W2) == type(None):
        np.random.seed(Seed)
        W2 = np.matlib.rand(q, K1+K2)

    if type(H) == type(None):
        np.random.seed(Seed)
        H = np.matlib.rand(K1+K2, n)

    if type(Z) == type(None):
        np.random.seed(Seed)
        Z = np.matlib.rand(n, n)

    if type(R) == type(None):
        np.random.seed(Seed)
        R = np.random.binomial(1, S, size=(n, n))

    nmf = NMF(n_components=K1, max_iter=1000, random_state=100)
    P = nmf.fit_transform(X1)     
        
    M1 = np.ones((q, K1))
    M2 = np.zeros((q, K2))
    M = np.hstack((M1,M2)) # mask W1
    N = np.ones((q,K1+K2)) - M # mask W2
    M3 = np.zeros((q, K2))
    M4 = np.zeros((q, K1))
    P1 = np.hstack((P, M3)) # Full P
    
    Maxiter_1 = 5
    X2tX2 = np.dot(X2.T, X2).A
    # Z = np.eye(n)
    Z = np.asarray(Z)
    
    nmf = NMF(n_components=K1+K2, max_iter=10000, random_state=Seed)
    W2 = nmf.fit_transform(X2)               
    H = nmf.components_                         
    W1 = P1
        
    H = matrix(H)
    lib = np.sum(H, axis=1)
    H = H /(np.tile(lib, (1,n))+ np.spacing(1))
    H = np.asarray(H)
        
    HHt = np.dot(H, H.T)
    
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

        
    Z = Z0
    Z = np.asarray(Z)
    
    J_record = np.zeros(500)
    J_min = float('inf')
    print("Updating...")
    X1 = X1.A
    X2 = X2.A

    obj = J_compute(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta)
    for iter in range(1, Maxiter_1):
               
        
        # normalize H
        H = matrix(H)
        lib = np.sum(H, axis=1)
        H = H /(np.tile(lib, (1,n))+ np.spacing(1))
        H = np.asarray(H)
        
        HHt = np.dot(H, H.T)
        
        
        # update Z

        Z, obj = Z_update(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta, obj)
        
        
        # update W2
        W2, obj = W2_update(X2, P1, Z, Z0, R, W1, W2, H,  Alpha, Lambda, Gamma, Gamma2, Theta, obj)
        
         # update H
        H, obj = H_update(X2, P1, Z, Z0, R, W1, W2, H, Alpha, Lambda, Gamma, Gamma2, Theta, obj)
        
        HHt = np.dot(H, H.T)
        Wm = W1 * M + W2 * N
        norm1 = pow(LA.norm(np.dot(X2, Z*R) - np.dot(W2, H), ord = 'fro'), 2)
        norm2 = pow(LA.norm(Z - np.dot(H.T, H), ord = 'fro'), 2)
        norm3 = pow(LA.norm(W2*N, ord = 'fro'), 2)
        norm4 = pow(LA.norm(H, ord = 'fro'), 2)
        norm5 = pow(LA.norm(Z - Z0, ord = 'fro'), 2)
        norm6 = pow(LA.norm(P1-W1*M, ord = 'fro'), 2)
        obj = norm1 + Lambda * norm2 + Gamma * norm3 + Gamma2 * norm4 + Theta * norm5 + Alpha * norm6
        #print('Iteration: '+str(iter)+'\t'+'J value: '+str(obj)+'\t'+'norm1: '+str(norm1)+'\t'+'norm2: '+str(norm2*Lambda)+'\t'+'norm3: '+str(norm3*Gamma) + '\t' + 'norm4: '+str(norm4*Gamma2) + '\t' + 'norm5: '+str(norm5*Theta) + '\t' + 'norm6: '+str(norm6*Alpha) )
    print('Finished')

    
    return np.dot(X2, Z),Z

def run_scCASE(X2, K, S=0.95, Alpha=1, Lambda=1000000, Gamma=1, Gamma2=1, Theta=0, Stop_rule=1, Seed=100, W2=None, H=None, Z=None, R=None,saveZ = False):
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
    
    print('Finished')


    return np.dot(X2, Z), W2, H, Z