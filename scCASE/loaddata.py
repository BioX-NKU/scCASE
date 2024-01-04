
import pandas as pd
import scanpy as sc
import numpy as np
import anndata
def load_data(data_path,ref_path = None,data_format = "h5ad",data_sep=",",ref_sep=","):
    if data_format == "count matrix":
        cell_data =  pd.read_csv(data_path,index_col=0,sep=data_sep)
        bulk_data =  np.array(pd.read_csv(ref_path,index_col=None,header = None,sep=ref_sep),dtype='float32') if ref_path != None else None
        index_name = cell_data.index.values
        cell_name = cell_data.columns.values
        cell_data = np.array(cell_data,dtype='float32')
    if data_format == "h5ad":
        if type(data_path) == str:
            cell_data = sc.read(data_path)
            cell_name = cell_data.obs
            index_name = cell_data.var
            cell_data = cell_data.X.T
        elif type(data_path) == anndata._core.anndata.AnnData:
            cell_data = data_path
            cell_name = cell_data.obs
            index_name = cell_data.var
            cell_data = cell_data.X.T    
    bulk_data = np.array(pd.read_csv(ref_path,index_col=None,header = None,sep=ref_sep),dtype='float32') if ref_path != None else None
        
    return cell_data,bulk_data,index_name,cell_name





def save_data(path,X,label,sep = ","):
    X_imputed = pd.DataFrame(X,columns = label[1],index = label[0])
    X_imputed.to_csv(path,sep = ",")
    
def save_other_matrixs(path,label,matrixs,save_W = True, save_H = True,save_Z = True):
    if save_W:
        W = matrixs[0]
        W = pd.DataFrame(W,index = label[0],columns = None)
        W.to_csv(path+"/scCASE_W.csv",sep = ",")
    if save_H:
        H = matrixs[1]
        H = pd.DataFrame(H,index = None,columns = label[1])
        H.to_csv(path+"/scCASE_H.csv",sep = ",",index = None)
    if save_Z:
        Z = matrixs[2]
        Z = pd.DataFrame(Z,columns = label[1],index = label[1])
        Z.to_csv(path+"/scCASE_Z.csv",sep = ",")