
import pandas as pd
import scanpy as sc
import numpy as np
def load_data(data_path,ref_path = None,data_format = "count matrix",data_sep=",",ref_sep=",",var_key = None,obs_key = None):
    if data_format == "count matrix":
        cell_data =  pd.read_csv(data_path,index_col=0,sep=data_sep)
        bulk_data =  np.array(pd.read_csv(ref_path,index_col=None,header = None,sep=ref_sep),dtype='float32') if ref_path != None else None
        index_name = cell_data.index.values
        cell_name = cell_data.columns.values
        cell_data = np.array(cell_data,dtype='float32')
    if data_format == "h5ad":
        cell_data = sc.read(data_path)
        bulk_data = sc.read(ref_path).X.T if ref_path != None else None
        if obs_key == None:
            cell_name = np.arange(cell_data.obs.shape[0]).astype(str)
        else:
            cell_name = np.array(cell_data.obs[obs_key])
        if obs_key == None:
            index_name = np.arange(cell_data.var.shape[0]).astype(str)
        else:
            index_name = np.array(cell_data.var[var_key])
        if type(cell_data.X) != np.matrix:
            cell_data.X = cell_data.X.todense()
        cell_data = np.array(cell_data.X.T,dtype='float32')
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