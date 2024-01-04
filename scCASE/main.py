from scCASE.preprocessing import *
from scCASE.loaddata import *
from scCASE.model import *
from scCASE.SeqDepth import *
from scCASE.BatchEffect import *
from scCASE.analysis import *

def run(data_path,ref_path = None,method = "scCASE",data_format = "h5ad",
         data_sep=",",ref_sep=",",type_number_range=range(3, 15),output_path = "./",
         batchlabel= "batch",threshold = 0.01,
         saveZ = False,K = "Auto",K_ref = "Auto",save_result = False):

    if method == "scCASE":
        raw_data, raw_ref,*label = load_data(data_path ,None,data_format,data_sep,ref_sep)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        if K == "Auto":
            K1 = Estimate_k(raw_data.A,method,type_number_range)
        else:
            K1 = int(K)

        Lambda = Lambda_calculate(raw_data, "p")
        X_imputed,*Matrices= run_scCASE(raw_data,K1,Lambda = Lambda,saveZ = saveZ)
        X_imputed = myscaler(X_imputed.T).T
        X_imputed = myscaler(X_imputed)
        X_imputed = pd.DataFrame(X_imputed)
        X_imputed = ad.AnnData(X_imputed.T)
        X_imputed.var = label[0]
        X_imputed.obs = label[1]
        X_imputed.varm["Projection"] = Matrices[0]
        X_imputed.obsm["Embeding"] = Matrices[1].T
        X_imputed.obsp["Similarity"] = Matrices[2]


    elif method == "Correct batch effect":
        raw_data, raw_ref,*label = load_data(data_path ,None,data_format,data_sep,ref_sep)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        if K == "Auto":
            K1 = Estimate_k(raw_data.A,method,type_number_range)
        else:
            K1 = int(K)
        Lambda = Lambda_calculate(raw_data, "p")
        
        
        
        Matrices= correct_batch_effect(raw_data,K1,Lambda = Lambda,saveZ = saveZ)
        
        batchname = label[1][batchlabel]
        batchset = label[1][batchlabel].unique()
        Z_fix = np.ones([label[1].shape[0], label[1].shape[0]])
        Z = Matrices[2]
        for i in range(label[1].shape[0]):    
            for j in batchset:

                a = np.mean(Z.T[i][batchname == j])
                b = np.mean(Z.T[i][~(batchname == j)])
                Z_fix[batchname == j,i] = b/(a+b)
                
        Z = Z*Z_fix
        X_imputed = np.dot(raw_data.A,Z)
        X_imputed = myscaler(X_imputed.T).T
        X_imputed = myscaler(X_imputed)
        X_imputed = pd.DataFrame(X_imputed)
        X_imputed = ad.AnnData(X_imputed.T)
        X_imputed.var = label[0]
        X_imputed.obs = label[1]
        X_imputed.varm["Projection"] = Matrices[0]
        X_imputed.obsm["Embeding"] = Matrices[1].T
        X_imputed.obsp["Similarity"] = Matrices[2]    
    
    
    elif method == "scCASER":
        raw_data, raw_ref,*label = load_data(data_path,ref_path,data_format,data_sep,ref_sep)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        if (K == "Auto") & (K_ref == "Auto"):
            K1 = Estimate_k(raw_data.A,method,type_number_range)
            K1 = round(K1)
            K2 = round(K1/3)
            K2 = round(K2)
        elif (K != "Auto") & (K_ref == "Auto"):
            K1 = int(K)
            K2 = round(K1/3)
            K2 = round(K2)
        else:  
            K1 = int(K)
            K2 = int(K_ref)
        Lambda = Lambda_calculate(raw_data, "r")
        X_imputed,*Matrices= run_scCASER(raw_data, raw_ref,K2,K1,Lambda = Lambda,saveZ = saveZ)
        X_imputed = myscaler(X_imputed.T).T
        X_imputed = myscaler(X_imputed)
        X_imputed = pd.DataFrame(X_imputed)
        X_imputed = ad.AnnData(X_imputed.T)
        X_imputed.var = label[0]
        X_imputed.obs = label[1]
        X_imputed.obsp["Similarity"] = Matrices[0]
        
    elif method == "Correct seq depth":
        raw_data, raw_ref,*label = load_data(data_path,ref_path,data_format,data_sep,ref_sep)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        depth = np.sum(raw_data,axis = 0)
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        if K == "Auto":
            K1 = Estimate_k(raw_data.A,method,type_number_range)
        else:
            K1 = int(K)
        depth = np.diag(np.array(depth)[0])
        Lambda = Lambda_calculate(raw_data, "p")
        X_imputed,*Matrices= correct_read_depth(raw_data,K1,Lambda = Lambda,depth = depth,saveZ = saveZ)
        X_imputed = myscaler(X_imputed.T).T
        X_imputed = myscaler(X_imputed)
        X_imputed = pd.DataFrame(X_imputed)
        X_imputed = ad.AnnData(X_imputed.T)
        X_imputed.var = label[0]
        X_imputed.obs = label[1]
        X_imputed.varm["Projection"] = Matrices[0]
        X_imputed.obsm["Embeding"] = Matrices[1].T
        X_imputed.obsp["Similarity"] = Matrices[2]
  
    else:
        print("Wrong method name!")
        return 0

    if save_result:    
        X_imputed.write(output_path)
    return X_imputed