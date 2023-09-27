from scCASE.preprocessing import *
from scCASE.loaddata import *
from scCASE.model import *
def run(data_path,ref_path = None,method = "scCASE",data_format = "count matrix",
         data_sep=",",ref_sep=",",type_number_range=range(3, 15),output_path = "./",
         save_other_matrixs_ = False,obs_key = None,var_key = None,threshold = 0.01,
         saveZ = False,changeK = "+0",changeK_ref = "+0",save_result = False):




    if method == "scCASE":
        raw_data, raw_ref,*label = load_data(data_path ,None,data_format,data_sep,ref_sep,var_key,obs_key)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        label[1] = np.array(pd.DataFrame(label[1])[0].str.replace(r'\.\d*','',regex =True))
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        K1 = Estimate_k(raw_data,method,type_number_range)
        K1 = eval("K1"+changeK)
        Lambda = Lambda_calculate(raw_data, "p")
        X_imputed,*Matrices= run_scAGp(raw_data,K1,Lambda = Lambda,saveZ = saveZ)

    
    elif method == "scCASER":
        raw_data, raw_ref,*label = load_data(data_path,ref_path,data_format,data_sep,ref_sep,var_key,obs_key)
        raw_data, raw_ref, peak_name = select_peak(raw_data, raw_ref, threshold )
        label[0] = label[0][peak_name]
        label[1] = np.array(pd.DataFrame(label[1])[0].str.replace(r'\.\d*','',regex =True))
        raw_data, raw_ref = tfidf(raw_data,raw_ref)
        K1 = Estimate_k(raw_data,method,type_number_range)
        K1 = eval("K1"+changeK)
        K1 = round(K1)
        K2 = round(K1/3)
        K2 = eval("K2"+changeK_ref)
        K2 = round(K2)
        Lambda = Lambda_calculate(raw_data, "r")
        X_imputed,*Matrices= run_scAGr(raw_data, raw_ref,K2,K1,Lambda = Lambda,saveZ = saveZ)

    else:
        print("Wrong method name")
        return 0
    
    X_imputed = myscaler(X_imputed.T).T
    X_imputed = myscaler(X_imputed)
    X_imputed = pd.DataFrame(X_imputed,columns = label[1],index = label[0])
    if save_result:    
        save_data(output_path+"result.csv",X_imputed,label)
    if save_other_matrixs_ == True:
        save_other_matrixs(output_path,label,Matrices,save_W = True, save_H = True,save_Z = True)
        return X_imputed,*Matrices
    return X_imputed