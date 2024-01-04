import episcanpy as epi
import numpy as np
import anndata as ad
from scipy.stats import wilcoxon
import pandas as pd


def tfidf_signac(data,inplace = True) : 
    try:
        count_mat = data.X.A.T
    except:
        count_mat = data.X.T

    tf_mat = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    tfidf_mat = np.log(1 + np.multiply(1e4*tf_mat, np.tile((1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1]))))
    tfidf_mat = np.nan_to_num(tfidf_mat)
    if inplace == True:
        data.X = tfidf_mat.T
        return data
    else: 
        _ = data.copy()
        _.X = tfidf_mat.T
        return _
    
def lazy(data):
    data = tfidf_signac(data)
    epi.pp.pca(data)
    epi.pp.neighbors(data)
    epi.tl.umap(data)
    epi.tl.tsne(data)
    
def identify_specific_peaks(data_enhanced,type_ ,obs = "cell_type",peak_name = "peak"):
    Projection_matrix = data_enhanced.varm["Projection"]
    Cell_embedding = ad.AnnData(data_enhanced.obsm["Embeding"])
    Cell_embedding.obs = data_enhanced.obs
    CompareCells = Cell_embedding[Cell_embedding.obs[obs] == type_]
    mean_value = np.mean(CompareCells.X,axis = 1)
    specific = pd.DataFrame([[wilcoxon(CompareCells.X[:,x],mean_value, alternative='greater').statistic,wilcoxon(CompareCells.X[:,x],mean_value, alternative='greater').pvalue] for x in range(Cell_embedding.shape[1])],columns = ["statistic","pvalue"])
    patterns_row = np.argmax(specific["statistic"])
    if specific.iloc[patterns_row,1]>0.05:
        print("Select cells may have not certain chromatin accessibility patterns.")
    patterns_value = Projection_matrix[:,patterns_row]
    patterns = pd.DataFrame(patterns_value,data_enhanced.var[peak_name].values,columns = ["Value"] )
    patterns.sort_values("Value",ascending=False,inplace = True)
    patterns = patterns/patterns.iloc[0,0]
    return patterns