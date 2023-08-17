#!/usr/bin/env python
import argparse
from scCASE import run

if __name__ == '__main__':
    def parse_args():
        parser = argparse.ArgumentParser(description='Parameters')
        parser.add_argument('--method', type=str, default='scCASE', help='scCASE or scCASER.')
        parser.add_argument('--data_path', type=str,default=None, help='The path and file name of the target data.')
        parser.add_argument('--ref_path', type=str, default=None, help='The path and file name of the reference data.')
        parser.add_argument('--data_format', type=str, default="count matrix", help='The format of data and reference data, including count matrix(csv/txt/tsv) and h5ad. Default:count matrix')
        parser.add_argument('--data_sep', type=str, default=",", help='The Separator of target data, only for count matrix file. Default:,')
        parser.add_argument('--ref_sep', type=str, default=",", help='The Separator of reference data, only for count matrix file. Default:,')
        parser.add_argument('--type_number_range', type=range, default=range(3, 15), help='The range of possible number of cell types. Default:(2:15)')
        parser.add_argument('--output_path', type=str, default="./", help='The path of the result. Defalut:./')
        parser.add_argument('--save_other_matrixs', type=bool, default=False, help='If save the matrices including Z,W and H. Default:False')
        parser.add_argument('--obs_key', type=str, default=None, help='The key of cell name in adata.obs. Default:celltype')
        parser.add_argument('--var_key', type=str, default=None, help='The key of peaks name in adata.var. Default:peaks')
        parser.add_argument('--threshold', type=float, default=0.01 , help='The threshold of preprocessing. Default:0.01')
        parser.add_argument('--saveZ', type=bool, default=False , help='If save the initialized matrix Z. If you need to use this method multiple times for the same data, with different parameters, please select True, which can significantly accelerate the speed')
        parser.add_argument('--changeK', type=str, default="+0" , help='Change the parameter K.')
        parser.add_argument('--changeK_ref', type=str, default="+0" , help='Change the parameter K of reference.')
        return parser.parse_args()
    args = parse_args()
    method = args.method
    data_path = args.data_path
    ref_path = args.ref_path
    data_format = args.data_format
    data_sep = args.data_sep
    ref_sep = args.ref_sep
    type_number_range = args.type_number_range
    output_path = args.output_path
    save_other_matrixs_ = True #save_other_matrixs_
    var_key = args.var_key
    obs_key = args.obs_key
    threshold = args.threshold
    saveZ = args.saveZ
    changeK = args.changeK
    changeK_ref = args.changeK_ref
    result = run(data_path,ref_path,method,data_format,data_sep,ref_sep,type_number_range,output_path,save_other_matrixs_,obs_key,var_key,threshold,saveZ,changeK,changeK_ref,True)