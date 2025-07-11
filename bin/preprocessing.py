#!/usr/bin/env python
# coding: utf-8

# problem in scanpy. solution: https://github.com/numba/numba/issues/4032 and https://github.com/dpeerlab/SEACells/issues/38
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'

import scanpy as sc
import scib
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import isspmatrix_csr

# import os

warnings.filterwarnings('ignore')


def run_preprocessing(input, output, batch, hvg, normalize, scale):
    """
    params:
        input: path of the anndata object
        output: path of the preprocessed file to be written
        batch: variable in cells metadata to be used as batch
        hvg: number of highly variable genes to use
        scale: set to true to activate scaling
    """

    adata = sc.read(input)
    # hvgs = adata.var.index

    # remove HVG if already precomputed
    if 'highly_variable' in adata.var:
        del adata.var['highly_variable']

    ## check the behaviour of this part
    print("Computing HVGs ...")
    hgv_total = find_highly_var_genes_scanpy(
        adata,
        total_hvgs = hvg,
        batch = batch,
        output_file = os.path.join(output, 'total-hvgs.csv'),
        verbose = True
    )
    adata = adata[:, hgv_total]

    adata.layers['counts'] = adata.X.copy()

    if normalize: 
        print("Normalizing and log-transforming data ...")
        sc.pp.normalize_total(adata, target_sum = 1e4)
        sc.pp.log1p(adata)

    if scale:
        print("Scaling data ...")
        adata = scib.preprocessing.scale_batch(adata, batch)

    # if seurat:
    #     print("Save as RDS")
    #     scib.preprocessing.saveSeurat(adata, output, batch, hvg)

    if not isspmatrix_csr(adata.X):
        print("X data format is incorrect. Changing sparse class to csr:")
        adata.X = adata.X.tocsr()


    print("Save as HDF5")
    sc.write(os.path.join(output, 'adata-preprocessed.h5ad'), adata)
 

# same function as in NMFusion
def find_highly_var_genes_scanpy(
    adata,
    total_hvgs = 2500,
    batch = None,
    output_file = None,
    verbose = True
):
    """
    Exported function in case user wants to do it step by step
    
    Compute highly variable genes for each batch and return the union of them. 
    """

    # computing HVGs
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    hvgs = sc.pp.highly_variable_genes(
        adata, 
        n_top_genes = total_hvgs, 
        flavor = 'seurat',
        batch_key = batch,
        inplace = False
    )

    union_hvg = hvgs.query('highly_variable').index.values
    ## adding hvgs in adata_dict: does not seem to be neceessary
    # for batch in adata_dict.keys():

    ## reset adata 
    adata.X = adata.layers['counts']
    del adata.var
    del adata.uns
    del adata.layers

    if verbose:
        print(f'\t - Number of total highly variable genes: {len(union_hvg)}')

    ## saving hvgs to disk 
    pd.DataFrame(union_hvg).to_csv(
        output_file, index = False, header = False
    )

    if verbose:
        print(f'\t - Saving highly variable genes: {output_file}')


    return union_hvg


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run data preprocessing')

    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    parser.add_argument('-b', '--batch', required = True, help = 'Batch variable')
    parser.add_argument('-v', '--hvgs', help = 'Number of highly variable genes', default = 2500)
    parser.add_argument('-n', '--normalize', action = 'store_true', help = 'Normalize and log-transform data')
    parser.add_argument('-s', '--scale', action = 'store_true', help = 'Scale the data per batch')

    args = parser.parse_args()
    file = args.input
    out = args.output
    batch = args.batch
    hvg = int(args.hvgs)
    normalize = args.normalize
    scale = args.scale

    print("Running preprocessing ...")
    print(f"\t - Input: {file}")
    print(f"\t - Output: {out}")
    print(f"\t - Batch variable: {batch}")
    print(f"\t - Number of highly variable genes: {hvg}")
    print(f"\t - Normalize: {normalize}")
    print(f"\t - Scale: {scale}")

    run_preprocessing(
        input = file, output = out, hvg = hvg, batch = batch, 
        normalize = normalize, scale = scale
    )
    
