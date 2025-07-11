#!/usr/bin/env python
# coding: utf-8

import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'

import scanpy as sc
import scib
import warnings
import pandas as pd

warnings.filterwarnings('ignore')


def save_seurat(output, batch, hvg):
    """
    params:
        input: path of the anndata object
        output: path of the preprocessed file to be written
        batch: variable in cells metadata to be used as batch
        hvg: number of highly variable genes to use
        scale: set to true to activate scaling
    """

    adata = sc.read(os.path.join(output, 'adata-preprocessed.h5ad'))
   
    print("Save as RDS")
    scib.preprocessing.save_seurat(
        adata = adata, 
        path = os.path.join(output, 'seurat-preprocessed.rds'), 
        batch = batch, hvgs = hvg
    )
 

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run data preprocessing')

    parser.add_argument('-o', '--output', required = True)
    parser.add_argument('-b', '--batch', required = True, help = 'Batch variable')
    parser.add_argument('-v', '--hvgs', help = 'Number of highly variable genes', default = 2000)

    args = parser.parse_args()
    out = args.output
    batch = args.batch
    hvg = int(args.hvgs)

    save_seurat(out, hvg, batch)
    
