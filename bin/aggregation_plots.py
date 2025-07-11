#!/usr/bin/env python
# coding: utf-8

import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'

import scanpy as sc
import anndata as ad
import glob
import warnings
import pandas as pd
import matplotlib.pyplot as plt


warnings.filterwarnings('ignore')


def agg_metrics(metrics_dir):
    '''
    Of note: this function only works on Module 3 output.
    '''
    
    csv_files = glob.glob(os.path.join(metrics_dir, '*_metrics.csv'))
    df = pd.concat([pd.read_csv(file) for file in csv_files], axis = 1)
    metrics = df.iloc[:, 0].values
    df_mod = df.loc[:, df.columns != 'Unnamed: 0']
    df_mod.index = metrics

    ## saving to disk
    df_mod.to_csv(os.join.file(metrics_dir, 'agg-metics.csv'))

    ## generating barplot
    df.plot(kind = 'bar', figsize = (20, 6))

    # Formatting the plot
    plt.xlabel('Metrics')
    plt.ylabel('Values')
    plt.title('Barplot of metrics across methods')
    plt.legend(title = 'Methods')
    plt.xticks(rotation = 45)  
    plt.savefig(os.join.path(metrics_dir, 'barplot-metrics.png'))

    return df_mod

def generate_umaps():



def save_seurat_rds():



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Module 4: Aggregation metrics and generation of UMAP plots'
    )
    
    parser.add_argument('-i', '--integrated', help = 'Directory with integrated adata objects', required = True)
    parser.add_argument('-m', '--metrics', help = 'Driectory of metrics. This location must contain the csv files output from Module 3: Computing metrics', required = True)
    parser.add_argument('-b', '--batch', required = True, help = 'Batch variable')
    parser.add_argument('-c', '--celltype', help = 'Cell type variable', default = None)
    parser.add_argument('-v', '--hvgs', help = 'Number of HVGs', default = 2500)

    args = parser.parse_args()
    integrated_dir = args.integrated
    metrics_dir = args.metrics
    batch = args.batch
    celltype = args.celltype
    hvgs = args.hvgs

    ## aggregate metrics into a single csv
    df = agg_metrics(metrics_dir = metrics_dir)
    

    
    
