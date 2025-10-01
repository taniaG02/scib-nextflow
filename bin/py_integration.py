#!/usr/bin/env python
# coding: utf-8

import argparse #Añadido
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'

import scanpy as sc
import scib
import warnings

warnings.filterwarnings('ignore')

def run_integration_py(
        input, output, batch, method, method_fun, 
        hvgs = 2500, celltype = None
):
    """
    params:
        input: path to the input file, adata already processed
        output: path to the output file, adata already with integration results
        batch: name of `adata.obs` column of the batch
        method: name of method
        hvgs: number of highly variable genes to be used
        celltype name of `adata.obs` column of cell types
    """

    adata = sc.read(input)

    if method == 'NMFusion-CPMs': 
        integrated = run_NMFusion_CPMs(input, output, batch, hvgs)
        ## out of the function to avoid compatibilities: format adata for metrics: NMF
        adata_NMF = integrated
        adata_NMF.obsm['X_emb'] = adata_NMF.obsm['NMF_k60_W_corr'].to_numpy()
        adata_NMF.write_h5ad(os.path.join(output, f'{method}-NMF-space-integrated.h5ad'))

    elif method == 'NMFusion-counts':
        integrated = run_NMFusion_counts(input, output, batch, hvgs)
        ## out of the function to avoid compatibilities: format adata for metrics: NMF
        adata_NMF = integrated
        adata_NMF.obsm['X_emb'] = adata_NMF.obsm['NMF_k60_W_corr'].to_numpy()
        adata_NMF.write_h5ad(os.path.join(output, f'{method}-NMF-space-integrated.h5ad'))
    else:
        if celltype is not None:
            integrated = method_fun(adata, batch, celltype)
        else:
            integrated = method_fun(adata, batch)

    integrated.write_h5ad(os.path.join(output, f'{method}-integrated.h5ad'))
    # sc.write(os.path.join(output, f'{method}-integrated.h5ad'), integrated)

def run_NMFusion_CPMs(input, output, batch, hvgs):
    """
    Function implementing the workflow to run NMFusion + additional steps for downstream analyses.
    """
    import sys
    # hardcoded, to be changed in future releases
    sys.path.append('/data_lab_DSM/Projects_3/Software/nmfusion')
    sys.path.append('/data_lab_DSM/Projects_3/Software/cnmf_new')
    try:
        from cnmf_new import cNMF_new
        import correction
        import factorize_data
        import utils
        
    except ImportError as e:
        print(f"Failed to import required modules")
        sys.exit(1)

    adata = sc.read(input)
    output_dir = os.path.join(output, 'NMFusion-CPMs-tmp')
    
    print(">>> Factorizing data")
    adata = factorize_data.decompose_adata(
        adata = adata, 
        batch = batch, 
        output_dir = output_dir,
        filter_perc = 0,
        k = 60, 
        total_hvgs = hvgs,
        keep_raw = True,
        remove_tmp_files = False,
        verbose = True,
        corr = 'cosine',
        warning_ignore = True
    )

    print(">>> Finding bio factors")
    adata = factorize_data.find_bio_factors(
        adata, 
        cutoff_agg = 0.6, 
        perc_trans = 0.6,
        perc_common = False,
        corr = 'spearman',
        n_clusters = 60
    )

    # Process bio factors with scanpy
    print(">>> Processing bio factors with scanpy")
    adata_NMF_bio_pca = adata.copy()
    
    sc.pp.normalize_total(adata_NMF_bio_pca, target_sum=1e4)
    sc.pp.log1p(adata_NMF_bio_pca)
    sc.pp.scale(adata_NMF_bio_pca, max_value=10)
    
    sc.tl.pca(adata_NMF_bio_pca, n_comps = 30)
    sc.pp.neighbors(adata_NMF_bio_pca, n_neighbors = 15)
    sc.tl.umap(adata_NMF_bio_pca)
    
    # Save bio factors result
    bio_factors_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_NMFusion_CPMs_bio_pca.h5ad')
    )
    print(f"Saving bio factors to: {bio_factors_output}")
    adata_NMF_bio_pca.write_h5ad(bio_factors_output)

    print(">>> Calculating correction factors")
    print('\tw/o H refit [just another varm is added]: ')
    adata = correction.corr_factors(
        adata = adata, 
        batch = batch, 
        n_iter = 500, 
        method = "nnls-02",
        transform = True, 
        stand = 'iter',
        copy = True, 
        verbose = True,
        tolerance = 0.00001,
        refit_H = True
    )

    print(">>> Running scanpy workflow")
    adata_lm = utils.workflow_scanpy_01(adata, k = 60)
    adata_lm_H = utils.workflow_scanpy_01(adata, k = 60, use_H_refit = True)
    
    # Save final results
    lm_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_lm_NMFusion_CPMs.h5ad')
    )
    lm_H_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_lm_H_refit_NMFusion_CPMs.h5ad')
    )
    final_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_NMFusion_CPMs.h5ad')
    )
    
    print(f"Saving final outputs to: {lm_output}, {lm_H_output} and {final_output}")
    adata_lm.write_h5ad(lm_output)
    adata_lm_H.write_h5ad(lm_H_output)
    adata.write_h5ad(final_output)

    adata = correction.get_corr_counts(adata, k = 60)

    ## running PCA on corrected counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps = 30)

    # format adata for metrics: PCA
    adata.obsm['X_emb'] = adata.obsm['X_pca']

    print('Workaround: to allow H_refit in the next module, the object is saved with the correct format here')
    adata_H = correction.get_corr_counts(adata, k = 60, use_H_refit = True)
    ## running PCA on corrected counts
    sc.pp.normalize_total(adata_H, target_sum=1e4)
    sc.pp.log1p(adata_H)
    sc.pp.scale(adata_H, max_value=10)
    sc.tl.pca(adata_H, n_comps = 30)
    # format adata for metrics: PCA
    adata_H.obsm['X_emb'] = adata_H.obsm['X_pca']
    adata_H.write_h5ad(os.path.join(output, 'NMFusion-CPMs-H-refit-integrated.h5ad'))

    return adata


def run_NMFusion_counts(input, output, batch, hvgs):
    """
    Function implementing the workflow to run NMFusion + additional steps for downstream analyses.
    """
    import sys
    # hardcoded, to be changed in future releases
    sys.path.append('/data_lab_DSM/Projects_3/Software/nmfusion')
    sys.path.append('/data_lab_DSM/Projects_3/Software/cnmf_new')
    try:
        from cnmf_new import cNMF_new
        import correction
        import factorize_data
        import utils
        
    except ImportError as e:
        print(f"Failed to import required modules")
        sys.exit(1)

    adata = sc.read(input)

    ## substitution of CPM by counts
    adata.X = adata.layers['counts']
    output_dir = os.path.join(output, 'NMFusion-counts-tmp')
    
    print(">>> Factorizing data")
    adata = factorize_data.decompose_adata(
        adata = adata, 
        batch = batch, 
        output_dir = output_dir,
        filter_perc = 0,
        k = 60, 
        total_hvgs = hvgs,
        keep_raw = True,
        remove_tmp_files = False,
        verbose = True,
        corr = 'cosine',
        warning_ignore = True
    )

    print(">>> Finding bio factors")
    adata = factorize_data.find_bio_factors(
        adata, 
        cutoff_agg = 0.6, 
        perc_trans = 0.6,
        perc_common = False,
        corr = 'spearman',
        n_clusters = 60
    )

    # Process bio factors with scanpy
    print(">>> Processing bio factors with scanpy")
    adata_NMF_bio_pca = adata.copy()
    
    sc.pp.normalize_total(adata_NMF_bio_pca, target_sum=1e4)
    sc.pp.log1p(adata_NMF_bio_pca)
    sc.pp.scale(adata_NMF_bio_pca, max_value=10)
    
    sc.tl.pca(adata_NMF_bio_pca, n_comps = 30)
    sc.pp.neighbors(adata_NMF_bio_pca, n_neighbors = 15)
    sc.tl.umap(adata_NMF_bio_pca)
    
    # Save bio factors result
    bio_factors_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_NMFusion_counts_bio_pca.h5ad')
    )
    print(f"Saving bio factors to: {bio_factors_output}")
    adata_NMF_bio_pca.write_h5ad(bio_factors_output)

    print(">>> Calculating correction factors")
    adata = correction.corr_factors(
        adata = adata, 
        batch = batch, 
        n_iter = 500, 
        method = "nnls-02",
        transform = True, 
        stand = 'iter',
        copy = True, 
        verbose = True,
        tolerance = 0.00001,
        refit_H = True
    )

    print(">>> Running scanpy workflow")
    adata_lm = utils.workflow_scanpy_01(adata, k = 60)
    adata_lm_H = utils.workflow_scanpy_01(adata, k = 60, use_H_refit = True)
    
    # Save final results
    lm_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_lm_NMFusion_counts.h5ad')
    )
    lm_H_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_lm_H_refit_NMFusion_counts.h5ad')
    )
    final_output = os.path.join(
        output, os.path.basename(input).replace('.h5ad', '_NMFusion_counts.h5ad')
    )
    
    print(f"Saving final outputs to: {lm_output} and {final_output}")
    adata_lm.write_h5ad(lm_output)
    adata_lm_H.write_h5ad(lm_H_output)
    adata.write_h5ad(final_output)

    adata = correction.get_corr_counts(adata, k = 60)

    ## running PCA on corrected counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps = 30)

    # format adata for metrics
    adata.obsm['X_emb'] = adata.obsm['X_pca']

    print('Workaround: to allow H_refit in the next module, the object is saved with the correct format here')
    adata_H = correction.get_corr_counts(adata, k = 60, use_H_refit = True)
    ## running PCA on corrected counts
    sc.pp.normalize_total(adata_H, target_sum=1e4)
    sc.pp.log1p(adata_H)
    sc.pp.scale(adata_H, max_value=10)
    sc.tl.pca(adata_H, n_comps = 30)
    # format adata for metrics: PCA
    adata_H.obsm['X_emb'] = adata_H.obsm['X_pca']
    adata_H.write_h5ad(os.path.join(output, 'NMFusion-counts-H-refit-integrated.h5ad'))

    return adata



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Run the integration methods implemented in python using scib'
    )
    
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    parser.add_argument('-b', '--batch', required = True, help = 'Batch variable')
    parser.add_argument('-m', '--method', required = True)
    parser.add_argument("-c", '--celltype', help='Cell type variable', default = None)
    parser.add_argument("-v", '--hvgs', help='Number of HVGs', default = 2500)

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    batch = args.batch
    celltype = args.celltype
    hvgs = args.hvgs
    method = args.method
    methods = {
        'scanorama': scib.integration.scanorama,
        'trvae': scib.integration.trvae,
        'trvaep': scib.integration.trvaep,
        'scgen': scib.integration.scgen,
        'mnn': scib.integration.mnn,
        'bbknn': scib.integration.bbknn,
        'scvi': scib.integration.scvi,
        'scanvi': scib.integration.scanvi,
        'combat': scib.integration.combat,
        'saucie': scib.integration.saucie,
        'desc': scib.integration.desc,
        #'NMFusion-CPMs': run_NMFusion_CPMs,
        #'NMFusion-counts': run_NMFusion_counts
    }

    if method not in methods.keys():
        raise ValueError(f'Method "{method}" does not exist. Please use one of '
                         f'the following:\n{list(methods.keys())}')
    
    run = methods[method]

    run_integration_py(
        input = input_file, 
        output = output_file, 
        method = method,
        method_fun = run, 
        batch = batch, 
        hvgs = hvgs,
        celltype = celltype
    )
    
    # #Añadido
    # def run_NMFusion_CPMs(input, output, batch, hvgs, nmfusion_path):
    #     sys.path.append(nmfusion_path)  # Usar parámetro
