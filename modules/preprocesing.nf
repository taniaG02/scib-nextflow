process PREPROCESSING {
    label 'high_memory'
    conda "conda/scib-cNMF-env.yml"
    
    input:
    path adata
    val batch
    val hvg
    
    output:
    path "adata-preprocessed.h5ad", emit: h5ad
    path "seurat_object.rds", emit: rds  # Asumido por R_integration.R
    
    script:
    """
    python bin/preprocessing.py \
        -i $adata \
        -o . \
        -b "$batch" \
        -v $hvg \
        --normalize
    
    # Script ficticio (adaptar a tu implementación real)
    python bin/save-seurat.py -i adata-preprocessed.h5ad -o seurat_object.rds
    """
}
