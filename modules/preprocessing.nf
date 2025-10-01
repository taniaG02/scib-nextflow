process PREPROCESSING {
    
    label 'high_memory'
	publishDir "/home/tgonzalos/documents/TGS25/PE_TGS/results/preprocessing", mode: 'copy'

    conda = '/home/tgonzalos/miniconda3/envs/scib-cNMF-env'

    input:
        path adata
        val batch
        val hvg
    
    output:
        path "adata-preprocessed.h5ad", emit: h5ad
    
    script:
    """
    python ${workflow.projectDir}/bin/preprocessing.py \
        -i $adata \
        -o . \
        -b "$batch" \
        -v $hvg \
        --normalize
    """
}

