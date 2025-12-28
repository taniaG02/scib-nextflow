process SAVE_SEURAT {
    
    label 'high_memory'
	publishDir "${params.output ?: 'results'}/preprocessing", mode: 'copy'

    conda = "${params.conda_base}/envs/r-integration-env"

    input:
        path adata_preprocessed
        val batch
        val hvg

    output:
        path "seurat-preprocessed.rds", emit: rds
        path "seurat-preprocessed.rds_hvg.RDS", emit: rds_hvg

    script:
    """
    python ${workflow.projectDir}/bin/save-seurat.py \
        -o . \
        -b "$batch" \
        -v $hvg
    """
}

