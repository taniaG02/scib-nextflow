process INTEGRATION_R {
    tag { method }
    publishDir "${params.output}/Integrated", mode: 'copy'
    
    conda = '/home/tgonzalos/miniconda3/envs/r-integration-env'

    input:
	tuple val(method), path(seurat_rds)
    val batch
    val hvg
    
    output:
    path "${method}-integrated.h5ad", emit: integrated
    
    script:
    """
    Rscript ${workflow.projectDir}/bin/R_integration.R \
        -i $seurat_rds \
        -o . \
        -b "$batch" \
        -m $method \
        -v $hvg \
        -s "${projectDir}/bin"
    """
}
