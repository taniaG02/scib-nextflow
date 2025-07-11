process INTEGRATION_R {
    tag { method }
    publishDir "${params.output}/Integrated", mode: 'copy'
    conda "conda/R-integration-env.yml"
    
    input:
    path seurat_rds
    val method
    val batch
    val hvg
    
    output:
    path "${method}-integrated.h5ad", emit: integrated
    
    script:
    """
    Rscript bin/R_integration.R \
        -i $seurat_rds \
        -o . \
        -b "$batch" \
        -m $method \
        -v $hvg \
        -s "${projectDir}/bin"
    """
}