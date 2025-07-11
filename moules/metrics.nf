process METRICS {
    tag { method }
    publishDir "${params.output}/Metrics", mode: 'copy'
    conda "conda/scib-metrics-env.yml"
    
    input:
    tuple val(method), path(integrated_adata)
    path uncorrected_adata
    val batch_key
    val label_key
    val organism
    val hvg
    
    output:
    path "${method}_metrics.csv", emit: metrics
    
    script:
    // Determinar tipo automáticamente (como en metrics.py)
    def type_arg = "full"
    if (method in ['harmony', 'conos', 'liger', 'fastmnn', 'NMFusion-CPMs-NMF-space', 'NMFusion-counts-NMF-space']) {
        type_arg = "embed"
    } else if (method == 'bbknn') {
        type_arg = "knn"
    }
    
    """
    python bin/metrics.py \
        -u $uncorrected_adata \
        -i $integrated_adata \
        -o ${method}_metrics.csv \
        -m $method \
        -b "$batch_key" \
        -l "$label_key" \
        --organism $organism \
        --type $type_arg \
        --hvgs $hvg
    """
}
