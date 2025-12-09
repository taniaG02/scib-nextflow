// Process to combine individual metric CSV files into a single comparison table

process COMBINE_METRICS {
    tag "Combining metrics"
    publishDir "${params.outdir}/Metrics", mode: 'copy'
    
    conda params.conda_metrics_env
    
    input:
    path metric_files  // Collects all *_metrics.csv files
    
    output:
    path "combined_metrics.csv", emit: combined_metrics
    
    script:
    """
    combine_metrics.py \\
        --input-dir . \\
        --output combined_metrics.csv \\
        --pattern '*_metrics.csv'
    """
}

process PLOT_METRICS {
    tag "Generating metric plots"
    publishDir "${params.outdir}/Plots", mode: 'copy'
    
    conda params.conda_metrics_env
    
    input:
    path combined_metrics
    path reference_metrics  // Optional reference for comparison
    
    output:
    path "*.png", emit: plots
    
    script:
    def reference_arg = reference_metrics.name != 'NO_FILE' ? "--reference ${reference_metrics}" : ""
    """
    plot_metrics.py \\
        --input ${combined_metrics} \\
        ${reference_arg} \\
        --output-dir . \\
        --format png \\
        --dpi 300
    """
}
