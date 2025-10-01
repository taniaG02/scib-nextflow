process METRICS {
    
    tag { method }
    publishDir "${params.output ?: 'results'}/Metrics", mode: 'copy'

    conda = '/home/tgonzalos/miniconda3/envs/scib-metrics-env'

    input:
    tuple val(method), path(integrated_adata), path(uncorrected_adata)
    val batch_key
    val label_key
    val organism
    val hvg


    output:
    path "${method}_metrics.csv", emit: metrics

    script:
    """
    export R_HOME=${HOME}/miniconda3/envs/scib-metrics-env/lib/R
    export PATH=\$R_HOME/bin:\$PATH
    export LD_LIBRARY_PATH=\$R_HOME/lib/R/lib:\${LD_LIBRARY_PATH:-}

    # Determinar tipo automáticamente
    type_arg="full"
    if [[ "${method}" == "harmony" || "${method}" == "conos" || "${method}" == "liger" || "${method}" == "fastmnn" || "${method}" == "NMFusion-CPMs-NMF-space" || "${method}" == "NMFusion-counts-NMF-space" ]]; then
        type_arg="embed"
    elif [[ "${method}" == "bbknn" ]]; then
        type_arg="knn"
    fi

    python ${workflow.projectDir}/bin/metrics.py \
        -u ${uncorrected_adata} \
        -i ${integrated_adata} \
        -o ${method}_metrics.csv \
        -m ${method} \
        -b "${batch_key}" \
        -l "${label_key}" \
        --organism ${organism} \
        --type \$type_arg \
        --hvgs ${hvg}
    """
}



