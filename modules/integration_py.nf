process INTEGRATION_PY {

    tag { method }
    publishDir "${params.output}/Integrated", mode: 'copy'

    conda = '/home/tgonzalos/miniconda3/envs/scib-cNMF-env'

    input:
    tuple val(method), path(h5ad_file)
    val batch
    val hvg

    output:
    path "*-integrated.h5ad", emit: integrated    

    script:
    """
    python ${workflow.projectDir}/bin/py_integration.py \
        -i "${h5ad_file}" \
        -o . \
        -b ${batch} \
        -m ${method} \
        -v ${hvg}
    """
}

