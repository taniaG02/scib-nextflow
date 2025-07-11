process INTEGRATION_PY {
    tag { method }
    publishDir "${params.output}/Integrated", mode: 'copy'
    conda "conda/scib-cNMF-env.yml"
    
    input:
    path adata
    val method
    val batch
    val hvg
    
    output:
    path "*-integrated.h5ad", emit: integrated  # Captura múltiples salidas
    
    script:
    def nmfusion_path = "${projectDir}/lib/nmfusion"  # Cambiar según tu instalación
    """
    python bin/py_integration.py \
        -i $adata \
        -o . \
        -b "$batch" \
        -m $method \
        -v $hvg \
        --nmfusion_path "$nmfusion_path"
    """
}