include PREPROCESSING from './modules/preprocessing'
include INTEGRATION_PY from './modules/integration_py'
include INTEGRATION_R from './modules/integration_r'
include METRICS from './modules/metrics'

params.methods_py = ['scanorama', 'bbknn', 'scvi', 'combat', 'NMFusion-CPMs', 'NMFusion-counts']
params.methods_r = ['Seurat-CCA', 'Seurat-RPCA', 'harmony', 'liger', 'fastmnn']

workflow {
    // Entrada de datos
    input_ad = Channel.fromPath(params.input)
    
    // 1. Preprocesamiento
    PREPROCESSING(input_ad, params.batch, params.hvg)
    
    // 2. Integración Python (paralela por método)
    methods_py = Channel.fromList(params.methods_py)
    INTEGRATION_PY(
        PREPROCESSING.out.h5ad,
        methods_py,
        params.batch,
        params.hvg
    )
    
    // 3. Integración R (paralela por método)
    methods_r = Channel.fromList(params.methods_r)
    INTEGRATION_R(
        PREPROCESSING.out.rds,
        methods_r,
        params.batch,
        params.hvg
    )
    
    // 4. Combinar resultados para métricas
    all_integrations = INTEGRATION_PY.out.integrated.mix(INTEGRATION_R.out.integrated)
    
    // 5. Cálculo de métricas (paralelo por método)
    METRICS(
        all_integrations,
        PREPROCESSING.out.h5ad,  // uncorrected
        params.batch,
        params.labelkey,
        params.organism,
        params.hvg
    )
}