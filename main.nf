include PREPROCESSING from './modules/preprocessing'
include INTEGRATION_PY from './modules/integration_py'
include INTEGRATION_R from './modules/integration_r'
include METRICS from './modules/metrics'

// Parámetros de ejecución condicional
params.run_preprocessing = true
params.run_integration_py = true
params.run_integration_r = true
params.run_metrics = true

params.methods_py = ['scanorama', 'bbknn', 'scvi', 'combat', 'NMFusion-CPMs', 'NMFusion-counts']
params.methods_r = ['Seurat-CCA', 'Seurat-RPCA', 'harmony', 'liger', 'fastmnn']

// Parámetros para entradas alternativas (opcionales, solo necesarios si se salta un paso)
params.input = null
params.input_h5ad = null
params.input_rds = null
params.input_integrated_py = null
params.input_integrated_r = null

workflow {
    // Validación: si se ejecuta preprocesamiento, se requiere params.input
    if (params.run_preprocessing) {
        if (!params.input) error "Se requiere --input cuando --run_preprocessing=true"
        input_ad = Channel.fromPath(params.input)
    } else {
        // Validación: si no se ejecuta preprocesamiento, se requieren los archivos preprocesados
        if (!params.input_h5ad) error "Se requiere --input_h5ad cuando --run_preprocessing=false"
        if (!params.input_rds) error "Se requiere --input_rds cuando --run_preprocessing=false"
    }

    // Inicialización explícita de canales
    preproc_h5ad = Channel.empty()
    preproc_rds  = Channel.empty()
    integrated_py = Channel.empty()
    integrated_r = Channel.empty()

    // 1. PREPROCESAMIENTO
    if (params.run_preprocessing) {
        Channel.fromPath(params.input).set { input_ad }
        PREPROCESSING(input_ad, params.batch, params.hvg)
        preproc_h5ad = PREPROCESSING.out.h5ad
        preproc_rds  = PREPROCESSING.out.rds
    } else {
        preproc_h5ad = Channel.fromPath(params.input_h5ad)
        preproc_rds  = Channel.fromPath(params.input_rds)
    }

    // 2. INTEGRACIÓN EN PYTHON
    if (params.run_integration_py) {
        methods_py = Channel.fromList(params.methods_py)
        INTEGRATION_PY(
            preproc_h5ad,
            methods_py,
            params.batch,
            params.hvg
        )
        integrated_py = INTEGRATION_PY.out.integrated
    } else if (params.input_integrated_py) {
        integrated_py = Channel.fromPath(params.input_integrated_py, checkIfExists: true)
    }

    // 3. INTEGRACIÓN EN R
    if (params.run_integration_r) {
        methods_r = Channel.fromList(params.methods_r)
        INTEGRATION_R(
            preproc_rds,
            methods_r,
            params.batch,
            params.hvg
        )
        integrated_r = INTEGRATION_R.out.integrated
    } else if (params.input_integrated_r) {
        integrated_r = Channel.fromPath(params.input_integrated_r, checkIfExists: true)
    }

    // 4. MÉTRICAS
    if (params.run_metrics) {
        all_integrations = integrated_py.mix(integrated_r)
        // Verificar que hay al menos un elemento en all_integrations
        if (all_integrations.count().toInteger() == 0) {
            log.warn "No hay integraciones para evaluar. El paso de métricas se omite."
        } else {
            METRICS(
                all_integrations,
                preproc_h5ad, // uncorrected
                // preproc_h5ad.first(), // Usar .first() para evitar duplicados
                params.batch,
                params.labelkey,
                params.organism,
                params.hvg
            )
        }
    }
}


// Ejemplo de uso
// nextflow run main.nf --run_preprocessing false \
//                     --input_h5ad path/to/preprocessed.h5ad \
//                     --input_rds path/to/preprocessed.rds \
//                     --run_integration_py true \
//                     --run_integration_r false
