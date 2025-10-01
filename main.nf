include { PREPROCESSING } from '/home/tgonzalos/documents/scib-nextflow/modules/preprocessing.nf'
include { SAVE_SEURAT } from '/home/tgonzalos/documents/scib-nextflow/modules/save_seurat.nf'
include { INTEGRATION_PY } from '/home/tgonzalos/documents/scib-nextflow/modules/integration_py.nf'
include { INTEGRATION_R } from '/home/tgonzalos/documents/scib-nextflow/modules/integration_r.nf'
include { METRICS } from '/home/tgonzalos/documents/scib-nextflow/modules/metrics.nf'

// Parámetros de ejecución condicional
params.run_preprocessing = true
params.run_save_seurat = true
params.run_integration_py = true
params.run_integration_r = true
params.run_metrics = true

params.methods = 'all'  // por defecto 'all'

// Métodos disponibles
def all_methods_py = ['scanorama', 'bbknn', 'scvi', 'combat']
def all_methods_r  = ['Seurat-CCA', 'Seurat-RPCA', 'harmony', 'liger', 'fastmnn']

// Convertir parámetro en lista
def methods_to_run = params.methods == 'all' ? 
    all_methods_py + all_methods_r :
    params.methods.split(',').collect{ it.trim() }

// Filtrar por entorno
def methods_py = methods_to_run.findAll { all_methods_py.contains(it) }
def methods_r  = methods_to_run.findAll { all_methods_r.contains(it) }

// Parámetros para entradas alternativas (opcionales)
params.input = null
params.input_h5ad = null
params.input_rds = null
params.input_integrated = null 


workflow {

    // PREPROCESAMIENTO
    if (params.run_preprocessing) {
        if (!params.input) error "Se requiere --input cuando --run_preprocessing=true"
        input_ad = Channel.fromPath(params.input)
    } else {
        if (!params.input_h5ad) error "Se requiere --input_h5ad cuando --run_preprocessing=false"
        if (!params.input_rds) error "Se requiere --input_rds cuando --run_preprocessing=false"
    }

    preproc_h5ad = Channel.empty()
    preproc_rds  = Channel.empty()
    integrated   = Channel.empty()

    if (params.run_preprocessing) {
        PREPROCESSING(input_ad, params.batch, params.hvg)
        preproc_h5ad = PREPROCESSING.out.h5ad
    } else {
        preproc_h5ad = Channel.fromPath(params.input_h5ad)
    }

    if (params.run_preprocessing || !params.input_rds) {
        SAVE_SEURAT(preproc_h5ad, params.batch, params.hvg)
        preproc_rds = SAVE_SEURAT.out.rds
    } else {
        preproc_rds = Channel.fromPath(params.input_rds)
    }

    // INTEGRACIÓN PYTHON
    if (params.run_integration_py && methods_py) {
        methods_py_ch = Channel.fromList(methods_py)
        py_inputs = methods_py_ch.combine(preproc_h5ad)
        INTEGRATION_PY(py_inputs, params.batch, params.hvg ?: 2000)
        integrated = INTEGRATION_PY.out.integrated
    }

    // INTEGRACIÓN R
    if (params.run_integration_r && methods_r) {
        methods_r_ch = Channel.fromList(methods_r)
        r_inputs = methods_r_ch.combine(preproc_rds)
        INTEGRATION_R(r_inputs, params.batch, params.hvg ?: 2000)
        integrated = integrated.mix(INTEGRATION_R.out.integrated)
    }

    // Entradas alternativas
    if (!params.run_integration_py && !params.run_integration_r && params.input_integrated) {
        def files_list = params.input_integrated.split(',').collect{ it.trim() }.collect { file(it) }
        integrated = Channel.fromList(files_list)
    }

    integrated_py = Channel.empty()
    integrated_r  = Channel.empty()

    // MÉTRICAS
    if (params.run_metrics) {
        println "[INFO] Iniciando cálculo de métricas"

        // Buscar todos los ficheros integrados disponibles
        integrated_files_ch = Channel.fromPath('results/Integrated/*-integrated.h5ad', checkIfExists: true)
            .map { f ->
                def name = f.getBaseName()
                def method = name.replaceFirst(/-integrated$/, '')
                tuple(method, f)
            }

        // Combinar con el preprocesado (uncorrected_adata)
        metrics_inputs_ch = integrated_files_ch.combine(preproc_h5ad)

        // Ejecutar METRICS por cada (method, integrated_file, preproc_h5ad)
        METRICS(
            metrics_inputs_ch,
            params.batch,
            params.label_key,
            params.organism,
            params.hvg
        )
    }

}
