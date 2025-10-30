#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { PREPROCESSING } from '/home/tgonzalos/Documents/scib-nextflow/modules/preprocessing.nf'
include { SAVE_SEURAT } from '/home/tgonzalos/Documents/scib-nextflow/modules/save_seurat.nf'
include { INTEGRATION_PY } from '/home/tgonzalos/Documents/scib-nextflow/modules/integration_py.nf'
include { INTEGRATION_R } from '/home/tgonzalos/Documents/scib-nextflow/modules/integration_r.nf'
include { METRICS } from '/home/tgonzalos/Documents/scib-nextflow/modules/metrics.nf'

//include { PREPROCESSING } from "${projectDir}/modules/preprocessing.nf"
//include { SAVE_SEURAT }   from "${projectDir}/modules/save_seurat.nf"
//include { INTEGRATION_PY } from "${projectDir}/modules/integration_py.nf"
//include { INTEGRATION_R }  from "${projectDir}/modules/integration_r.nf"
//include { METRICS }        from "${projectDir}/modules/metrics.nf"


// Parameter validation
if (!params.batch) {
    error "Please specify --batch parameter"
}
if (!params.label_key) {
    error "Please specify --label_key parameter"
}

// Conditional execution parameters
params.run_preprocessing = true
params.run_save_seurat = true
params.run_integration_py = true
params.run_integration_r = true
params.run_metrics = true
params.methods = 'all'  // por defecto 'all'

// Available methods
def all_methods_py = ['scanorama', 'bbknn', 'scvi', 'combat']
def all_methods_r  = ['Seurat-CCA', 'Seurat-RPCA', 'harmony', 'liger', 'fastmnn']

// Convert parameter to list
def methods_to_run = params.methods == 'all' ? 
    all_methods_py + all_methods_r :
    params.methods.split(',').collect{ it.trim() }

// Filter by environment
def methods_py = methods_to_run.findAll { all_methods_py.contains(it) }
def methods_r  = methods_to_run.findAll { all_methods_r.contains(it) }

// Parameters for alternative entries (optional)
params.input = null
params.input_h5ad = null
params.input_rds = null
params.input_integrated = null 

// Main workflow
workflow {

    // ----------------------------
    // PRE-PROCESSING
    // ----------------------------
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

    // ----------------------------
    // PYTHON INTEGRATION
    // ----------------------------
    if (params.run_integration_py && methods_py) {
        methods_py_ch = Channel.fromList(methods_py)
        py_inputs = methods_py_ch.combine(preproc_h5ad)
        INTEGRATION_PY(py_inputs, params.batch, params.hvg ?: 2000)
        integrated = INTEGRATION_PY.out.integrated
    }

    // ----------------------------
    // R INTEGRATION
    // ----------------------------
    if (params.run_integration_r && methods_r) {
        methods_r_ch = Channel.fromList(methods_r)
        r_inputs = methods_r_ch.combine(preproc_rds)
        INTEGRATION_R(r_inputs, params.batch, params.hvg ?: 2000)
        integrated = integrated.mix(INTEGRATION_R.out.integrated)
    }

    // ----------------------------
    // Alternative entries
    // ----------------------------
    if (!params.run_integration_py && !params.run_integration_r && params.input_integrated) {
        def files_list = params.input_integrated.split(',').collect{ it.trim() }.collect { file(it) }
        integrated = Channel.fromList(files_list)
    }

    // ----------------------------
    // MÉTRICS
    // ----------------------------
    if (params.run_metrics) {

        println "[INFO] Iniciando cálculo de métricas (esperando integración completada)"

        // Mapea los archivos integrados a su método
        integrated_files_ch = integrated.map { f ->
            def name = f.getBaseName()
            def method = name.replaceFirst(/-integrated$/, '')
            tuple(method, f)
        }

        // Combina con el archivo preprocesado (no corregido)
        metrics_inputs_ch = integrated_files_ch.combine(preproc_h5ad)

        // Ejecuta métricas solo cuando los inputs estén disponibles
        METRICS(
            metrics_inputs_ch,
            params.batch,
            params.label_key,
            params.organism,
            params.hvg
        )
    }
}
