#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { PREPROCESSING } from './modules/preprocessing'
include { SAVE_SEURAT } from './modules/save_seurat'
include { INTEGRATION_PYTHON } from './modules/integration_py'
include { INTEGRATION_R } from './modules/integration_r'
include { METRICS } from './modules/metrics'

// Parameter validation
if (!params.batch) {
    error "Please specify --batch parameter"
}
if (!params.label_key) {
    error "Please specify --label_key parameter"
}

// Main workflow
workflow {
    // Initialize channels
    if (params.run_preprocessing && params.input) {
        input_ch = Channel.fromPath(params.input, checkIfExists: true)
        
        PREPROCESSING(input_ch, params.batch, params.label_key, params.hvg, params.organism)
        preprocessed_h5ad = PREPROCESSING.out.h5ad
        
        SAVE_SEURAT(preprocessed_h5ad, params.batch, params.label_key)
        preprocessed_rds = SAVE_SEURAT.out.rds
        
    } else {
        // Use provided preprocessed files
        preprocessed_h5ad = Channel.fromPath(params.input_h5ad, checkIfExists: true)
        preprocessed_rds = Channel.fromPath(params.input_rds, checkIfExists: true)
    }
    
    // Python integration
    if (params.run_integration_py) {
        INTEGRATION_PYTHON(preprocessed_h5ad, params.batch, params.label_key, params.methods)
        integrated_py = INTEGRATION_PYTHON.out.integrated
    } else {
        integrated_py = Channel.empty()
    }
    
    // R integration
    if (params.run_integration_r) {
        INTEGRATION_R(preprocessed_rds, params.batch, params.label_key, params.methods)
        integrated_r = INTEGRATION_R.out.integrated
    } else {
        integrated_r = Channel.empty()
    }
    
    // Combine all integrated results
    all_integrated = integrated_py.mix(integrated_r)
    
    // Compute metrics
    if (params.run_metrics) {
        METRICS(all_integrated, preprocessed_h5ad, params.batch, params.label_key)
    }
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}