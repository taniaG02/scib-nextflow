// main.nf
nextflow.enable.dsl = 2

// Parámetros para controlar qué pasos ejecutar
params.run_step1 = true
params.run_step2 = true
params.run_step3 = true
params.run_step4 = true

// Parámetros de entrada
params.reads = "data/*_{1,2}.fastq.gz"
params.genome = "reference/genome.fa"
params.annotation = "reference/genes.gtf"

// Definición de procesos (manteniendo tu lógica original)
process PREPROCESSING {
    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("trimmed_*"), emit: out
    
    script:
    """
    echo "Preprocesando muestra: $sample_id"
    # Tu comando de preprocesamiento aquí
    touch trimmed_${sample_id}_R1.fq
    touch trimmed_${sample_id}_R2.fq
    """
}

process ASSEMBLY {
    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("assembly_*"), emit: out
    
    script:
    """
    echo "Ensamblando muestra: $sample_id"
    # Tu comando de ensamblaje aquí
    touch assembly_${sample_id}.fa
    """
}

process MAPPING {
    input:
        tuple val(sample_id), path(reads)
        path genome
    
    output:
        tuple val(sample_id), path("mapped_*"), emit: out
    
    script:
    """
    echo "Mapeando muestra: $sample_id"
    # Tu comando de mapeo aquí
    touch mapped_${sample_id}.bam
    """
}

process QUANTIFICATION {
    input:
        tuple val(sample_id), path(bam)
        path annotation
    
    output:
        tuple val(sample_id), path("quant_*"), emit: out
    
    script:
    """
    echo "Cuantificando muestra: $sample_id"
    # Tu comando de cuantificación aquí
    touch quant_${sample_id}.txt
    """
}

workflow {
    // Canal de entrada principal
    Channel
        .fromFilePairs(params.reads, size: 2)
        .set { reads_ch }
    
    // Canal para genoma y anotación
    Channel
        .fromPath(params.genome)
        .first()
        .set { genome_ch }
    
    Channel
        .fromPath(params.annotation)
        .first()
        .set { annotation_ch }
    
    // Variables para almacenar salidas
    out1 = Channel.empty()
    out2 = Channel.empty()
    out3 = Channel.empty()
    
    // PASO 1: Preprocesamiento
    if (params.run_step1) {
        PREPROCESSING(reads_ch)
        out1 = PREPROCESSING.out
    } else {
        // Si no se ejecuta el paso 1, usamos los datos crudos directamente
        out1 = reads_ch
    }
    
    // PASO 2: Ensamblaje
    if (params.run_step2) {
        if (params.run_step1) {
            ASSEMBLY(out1)
            out2 = ASSEMBLY.out
        } else {
            // Permitir ejecutar paso 2 sin paso 1 usando datos crudos
            ASSEMBLY(reads_ch)
            out2 = ASSEMBLY.out
        }
    } else {
        out2 = Channel.empty()
    }
    
    // PASO 3: Mapeo
    if (params.run_step3) {
        if (params.run_step1 || params.run_step2) {
            // Usar salida del paso 1 si está disponible, sino del paso 2
            input_mapeo = params.run_step1 ? out1 : out2
            MAPPING(input_mapeo, genome_ch)
            out3 = MAPPING.out
        } else {
            // Permitir ejecutar paso 3 con datos crudos si no hay pasos previos
            MAPPING(reads_ch, genome_ch)
            out3 = MAPPING.out
        }
    } else {
        out3 = Channel.empty()
    }
    
    // PASO 4: Cuantificación
    if (params.run_step4) {
        if (params.run_step3) {
            QUANTIFICATION(out3, annotation_ch)
        } else if (params.run_step1 || params.run_step2) {
            // Crear canal dummy para mapeo si no se ejecutó
            QUANTIFICATION(Channel.empty(), annotation_ch)
        } else {
            error "Paso 4 requiere salida de mapeo o datos preprocesados"
        }
    }
}


// Ejemplos de uso
// Ejecutar todo:
// nextflow run main.nf
//Ejecutar solo preprocesamieto:
// nextflow run main.nf --run_step1 true --run_step2 false --run_step3 false --run_step4 false