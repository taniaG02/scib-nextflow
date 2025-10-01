# ScIB Nextflow pipelines for benchmarking single-cell RNA-seq data integration methods

This repository contains a Nextflow pipeline for processing and integrating single-cell RNA-seq data, with extended capabilities to run individual steps and considerate dependencies.

Pipeline inspired in <https://github.com/theislab/scib-pipeline.git> from Luecken et al., 2022. For methods implemented in Python, it uses the `scib` Python module (<https://github.com/theislab/scib.git>). For methods implemented in R, check `src/R-integration-func.R`.


### Pipeline design

There are three main steps: 

1. Preprocessing: done in Python, the idea is to compute HVGs using the batch-aware `scib` implementation, normalize data using CPMs, and perform a batch-aware scaling. The last step is not actually necessary.  
2. Integration: this part is divided into two parts: 
    * Methods implemented in R. 
    * Methods implemented in Python. 
3. Computing metrics: take the results and compute metrics using `scib`. 

Regarding how the pipeline is organized: 

* `src`: contains the scripts calling methods from `scib` and the R functions implemented to run each method. 
* `main.sh`: this script will call the functions for each step. This is probably not very optimal and it would be better to implement a Nextflow workflow. Maybe in future releases. 
    **Parameters:**
    * input: path to an AnnData
    * output: directory used for output
    * batch: batch variable
    * hvg: 2000
    * methods: methods to run. It can be an array with the specific methods or just all. 
    In those cases where job arrays can be used, `main.sh` will call another sh script which will be in charge of calling the Python/R script.
* `integration-R.sh`: it calls the R script. 
* `integration-py.sh`: it calls the Python script. 
* `metrics.sh`: it can also be parallelized using jobarrays. 

### Installation 

I have created some conda environments that in principle should work out . The yml files to create them are in `conda-envs`. As some packages show conflicts in the versions that can be installed, the conda envs are two (for now): 

* `scib-cNMF-env`: it contains the python integration methods, `scib` and `NMFusion`. Not all the integration methods available in `scib` have been used. See below the comments on each of them: 
    * `bbknn`: OK
    * `combat`: uses the `scanpy` implementation: `scanpy.pp.combat`.
    * `mnnpy`: this is an implementation of the `scran` integration method. Does not work, `fastmnn` should be enough. 
    * `scanorama`: `pip install scanorama`: OK
    * `scanvi`: OK, is part of scvi.
    * `scgen`: OK
    * `scvi`: OK
    * `trvae`: seems to depend on python 3.6. They have an implementation in pytorch that shouldn't conflict with this env. 
    * `trvaep`: pytorch version. Does not work. 

* `R-integration-env`: it contains the dependencies for those methods implemented in R. In particular: 
    * `rliger`: OK
    * `Seurat` methods: CCA and RPCA: OK
    * `fastMNN`: OK
    * `harmony`: OK
    * `conos`: OK

* `scib-metrics-env`: it contains the dependencies for computing the metrics. 

### Making conda environments portable 

It worked once with Carlos, so let's see what happens. The idea is to change all the paths in each environment so that it can be used independently of the user. The workflow would be as follows: 

**Note:** It seems he created a complete miniconda, not just the env. Makes sese: there must be other files required in other locations. What I don't get is that I created a new env in my preexisting conda, bt anyway. 

1. Copy my miniconda changing its name. 
2. Remove all envs not related to scib
3. Create a txt with all files with the new path. 
4. Run the following line of code to substitute these paths with a new one (by CTorroja): 

```bash
cat filesWithPathsToChange.txt | xargs sed -i -e 's|/ data3/PipeLines_scripts/Rhapsody/miniconda3.v4.12|EL_NUEVO_PATH/miniconda3.v4.12|g'
```

This step must be done by Lucía. In principle, this should work out. In addition, check the `` file, which loads this new conda for its use. 

1. Checking files with my path

```bash
grep -r "/home/dmananesc" scib-metrics-env > scib-metrics-env.txt
grep -r "/home/dmananesc" scib-cNMF-env > scib-cNMF-env.txt
grep -r "/home/dmananesc" R-integration-env > R-integration-env-path.txt
```

Carlos sent me this code, maybe its better: 

```bash
grep -rP "dmananesc" scib-metrics-env | cut -f1 -d":" | sort | uniq > scib-metrics-env.txt
grep -rP "dmananesc" scib-cNMF-env | cut -f1 -d":" | sort | uniq > scib-cNMF-env.txt
grep -rP "dmananesc" R-integration-env | cut -f1 -d":" | sort | uniq > R-integration-env.txt
```

Then, to remove those binary files where a match has been found: 

```bash 
grep -v "Binary" scib-metrics-env.txt > scib-metrics-env-no-bin.txt
grep -v "Binary" scib-cNMF-env.txt > scib-cNMF-env-no-bin.txt
grep -v "Binary" R-integration-env.txt > R-integration-env-no-bin.txt
```

2. Modifying matches

```bash
cat filesWithPathsToChange.txt | xargs sed -i -e 's|/ data3/PipeLines_scripts/Rhapsody/miniconda3.v4.12|EL_NUEVO_PATH/miniconda3.v4.12|g'
```

## References

<table>
  <tr><td> Luecken, M.D., Büttner, M., Chaichoompu, K. et al. (2022). Benchmarking atlas-level data integration in single-cell genomics.
  <i>Nat Methods</i>
   <b>19</b> 41-50
  <a href='https://doi.org/10.1038/s41592-021-01336-8'>doi:10.1038/s41592-021-01336-8</a>
  </td></tr>
</table>


# Implementation scIB-Nextflow

Este repositorio contiene una implementación de **[scIB](https://scib.readthedocs.io/en/latest/)** en **Nextflow**, que permite ejecutar de manera reproducible y modular un flujo de análisis para la **integración de datos de scRNA-seq** y el cálculo de métricas de calidad.

El pipeline incluye las siguientes etapas:  

1. **Preprocesamiento** de los datos (`.h5ad`)  
2. **Conversión a formato Seurat** (`.rds`)  
3. **Integración con métodos en Python** (Scanorama, BBKNN, scVI, ComBat)  
4. **Integración con métodos en R** (Seurat-CCA, Seurat-RPCA, Harmony, LIGER, FastMNN)  
5. **Cálculo de métricas** para evaluar la calidad de la integración  


## Requisitos

- [Nextflow](https://www.nextflow.io/) ≥ 22.10  
- [Conda](https://docs.conda.io/en/latest/) (o Mamba)  
- Python ≥ 3.8 y R ≥ 4.0  

Los módulos usan distintos entornos de conda:  
- `scib-cNMF-env` → preprocesamiento e integración en Python  
- `r-integration-env` → integración en R y conversión a Seurat  
- `scib-metrics-env` → cálculo de métricas  

## Estructura
scib-nextflow/
├── main.nf
├── nextflow.config
├── modules/
│   ├── preprocessing.nf
│   ├── save_seurat.nf
│   ├── integration_py.nf
│   ├── integration_r.nf
│   └── metrics.nf
├── bin/                  # Scripts auxiliares en Python/R
│   ├── preprocessing.py
│   ├── save-seurat.py
│   ├── py_integration.py
│   ├── R_integration.R
│   └── metrics.py
└── results/              # Resultados organizados por etapa



## Ejecución

### 1. Preprocesamiento completo + integración + métricas

```bash
nextflow run main.nf \
  --input data/mydata.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --organism human \
  --output results
```

### 2. Usar datos ya preprocesadosas

nextflow run main.nf \
  --input_h5ad results/preprocessing/adata-preprocessed.h5ad \
  --input_rds results/preprocessing/seurat-preprocessed.rds \
  --batch batch_column \
  --label_key cell_type \
  --organism mouse

### 3. Ejecutar solo métricas con integraciones existentes

nextflow run main.nf \
  --input_integrated results/Integrated/scanorama-integrated.h5ad,results/Integrated/harmony-integrated.h5ad \
  --input_h5ad results/preprocessing/adata-preprocessed.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --organism human \
  --run_preprocessing false \
  --run_integration_py false \
  --run_integration_r false \
  --run_metrics true


### Parámetros principales
| Parámetro              | Descripción                                                           |
| ---------------------- | --------------------------------------------------------------------- |
| `--input`              | Archivo `.h5ad` de entrada (solo si se ejecuta el preprocesamiento).  |
| `--input_h5ad`         | Archivo `.h5ad` ya preprocesado.                                      |
| `--input_rds`          | Archivo `.rds` de Seurat ya preprocesado.                             |
| `--input_integrated`   | Lista de archivos `.h5ad` integrados (para solo métricas).            |
| `--batch`              | Columna con la información de lote.                                   |
| `--label_key`          | Columna con las etiquetas biológicas (ej. tipo celular).              |
| `--organism`           | Organismo (`human` o `mouse`).                                        |
| `--hvg`                | Número de genes variables a usar (default: 2000).                     |
| `--methods`            | Métodos a ejecutar, separados por coma (`scanorama,harmony`) o `all`. |
| `--output`             | Directorio de salida (default: `results`).                            |
| `--run_preprocessing`  | `true/false` para ejecutar la etapa de preprocesamiento.              |
| `--run_integration_py` | Ejecutar integración en Python.                                       |
| `--run_integration_r`  | Ejecutar integración en R.                                            |
| `--run_metrics`        | Ejecutar cálculo de métricas.                                         |


## Resultados

Los resultados se organizan automáticamente en subdirectorios:

* results/preprocessing/ → archivos preprocesados (.h5ad, .rds)
* results/Integrated/ → integraciones por método (*-integrated.h5ad)
* results/Metrics/ → métricas por método (*_metrics.csv)

## Funcionamiento

- Preprocessing: toma los archivos .h5ad originales y los normaliza, filtrando genes y células según los parámetros.
- Save Seurat: convierte los datos preprocesados a .rds para los métodos de integración en R.
- Integration Python: ejecuta métodos como Scanorama, BBKNN, scVI o ComBat.
- Integration R: ejecuta métodos como Seurat-CCA, Seurat-RPCA, Harmony, LIGER o FastMNN.
- Metrics Calculation: compara los datos integrados con los originales para evaluar la calidad de la integración.

## Ejemplo de uso

nextflow run main.nf \
  --input data/pbmc.h5ad \
  --batch sample \
  --label_key cell_type \
  --organism human \
  --methods scanorama,harmony
