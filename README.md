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


## Recently implemented changes

### 1. Individual step execution: 

Each module can be run independently using control parameters:

``` bash
# Run only preprocessing
nextflow run main.nf --run_preprocessing true --run_integration_py false --run_integration_r false --run_metrics false

# Run only Python integration
nextflow run main.nf --run_preprocessing false --run_integration_py true --run_integration_r false --run_metrics false

# Run only metrics with pre-existing data
nextflow run main.nf --run_preprocessing false --run_integration_py false --run_integration_r false --run_metrics true \
  --input_integrated_py "results/python_integrations/*.h5ad" \
  --input_integrated_r "results/r_integrations/*.rds"

```

### 1. Functional Improvements
#### Modular System
- Toggle steps with run_* parameters
- Alternative inputs for each module
- Automatic dependency validation

#### Robust Data Handling
- Load intermediate results
- Validate input files
- Support for h5ad (Python) and rds (R) formats

#### Enhanced Configuration
- Method-specific parameters
- Customizable analysis options

## Installation
```bash
# 1. Clone repository
git clone https://github.com/taniaG02/scib-nextflow.git
cd scib-nextflow

# 2. Install Nextflow
curl -s https://get.nextflow.io | bash

# 3. Create and activate Conda environment
conda env create -f environment.yml
conda activate scib-nextflow

```

## Basic usage
```bash
nextflow run main.nf \
  --input "data/*.h5ad" \
  --batch batch_key \
  --labelkey cell_type \
  --organism human
```
## Parameters

| Parameter              | Description                          | Default Value  |
|------------------------|--------------------------------------|---------------------|
| `--input`              | Path to input files       | `data/*.h5ad`       |
| `--batch`              | Batch/metadata column             | Required           |
| `--labelkey`           | Cell type column              | Required           |
| `--organism`           | Organism (human/mouse)             | `human`             |
| `--hvg`                | Numbr of highly variable genes  | `5000`              |
| `--run_preprocessing`  | Run preprocessgin step            | `true`              |
| `--run_integration_py` | Run Python integration       | `true`              |
| `--run_integration_r`  | Run R integration            | `true`              |
| `--run_metrics`        | Compute evaluation metrics                    | `true`              |


## Pipeline structure:



## Supported integration methods

Python: 
``` bash
--methods_py ['scanorama', 'bbknn', 'scvi', 'combat', 'NMFusion-CPMs', 'NMFusion-counts']
```

R:
``` bash
--methods_r ['Seurat-CCA', 'Seurat-RPCA', 'harmony', 'liger', 'fastmnn']
```

## Advanced execution
Use preprocessed data

``` bash
nextflow run main.nf \
  --run_preprocessing false \
  --input_h5ad "preprocessed/adata.h5ad" \
  --input_rds "preprocessed/data.rds" \
  --run_integration_py true \
  --run_metrics false
```

Combine specific integrations

``` bash
nextflow run main.nf \
  --run_preprocessing false \
  --run_integration_py false \
  --run_integration_r false \
  --run_metrics true \
  --input_integrated_py "integrations/python/scanorama.h5ad" \
  --input_integrated_r "integrations/r/harmony.rds"
```

Expected output structure

``` text
results/
├── preprocessing/
│   ├── adata.h5ad
│   └── data.rds
├── python_integrations/
│   ├── scanorama.h5ad
│   ├── bbknn.h5ad
│   └── ...
├── r_integrations/
│   ├── Seurat-CCA.rds
│   ├── harmony.rds
│   └── ...
└── metrics/
    ├── combined_metrics.csv
    ├── batch_correction/
    └── biological_conservation/
```
