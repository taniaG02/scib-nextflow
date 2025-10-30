# scib-nextflow

A reproducible Nextflow pipeline for single-cell RNA-seq data integration and benchmarking, implementing the scIB (single-cell Integration Benchmarking) framework.

## Overview

This pipeline provides a modular, reproducible implementation of [scIB](https://scib.readthedocs.io/en/latest/) in Nextflow, enabling systematic benchmarking of single-cell integration methods across different datasets and batch correction strategies.

**Pipeline inspired by:** [theislab/scib-pipeline](https://github.com/theislab/scib-pipeline.git) (Luecken et al., 2022)

### Key Features

- **Modular architecture**: Run complete workflow or individual steps independently
- **Multiple integration methods**: Python-based (Scanorama, BBKNN, scVI, ComBat) and R-based (Seurat-CCA, Seurat-RPCA, Harmony, LIGER, FastMNN)
- **Comprehensive metrics**: Automated calculation of batch correction and biological conservation metrics
- **Environment isolation**: Separate conda environments for Python, R, and metrics computation
- **Flexible execution**: Toggle steps on/off, skip preprocessing, or compute only metrics
- **HPC-ready**: Configured for execution on local machines, HPC clusters, or cloud environments

## Pipeline Stages

The pipeline consists of four main stages:

1. **Preprocessing** (Python/Scanpy)
   - Batch-aware highly variable gene (HVG) selection using scIB implementation
   - CPM normalization
   - Batch-aware scaling (optional)

2. **Seurat Conversion** (Python)
   - Convert preprocessed AnnData (.h5ad) to Seurat format (.rds) for R-based methods

3. **Integration**
   - **Python methods**: Scanorama, BBKNN, scVI, scANVI, scGen, Combat, NMFusion
   - **R methods**: Seurat-CCA, Seurat-RPCA, Harmony, LIGER, FastMNN, Conos

4. **Metrics Computation** (Python/scIB)
   - Batch correction metrics: kBET, LISI, ASW (batch)
   - Biological conservation metrics: ARI, NMI, ASW (label), silhouette score
   - Combined scores and visualizations

## Installation

### Requirements

- Nextflow >= 21.04.0
- Conda/Mamba or Docker/Singularity
- For HPC: SLURM, SGE, or PBS job scheduler

### Quick Start

```bash
# Clone the repository
git clone https://github.com/taniaG02/scib-nextflow.git
cd scib-nextflow

# Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash

# Test with example data
nextflow run main.nf --input data/test.h5ad --batch batch --label_key celltype --organism human
```

## Usage

### Basic Usage

Run the complete pipeline:

```bash
nextflow run main.nf \
  --input data/mydata.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --organism human \
  --output results
```

### Advanced Usage

#### Run only preprocessing and Python integration:

```bash
nextflow run main.nf \
  --input data/mydata.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --organism human \
  --run_integration_r false \
  --run_metrics false
```

#### Skip preprocessing (use pre-processed data):

```bash
nextflow run main.nf \
  --input_h5ad results/preprocessing/adata-preprocessed.h5ad \
  --input_rds results/preprocessing/seurat-preprocessed.rds \
  --batch batch_column \
  --label_key cell_type \
  --run_preprocessing false
```

#### Compute metrics only (from existing integrated data):

```bash
nextflow run main.nf \
  --input_integrated "results/Integrated/scanorama-integrated.h5ad,results/Integrated/harmony-integrated.h5ad" \
  --input_h5ad results/preprocessing/adata-preprocessed.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --run_preprocessing false \
  --run_integration_py false \
  --run_integration_r false \
  --run_metrics true
```

#### Run specific methods only:

```bash
nextflow run main.nf \
  --input data/mydata.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --methods scanorama,harmony,scvi
```

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Input AnnData file (.h5ad) for preprocessing | `data/mydata.h5ad` |
| `--batch` | Column name containing batch information | `batch` |
| `--label_key` | Column name with biological labels (e.g., cell type) | `cell_type` |
| `--organism` | Organism for the dataset | `human` or `mouse` |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_h5ad` | Pre-processed AnnData file (skips preprocessing) | `null` |
| `--input_rds` | Pre-processed Seurat file (skips conversion) | `null` |
| `--input_integrated` | Comma-separated list of integrated .h5ad files | `null` |
| `--hvg` | Number of highly variable genes to select | `2000` |
| `--methods` | Integration methods to run (`all` or comma-separated) | `all` |
| `--output` | Output directory for results | `results` |

### Execution Control Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--run_preprocessing` | Execute preprocessing step | `true` |
| `--run_integration_py` | Execute Python integration methods | `true` |
| `--run_integration_r` | Execute R integration methods | `true` |
| `--run_metrics` | Execute metrics calculation | `true` |

## Output Structure

Results are automatically organized into subdirectories:

```
results/
├── preprocessing/
│   ├── adata-preprocessed.h5ad    # Preprocessed AnnData
│   └── seurat-preprocessed.rds    # Preprocessed Seurat object
├── Integrated/
│   ├── scanorama-integrated.h5ad  # Scanorama integration
│   ├── harmony-integrated.h5ad    # Harmony integration
│   ├── scvi-integrated.h5ad       # scVI integration
│   └── ...                        # Other methods
└── Metrics/
    ├── scanorama_metrics.csv      # Metrics for Scanorama
    ├── harmony_metrics.csv        # Metrics for Harmony
    └── ...                        # Metrics for other methods
```

## Repository Structure

```
scib-nextflow/
├── main.nf                    # Main workflow entry point
├── nextflow.config            # Pipeline configuration
├── modules/                   # Nextflow process modules
│   ├── preprocessing.nf       # Preprocessing module
│   ├── save_seurat.nf         # AnnData to Seurat conversion
│   ├── integration_py.nf      # Python integration methods
│   ├── integration_r.nf       # R integration methods
│   └── metrics.nf             # Metrics computation
├── bin/                       # Executable scripts
│   ├── preprocessing.py       # Preprocessing script
│   ├── save-seurat.py         # Conversion script
│   ├── py_integration.py      # Python integration wrapper
│   ├── R_integration.R        # R integration wrapper
│   └── metrics.py             # Metrics calculation script
├── conda-envs/                # Conda environment specifications
│   ├── scib-cNMF-env.yml      # Python methods environment
│   ├── R-integration-env.yml  # R methods environment
│   └── scib-metrics-env.yml   # Metrics environment
└── README.md                  # This file
```

## Conda Environments

The pipeline uses three separate conda environments to avoid dependency conflicts:

### 1. `scib-cNMF-env` (Python Integration)
Contains Python-based integration methods:
- scIB module
- Scanorama
- BBKNN
- scVI/scANVI
- scGen
- Combat (via Scanpy)

### 2. `R-integration-env` (R Integration)
Contains R-based integration methods:
- Seurat (CCA and RPCA)
- Harmony
- LIGER
- FastMNN
- Conos

### 3. `scib-metrics-env` (Metrics)
Contains dependencies for computing scIB metrics:
- scIB Python module
- Metrics computation tools

## Execution

The pipeline supports local execution environments through Nextflow:

### Local Execution (default)

```bash
nextflow run main.nf --input data/mydata.h5ad [other params]
```

## Configuration

### Custom Resource Allocation

Modify `nextflow.config` to adjust resources per process:

```groovy
process {
    withName: PREPROCESSING {
        cpus = 4
        memory = '16 GB'
        time = '2h'
    }
    withName: INTEGRATION_PYTHON {
        cpus = 8
        memory = '32 GB'
        time = '4h'
    }
}
```

## Integration Methods

### Python-based Methods

| Method | Description | Reference |
|--------|-------------|-----------|
| **Scanorama** | Panorama-based integration for large datasets | Hie et al., 2019 |
| **BBKNN** | Batch balanced k-nearest neighbors | Polański et al., 2020 |
| **scVI** | Variational inference-based deep learning | Lopez et al., 2018 |
| **scANVI** | Semi-supervised variant of scVI | Xu et al., 2021 |
| **Combat** | Classical batch correction method | Johnson et al., 2007 |
| **scGen** | Generative modeling for perturbation prediction | Lotfollahi et al., 2019 |

### R-based Methods

| Method | Description | Reference |
|--------|-------------|-----------|
| **Seurat CCA** | Canonical correlation analysis | Stuart et al., 2019 |
| **Seurat RPCA** | Reciprocal PCA integration | Stuart et al., 2019 |
| **Harmony** | Fast integration with linear models | Korsunsky et al., 2019 |
| **LIGER** | Integrative non-negative matrix factorization | Welch et al., 2019 |
| **FastMNN** | Fast mutual nearest neighbors | Haghverdi et al., 2018 |
| **Conos** | Joint graph-based integration | Barkas et al., 2019 |

## Metrics

The pipeline computes comprehensive metrics for evaluating integration quality:

### Batch Correction Metrics
- **kBET**: k-nearest neighbor batch effect test
- **LISI**: Local inverse Simpson's index
- **ASW (batch)**: Average silhouette width for batch mixing

### Biological Conservation Metrics
- **ARI**: Adjusted Rand index
- **NMI**: Normalized mutual information
- **ASW (label)**: Average silhouette width for biological labels
- **Isolated labels**: F1 score for rare cell type detection

## Contact

For questions, issues, or suggestions:
- Open an issue on [GitHub](https://github.com/taniaG02/scib-nextflow/issues)
- Contact: tania.gonzalo@externo.cnic.es

## References

Complete list of method references available in the [original scIB publication](https://www.nature.com/articles/s41592-021-01336-8).
Hie et al., 2019
Polański et al., 2020
... complete with link...