# scib-nextflow

A reproducible Nextflow pipeline for single-cell RNA-seq data integration and benchmarking, implementing the scIB (single-cell Integration Benchmarking) framework.

## Overview

This pipeline provides a modular, reproducible implementation of [scIB](https://scib.readthedocs.io/en/latest/) in Nextflow, enabling systematic benchmarking of single-cell integration methods across different datasets and batch correction strategies.

**Pipeline inspired by:** [theislab/scib-pipeline](https://github.com/theislab/scib-pipeline.git) (Luecken et al., 2022)

### Key Features

- **Modular architecture**: Independent processes for preprocessing, integration, and metrics calculation
- **Multiple integration methods**: 
  - Python-based: Scanorama, BBKNN, scVI, scANVI, Combat, scGen
  - R-based: Seurat CCA/RPCA, Harmony, LIGER, FastMNN
- **Comprehensive metrics**: ARI, NMI, ASW, kBET, LISI, and other scIB metrics
- **Flexible execution**: Toggle individual pipeline stages via command-line parameters
- **Reproducibility**: Conda environment management for consistent dependency resolution
- **Automated reporting**: Combines metrics across methods and generates comparison plots

## Pipeline Stages

The pipeline consists of four main stages:

1. **Preprocessing** (Python/Scanpy)
   - Batch-aware highly variable gene (HVG) selection using scIB implementation
   - CPM normalization
   - Batch-aware scaling (optional)

2. **Seurat Conversion** (Python)
   - Convert preprocessed AnnData (.h5ad) to Seurat format (.rds) for R-based methods

3. **Integration**
   - **Python methods**: Scanorama, BBKNN, scVI, scANVI, scGen, Combat
   - **R methods**: Seurat-CCA, Seurat-RPCA, Harmony, LIGER, FastMNN, Conos

4. **Metrics Computation** (Python/scIB)
   - Batch correction metrics: kBET, LISI, ASW (batch)
   - Biological conservation metrics: ARI, NMI, ASW (label), silhouette score
   - Combined scores and visualizations

## Installation

### Requirements

- [Nextflow](https://www.nextflow.io/) (>=21.04.0)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- Git

### Quick Start

1. Clone the repository:
```bash
git clone https://github.com/taniaG02/scib-nextflow.git
cd scib-nextflow
```

2. Ensure conda environments are available (they will NOT be created automatically - limitation):
```bash
# Nextflow will not create environments from YAML files in conda-envs/
# This is a limitation that must be resolve in the future

# Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash
```

3. Make scripts in `bin/` executable:
```bash
chmod +x bin/*.py
```

4. Test with example data
```bash
nextflow run main.nf --input data/test.h5ad --batch batch --label_key celltype --organism human
```

## Usage

### Basic Usage

Run the complete pipeline from a directory containing your data:

```bash
nextflow run main.nf \
  --input data/mydata.h5ad \
  --batch batch_column \
  --label_key cell_type \
  --organism human \
  --outdir ./results
```

## Parameters

### Required Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Input AnnData file (.h5ad) for preprocessing | Required |
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
| `--outdir` | Output directory for results | `./scib_results` |
| `--generate_plots` | Generate comparison plots | `false` |
| `--reference_metrics` | Reference metrics CSV for comparison | `null` |

### Execution Control Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--run_preprocessing` | Execute preprocessing step | `true` |
| `--run_integration_py` | Execute Python integration methods | `true` |
| `--run_integration_r` | Execute R integration methods | `true` |
| `--run_metrics` | Execute metrics calculation | `true` |


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

#### Execution from different directory

The pipeline is designed to be executed from any working directory. Results and logs will be saved in your current location:

```bash
# Navigate to your desired working directory
cd /path/to/your/project/data

# Run the pipeline (repository can be anywhere)
nextflow run /path/to/scib-nextflow/main.nf \
    --input dataset.h5ad \
    --outdir ./results
```
#### Compute all steps and generate plots:

```bash
nextflow run main.nf \
  --input data/*.h5ad \
  --outdir ./results \
  --batch batch \
  --label_key final_annotation \
  --generate_plots true
```

#### Compare with reference results:
```bash
nextflow run main.nf \
    --input data/*.h5ad \
    --reference_metrics /path/to/original_combined_metrics.csv \
    --generate_plots true
```

#### Resume a previous run:
```bash
nextflow run /path/to/scib-nextflow/main.nf -resume
```


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
│   ├── scanorama_metrics.csv      # Metrics for Scanorama
│   ├── harmony_metrics.csv        # Metrics for Harmony
│   └── ...                        # Metrics for other methods
└── Plots/
    ├── barplot_ARI.png
    ├── barplot_all_metrics.png
    ├── heatmap_diff_absolute.png
    └── ...
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

## Integration Methods

### Python-based Methods

| Method | Description | Reference |
|--------|-------------|-----------|
| **Scanorama** | Panorama-based integration for large datasets | [Hie et al., 2019](https://doi.org/10.1038/s41587-019-0113-3) |
| **BBKNN** | Batch balanced k-nearest neighbors | [Polański et al., 2020](https://doi.org/10.1093/bioinformatics/btz625) |
| **scVI** | Variational inference-based deep learning | [Gayoso et al., 2022](https://doi.org/10.1038/s41587-021-01206-w) |
| **scANVI** | Semi-supervised variant of scVI | Xu et al., 2021 |
| **Combat** | Classical batch correction method (lineal) | Johnson et al., 2007 |
| **scGen** | Generative modeling for batch correction | Lotfollahi et al., 2019 |

### R-based Methods

| Method | Description | Reference |
|--------|-------------|-----------|
| **Seurat CCA** | Canonical correlation analysis | [Stuart et al., 2019](https://doi.org/10.1016/j.cell.2019.05.031)|
| **Seurat RPCA** | Reciprocal PCA integration | [Stuart et al., 2019](https://doi.org/10.1016/j.cell.2019.05.031) |
| **Harmony** | Fast integration with linear models | [Korsunsky et al., 2019](https://doi.org/10.1038/s41592-019-0619-0) |
| **LIGER** | Integrative non-negative matrix factorization | [Welch et al., 2019](https://doi.org/10.1016/j.cell.2019.05.006) |
| **FastMNN** | Fast mutual nearest neighbors | [Haghverdi et al., 2018](https://doi.org/10.1038/nbt.4091) |
| **Conos** | Joint graph-based integration | Barkas et al., 2019 |

## Benchmarking Metrics

The pipeline calculates metrics for evaluating integration quality defined by the scIB framework ([Luecken et al., 2022](https://doi.org/10.1038/s41592-021-01336-8)):

**Batch correction metrics:**
- kBET: k-nearest neighbor batch effect test
- LISI: Local Inverse Simpson's Index
- ASW (batch): Average Silhouette Width for batch mixing

**Biological conservation metrics:**
- ARI: Adjusted Rand Index
- NMI: Normalized Mutual Information
- ASW (label): Average Silhouette Width for biological labels

## Current limitations:
- Conda environments must be created and managed manually
- Full portability across all systems may require environment adjustments

## Contact

For questions, issues, or suggestions:
- Open an issue on [GitHub](https://github.com/taniaG02/scib-nextflow/issues)
- Contact: tania.gonzalo@externo.cnic.es

## Acknowledgments

- [Theis Lab](https://github.com/theislab) for the original scIB framework
- Centro Nacional de Investigaciones Cardiovasculares Carlos III (CNIC)
- MSc in Bioinformatics and Computational Biology program

## References

- Ewels, P.A., Peltzer, A., Fillinger, S. et al. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol* **38**, 276–278 (2020). [https://doi.org/10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

- Gayoso, A., Lopez, R., Xing, G. et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol* **40**, 163–166 (2022). [https://doi.org/10.1038/s41587-021-01206-w](https://doi.org/10.1038/s41587-021-01206-w)

- Haghverdi, L., Lun, A.T.L., Morgan, M.D. et al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nat Biotechnol* **36**, 421–427 (2018). [https://doi.org/10.1038/nbt.4091](https://doi.org/10.1038/nbt.4091)

- Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. *Nat Biotechnol* **37**, 685–691 (2019). [https://doi.org/10.1038/s41587-019-0113-3](https://doi.org/10.1038/s41587-019-0113-3)

- Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* **16**, 1289–1296 (2019). [https://doi.org/10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)

- Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods* **19**, 41–50 (2022). [https://doi.org/10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8)

- Polański, K., Young, M.D., Miao, Z. et al. BBKNN: fast batch alignment of single cell transcriptomes. *Bioinformatics* **36**(3), 964–965 (2020). [https://doi.org/10.1093/bioinformatics/btz625](https://doi.org/10.1093/bioinformatics/btz625)

- Stuart, T., Butler, A., Hoffman, P. et al. Comprehensive Integration of Single-Cell Data. *Cell* **177**, 1888–1902.e21 (2019). [https://doi.org/10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)

- Welch, J.D., Kozareva, V., Ferreira, A. et al. Single-Cell Multi-omic Integration Compares and Contrasts Features of Brain Cell Identity. *Cell* **177**, 1873–1887.e17 (2019). [https://doi.org/10.1016/j.cell.2019.05.006](https://doi.org/10.1016/j.cell.2019.05.006)

- Zappia, L., Ramírez-Suástegui, C., Kfuri-Rubens, R. et al. Feature selection methods affect the performance of scRNA-seq data integration and querying. *Nat Methods* **22**, 834–844 (2025). [https://doi.org/10.1038/s41592-025-02624-3](https://doi.org/10.1038/s41592-025-02624-3)
