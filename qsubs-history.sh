## qsub history

#### Human PBMCs [Run 2]

### 1 - Module 1: preprocessing <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h '2500' -m 1 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-1" -l h_vmem=3G -pe pthreads 4


echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h '2500' -m 1 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-1-Correcting-Seurat" -l h_vmem=3G -pe pthreads 4

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-2"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-2-NMFusion"

## checking H refit
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human' -t NMFusion-CPMs,NMFusion-counts" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-2-NMFusion-H-refit"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-NMFusion"

### 3 - Module 3: Computing metrics <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-3"


echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-3-NMFusion"

## testing H-refit
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-Mod-3-H-refit"

## NMFusion metrics
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/input/PBMCs-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-NMFusion"

#### Human PBMCs and BM

### 1 - Module 1: preprocessing <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h '2500' -m 1 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMC-BM-Mod-1" -l h_vmem=3G -pe pthreads 4

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-Mod-2"

# H refit
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human' -t NMFusion-CPMs,NMFusion-counts" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-Mod-2-H-refit"


## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-NMFusion"

### 3 - Module 3: Computing metrics <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-Mod-3"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-Mod-3-H-refit"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/input/PBMCs-BM-Human-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-PBMCs-BM-NMFusion-Metrics"


#### Mouse BM

### 1 - Module 1: preprocessing <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h '2500' -m 1 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-Mod-1" -l h_vmem=6G -pe pthreads 4

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-Mod-2"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human' -t NMFusion-CPMs,NMFusion-counts" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-Mod-2-NMFusion-H-refit"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-NMFusion"

### 3 - Module 3: Computing metrics <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-Mod-3"

## H refit
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-Mod-3-H-refit"

## NMFusion metrics
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/input/BM-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Mouse-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Mouse-BM-NMFusion-Metrics"


#### Human and mouse PBMCs and BM

### 1 - Module 1: preprocessing -> ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h '2500' -m 1 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMC-BM-Mod-1" -l h_vmem=5G -pe pthreads 6

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-2"

# H-refit
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human' -t NMFusion-CPMs,NMFusion-counts" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-2-H-refit"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion"

### 3 - Module 3: Computing metrics
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-3"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-3-H-refit"

## NMFusion metrics
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion-Metrics"




#### CD8 viral infection atlas Carmona lab

### 1 - Module 1: preprocessing -> ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/input/seurat.def.LCMV.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/output -b 'Study' -h '2500' -m 1 -l 'functional_cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 6 -N "CD8Tcell-Viral-Mod-1" -l h_vmem=5G

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/input/seurat.def.LCMV.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/output -b 'Study' -h '2500' -m 2 -l 'functional_cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "CD8Tcell-Viral-Mod-2"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion"

### 3 - Module 3: Computing metrics
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/input/seurat.def.LCMV.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-Viral/output -b 'Study' -h 2500 -m 3 -l 'functional_cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "CD8Tcell-Viral-Mod-3"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-3-H-refit"

## NMFusion metrics
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion-Metrics"




#### CD8 TILs atlas Carmona lab

### 1 - Module 1: preprocessing -> ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/input/seurat.def.TILs.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/output -b 'Study' -h '2500' -m 1 -l 'funcitonal.cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 6 -N "CD8Tcell-TILs-Mod-1" -l h_vmem=5G

### 2 - Module 2: Integration of data <- ok
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/input/seurat.def.TILs.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/output -b 'Study' -h '2500' -m 2 -l 'funcitonal.cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "CD8Tcell-TILs-Mod-2"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/input/seurat.def.TILs.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/output -b 'Study' -h '2500' -m 2 -l 'funcitonal.cluster' -s 'mouse' -t NMFusion-CPMs,NMFusion-counts" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "CD8Tcell-TILs-Mod-2-NMF"


## reminder: probably I will need to change the class of the X matrix: 
# from scipy.sparse import csc_matrix
# csr = csc.tocsr()

# relaunching R methods

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/input/seurat.def.TILs.mouse.counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/CD8Tcell-Map-TILs/output -b 'Study' -h '2500' -m 3 -l 'functional.cluster' -s 'mouse'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "CD8Tcell-TILs-Mod-3"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human' -t Seurat-CCA,Seurat-RPCA,harmony,liger,fastmnn" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-2-R-methods"

## NMFusion [implemented after the previous call] 
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 2 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion"

### 3 - Module 3: Computing metrics
echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-3"

echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human' -t NMFusion-CPMs-H-refit,NMFusion-counts-H-refit" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-Mod-3-H-refit"

## NMFusion metrics
# echo "bash /data_lab_DSM/Projects_3/scib-integration-pipeline/main.sh -i /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/input/PBMCs-BM-Human-Mouse-raw-counts.h5ad -o /data_lab_DSM/Projects_3/Results-Pipeline/Human-Mouse-PBMCs-BM/output -b 'batch' -h 2500 -m 3 -l 'final_annotation' -s 'human'" | qsub -P DSM -A "LAB_DSM" -pe pthreads 1 -N "Human-Mouse-PBMCs-BM-NMFusion-Metrics"


