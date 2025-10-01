set -e

echo ">>>> Module 2: Integration of data - Python methods"

## hardcoded, this could have been passed as an argument
# conos has been removed: the AnnData object is incorrect
# METHODS_ALL=('scanorama' 'bbknn' 'scvi' 'combat' 'Seurat-CCA' 'Seurat-RPCA' 'harmony' 'liger' 'fastmnn' 'NMFusion-CPMs' 'NMFusion-CPMs-NMF-space' 'NMF-counts' 'NMF-counts-NMF-space')
# METHODS_ALL=('NMFusion-CPMs' 'NMFusion-CPMs-NMF-space' 'NMFusion-counts' 'NMFusion-counts-NMF-space')
# METHODS_ALL=('scanorama' 'bbknn' 'scvi' 'combat' 'Seurat-CCA' 'Seurat-RPCA' 'harmony' 'liger' 'fastmnn' 'NMFusion-CPMs' 'NMFusion-CPMs-NMF-space' 'NMFusion-counts' 'NMFusion-counts-NMF-space')

METHODS_STRING=$1

if [ -z "$METHODS_STRING" ]; then
    echo "Error: Methods string not provided"
    exit 1
fi

IFS=',' read -ra METHODS_ARRAY <<< "$METHODS_STRING"

taskid=$SGE_TASK_ID
let "taskid -= 1" || taskid=0

INTEGRATED_OUTPUT_FILE="$INTEGRATED_OUTPUT/${METHODS_ARRAY[$taskid]}-integrated.h5ad"
PREPROCESSED_OUTPUT_FILE="$PREPROCESSED_OUTPUT/adata-preprocessed.h5ad"

echo ">>> Task ID: $SGE_TASK_ID / $taskid"
echo ">>> Method: ${METHODS_ARRAY[$taskid]}"
date +"%F %X" || exit 100 &&

echo ">>> COMMAND:" || exit 100 &&

## harcoded
echo "source $HOME/miniconda3.cluster.source || exit 100 && 
conda activate scib-metrics-env || exit 100 &&
python $SOFTWARE_DIR/src/metrics.py -u $PREPROCESSED_OUTPUT_FILE -i $INTEGRATED_OUTPUT_FILE -o $METRICS_OUTPUT/${METHODS_ARRAY[$taskid]}_metrics.csv -m ${METHODS_ARRAY[$taskid]} -b $BATCH -l $LABEL_KEY --organism $ORGANISM --type 'full' --hvgs $HVGS --verbose || exit 100 &&" || exit 100 &&


source $HOME/miniconda3.cluster.source || exit 100 && 
conda activate scib-metrics-env || exit 100 &&
python $SOFTWARE_DIR/src/metrics.py -u $PREPROCESSED_OUTPUT_FILE -i $INTEGRATED_OUTPUT_FILE -o $METRICS_OUTPUT/${METHODS_ARRAY[$taskid]}_metrics.csv -m ${METHODS_ARRAY[$taskid]} -b $BATCH -l $LABEL_KEY --organism $ORGANISM --type 'full' --hvgs $HVGS --verbose || exit 100
