set -e

echo ">>>> Module 2: Integration of data - Python methods"

## hardcoded, this could have been passed as an argument
# METHODS_PY=('scanorama' 'bbknn' 'scvi' 'combat' 'NMFusion-CPMs' 'NMFusion-CPMs-NMF-space' 'NMF-counts' 'NMF-counts-NMF-space')
# METHODS_PY=('scanorama' 'bbknn' 'scvi' 'combat' 'NMFusion-CPMs' 'NMFusion-counts')

# Get parameter file from command line argument
METHODS_STRING=$1

if [ -z "$METHODS_STRING" ]; then
    echo "Error: Methods string not provided"
    exit 1
fi

IFS=',' read -ra METHODS_ARRAY <<< "$METHODS_STRING"


echo "Methods to be processed: ${METHODS_ARRAY[*]}"


taskid=$SGE_TASK_ID
let "taskid -= 1" || taskid=0

echo ">>> Task ID: $SGE_TASK_ID / $taskid"
echo ">>> Method: ${METHODS_ARRAY[$taskid]}"

date +"%F %X" || exit 100 &&

echo ">>> COMMAND:" || exit 100 &&
echo "source $HOME/miniconda3.cluster.source || exit 100 &&
conda activate scib-cNMF-env || exit 100 &&
python $SOFTWARE_DIR/src/py_integration.py -i $PREPROCESSED_OUTPUT/adata-preprocessed.h5ad -o $INTEGRATED_OUTPUT -b $BATCH -m ${METHODS_ARRAY[$taskid]} -v $HVGS || exit 100 &&" || exit 100 &&

source $HOME/miniconda3.cluster.source || exit 100 &&
conda activate scib-cNMF-env || exit 100 &&
python $SOFTWARE_DIR/src/py_integration.py -i $PREPROCESSED_OUTPUT/adata-preprocessed.h5ad -o $INTEGRATED_OUTPUT -b $BATCH -m ${METHODS_ARRAY[$taskid]} -v $HVGS || exit 100 &&
echo "DONE" || exit 100
