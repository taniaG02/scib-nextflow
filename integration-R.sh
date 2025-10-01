set -e

echo ">>>> Module 2: Integration of data - R methods"

## hardcoded, this could have been passed as an argument
# METHODS_R=('Seurat-CCA' 'Seurat-RPCA' 'harmony' 'liger' 'fastmnn')

METHODS_STRING=$1

if [ -z "$METHODS_STRING" ]; then
    echo "Error: Methods string not provided"
    exit 1
fi

IFS=',' read -ra METHODS_ARRAY <<< "$METHODS_STRING"

taskid=$SGE_TASK_ID
let "taskid -= 1" || taskid=0

echo ">>> Task ID: $SGE_TASK_ID / $taskid"
echo ">>> Method: ${METHODS_ARRAY[$taskid]}"
date +"%F %X" || exit 100 &&

echo ">>> COMMAND:" || exit 100 &&
echo "source $HOME/miniconda3.cluster.source || exit 100 &&
conda activate R-integration-env || exit 100 &&
Rscript $SOFTWARE_DIR/src/R_integration.R -i $PREPROCESSED_OUTPUT/seurat-preprocessed.rds -o $INTEGRATED_OUTPUT -b $BATCH -m ${METHODS_ARRAY[$taskid]} -v $HVGS -s $SOFTWARE_DIR/src/R_integration_func.R || exit 100 &&" || exit 100 &&

source $HOME/miniconda3.cluster.source || exit 100 &&
conda activate R-integration-env || exit 100 &&
Rscript $SOFTWARE_DIR/src/R_integration.R -i $PREPROCESSED_OUTPUT/seurat-preprocessed.rds -o $INTEGRATED_OUTPUT -b $BATCH -m ${METHODS_ARRAY[$taskid]} -v $HVGS -s $SOFTWARE_DIR/src/R_integration_func.R || exit 100 &&
echo "DONE" || exit 100
