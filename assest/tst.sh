nextflow run main.nf \
  --input "data/input.h5ad" \
  --output "results_test" \
  --batch "batch_var" \
  --labelkey "cell_type" \
  --organism "human" \
  -with-conda
