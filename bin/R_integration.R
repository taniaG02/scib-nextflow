#!/usr/bin/env Rscript

## Descripton: Callin methods for scRNA-seq data integration implemented in R
# Date: 2025-02-03


library("optparse")
library("Seurat")
library("SingleCellExperiment")
library("zellkonverter")
# library("rlang")

# Change directory to script dir
# getScriptPath <- function() {
#   cmd.args <- commandArgs()
#   m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#   script.dir <- dirname(regmatches(cmd.args, m))
#   if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
#   if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
#   return(script.dir)
# }
# 
# setwd(getScriptPath())


options.list <- list(
  make_option(
    c("-i", "--input"), type = "character", default = NA, 
    help = "Input Seurat object (RDS file)"
  ),
  make_option(
    c("-o", "--output"), type = "character", 
    default = NA, help = "Output Seurat object"
  ),
  make_option(
    c("-b", "--batch"), type = "character", 
    default = NA, help = "Batch variable"
  ),
  make_option(
    c("-m", "--method"), 
    type = "character", 
    default = NULL, 
    help = "Integration method to use"
  ),
  make_option(
    c("-v", "--hvg"), type = "character", default = NA, 
    help = "hvg list for Seurat"
  ),
  make_option(
    c("-s", "--functions"), type = "character", default = NA, 
    help = "Dir path for integration functions"
  )
)

opts <- OptionParser(
  usage = "usage: %prog [options]",
  option_list = options.list, 
  add_help_option = TRUE,
  prog = "R-integration-call.R", 
  description = paste(
    "R script for calling integration methods implemented in R."
  ),
  epilogue = ""
)

opt <- parse_args(
  opts,  
  args = commandArgs(trailingOnly = TRUE),
  print_help_and_exit = TRUE,
  positional_arguments = FALSE
)

message("Running R-based integration methods...")
message("\t - Input: ", opt$input)
message("\t - Output: ", opt$output)
message("\t - Batch: ", opt$batch)
message("\t - Method: ", opt$method)
message("\t - HVG: ", opt$hvg)
message("\t - Functions: ", opt$functions)

## from a parameter, not optimal, but anyway
source(file.path(opt$functions, "R_integration_func.R"))

seurat.obj <- readRDS(opt$input)
# seurat.obj@assays$RNA <- seurat.obj@assays$originalexp
# seurat.obj@assays$originalexp <- NULL

if (opt$method == "Seurat-CCA") {
  
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names = FALSE)
  # } else {
  hvg <- rownames(seurat.obj@assays$originalexp)
  # }
  out <- runSeuratCCA(seurat.obj, opt$b, hvg)
  
} else if (opt$method == "Seurat-RPCA") {
  
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names = FALSE)
  # } else {
  hvg <- rownames(seurat.obj@assays$originalexp)
  # }
  out <- runSeuratRPCA(seurat.obj, opt$batch, hvg)
  
} else if (opt$method == "conos") {
  
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names=FALSE)
  #   seurat.obj <- subset(seurat.obj, features = hvg)
  # } else {
  hvg <- rownames(seurat.obj@assays$originalexp)
  # }
  out <- runConos(seurat.obj, opt$batch)
  
  saveConos(out, opt$output)
  out <- as.Seurat(out)
  

} else if(opt$method == 'harmony'){
  
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names = FALSE)
  #   seurat.obj <- subset(seurat.obj, features = hvg)
  # } else {
  hvg <- rownames(seurat.obj@assays$originalexp)
  # }
  out <- runHarmonyR(seurat.obj, opt$batch)
  
} else if (opt$method == "liger") {
  
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names=FALSE)
  # } else {
  hvg <- rownames(seurat.obj)
  # }
  out <- runLiger(seurat.obj, opt$batch, hvg)
  
} else if (opt$method == "fastmnn") {
  # if (!is.na(opt$hvg)) {
  #   hvg <- unlist(readRDS(opt$hvg), use.names = FALSE)
  #   seurat.obj <- subset(seurat.obj, features = hvg)
  # }
  out <- runFastMNN(seurat.obj, opt$batch)
  
}

## saving results
saveRDS(out, file.path(opt$output, paste0(opt$method, "-integrated.rds")))

#TODO: in case it is a Seurat obejct (all functions should return a Seruat object)
# there is some code to do it using scib, in case it could be useful: 
# <https://github.com/theislab/scib-pipeline/blob/e97631a478e063883cee5db0dcebbf42d8268ab0/scripts/integration/runPost.py>


# DefaultAssay(out) <- "RNA"
out.obj <- as.SingleCellExperiment(out)
zellkonverter::writeH5AD(
  out.obj, 
  file = file.path(opt$output, paste0(opt$method, "-integrated.h5ad"))
  #   file.path(
  #   paste0(tools::file_path_sans_ext(opt$output), ".h5ad")
  # )
)
