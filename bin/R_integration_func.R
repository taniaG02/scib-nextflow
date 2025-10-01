## Description: Functions for integration methods implemented in R
# Date: 2025-02-03

## loading packages [assuming they are available]
library("Seurat")
library("SeuratWrappers")
library("dplyr")
library("conos")
library("harmony")
library("rliger")
library("batchelor")

runSeuratCCA <- function(
	seurat.obj, 
	batch, 
	hvg = 2500
) {
	batch_list <- SplitObject(seurat.obj, split.by = batch)

	anchors <- FindIntegrationAnchors(
		object.list = batch_list,
		anchor.features = hvg,
		scale = T,
		l2.norm = T,
		dims = 1:30,
		k.anchor = 5,
		k.filter = 200,
		k.score = 30,
		max.features = 200,
		eps = 0
	)
	seurat.integrated <- IntegrateData(
		anchorset = anchors,
		new.assay.name = "integrated",
		features = NULL,
		features.to.integrate = NULL,
		dims = 1:30,
		k.weight = 80,
		weight.reduction = NULL,
		sd.weight = 1,
		sample.tree = NULL,
		preserve.order = F,
		eps = 0,
		verbose = T
	)

	return(seurat.integrated)
}

runSeuratRPCA <- function(
	seurat.obj, 
	batch, 
	hvg = 2500
) {
	
	batch_list <- SplitObject(seurat.obj, split.by = batch)
	#features <- SelectIntegrationFeatures(batch_list)
	batch_list <- lapply(
			X = batch_list, 
			FUN = function(x) {
				x <- ScaleData(x, features = hvg)
				x <- RunPCA(x, features = hvg)
				return(x)
		}
	)
	anchors <- FindIntegrationAnchors(
		object.list = batch_list,
		anchor.features = hvg,
		scale = T,
		l2.norm = T,
		dims = 1:30,
		k.anchor = 5,
		k.filter = 200,
		k.score = 30,
		reduction = "rpca",
		max.features = 200,
		eps = 0
	)

	seurat.integrated <- IntegrateData(
		anchorset = anchors,
		new.assay.name = "integrated",
		features = NULL,
		features.to.integrate = NULL,
		dims = 1:30,
		k.weight = 80,
		weight.reduction = NULL,
		sd.weight = 1,
		sample.tree = NULL,
		preserve.order = F,
		eps = 0,
		verbose = T
	)

	return(seurat.integrated)
}


# func_profiler <- function(expr, chunksize=20000, filename='timing.out', prof.interval=0.02) {
# 	      Rprof(filename, memory.profiling=T, interval=prof.interval)
# 	      res = expr
# 	      Rprof(NULL)
# 	      t = summaryRprof(filename, chunksize=chunksize, memory="both")$sampling.time
# 	      mem = max(summaryRprof(filename, chunksize=chunksize, memory="both")$by.total$mem.total)
# 	      return(list(results=res, time=t, memory=mem))
# }
# Example call:
#   sobj = load_seurat_object('small_test.RDS')
#   out = func_profiler(runSeurat(sobj, "batch"))
#   out$results is results
#   out$time is timing
#   out$memory is memory use

seuratWorkflow <- function(
	seurat.obj, 
	hvg = 2500,
	vars.to.regress = NULL, 
	n.pcs = 50,
	verbose = TRUE
) {
#   if (verbose) 
#     message("Running Seurat v3 workflow")

  seurat.obj <- seurat.obj %>% FindVariableFeatures(nfeatures = hvg, verbose = verbose) %>%
    ScaleData(verbose = verbose) %>% 
    RunPCA(npcs = n.pcs, verbose = verbose)

  return(seurat.obj)
}

runConos <- function(seurat.obj, batch) {

  batch_list <- SplitObject(seurat.obj, split.by=batch)
 	pp <- lapply(batch_list, seuratWorkflow)

	con <- Conos$new(pp)
	con$buildGraph(space = "genes")
	con$findCommunities()
	con$embedGraph(method = "UMAP")

	return(con)
}

saveConos <- function(con, outdir) {
	# dir.create(outdir)
	hdf5_filename <- paste0(tools::file_path_sans_ext(basename(outdir)), ".h5")
	outdir.path <- dirname(outdir)
	saveConosForScanPy(
	  con, 
	  output.path = outdir.path,
	  hdf5_filename = hdf5_filename,
	  pseudo.pca = TRUE, pca = TRUE,
	  verbose = TRUE
	)
}

runHarmonyR <- function(seurat.obj, batch) {
  seurat.obj <- seurat.obj %>% ScaleData() %>% 
    RunPCA(features = rownames(seurat.obj)) %>% 
    RunHarmony(batch)
  
  seurat.obj[['X_emb']] <- seurat.obj[['harmony']]

  return(seurat.obj)
}

## TODO: change according to this tutorial: <https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html>
runLiger <- function(
  seurat.obj, 
  batch, 
  hvg, 
  k = 30, 
  res = 0.5, 
  small.clust.thresh = 20
) {
  # only counts is converted to liger object. To pass our own normalized data,
  # store it in the "counts" slot
#   if (is.null(seurat.obj@assays$RNA)) {
#     # Seurat v4
#     data <- GetAssayData(seurat.obj, slot = "data")
#     SetAssayData(seurat.obj, slot = "counts", new.data = as.matrix(data))
#   } 

  # Create Liger object: it does not work in this vertion (2.1.0). Let's do it manually
#   liger.obj <- seuratToLiger(
#     seurat.obj,
# 	assays.use = "originalexp",
#     combined.seurat = T,
#     meta.var = batch,
#     renormalize = F,
#     remove.missing = F
#   )

  list.data <- SplitObject(seurat.obj, split.by = batch) %>% lapply(
	function(x) {
	  return(x@assays$originalexp$data)
	}
  )

  liger.obj <- createLiger(list.data)
  liger.obj@datasets <- lapply(
	X = liger.obj@datasets,
	FUN = function(x) {
		x@normData <- x@rawData
		return(x)
	}
  )
  # We only pass nomarlized data, so store it as such
#   liger.obj@norm.data <- liger.obj@raw.data

  # Assign hvgs
#   liger.obj@var.genes <- hvg
  liger.obj <- liger.obj %>% selectGenes(nGenes = 2500)
  liger.obj <- scaleNotCenter(liger.obj) # Can't do our own scaling atm

  # Use tutorial coarse k suggests of 20.
  liger.obj <-  runIntegration(liger.obj, k = k)
  liger.obj <- alignFactors(liger.obj, method = "centroidAlign")
  
  ## change this accordingly: this is an old version of the package
#   liger.obj <- quantileAlignSNF(liger.obj, resolution = res, small.clust.thresh = small.clust.thresh)

  # Store embedding in initial Seurat object
  # Code taken from ligerToSeurat() function from LIGER
  colnames(liger.obj@H.norm) <- paste0("X-emb-", 1:ncol(liger.obj@H.norm))
  colnames(liger.obj@W) <- paste0("X-emb-", 1:ncol(liger.obj@W))
#   inmf.obj <- new(
#     Class = "DimReduc", feature.loadings = t(liger.obj@W),
#     cell.embeddings = liger.obj@H.norm, key = "X_emb"
#   )

  seurat.obj@reductions['X_emb'] <- CreateDimReducObject(
	embeddings = liger.obj@H.norm, key = "X_emb_", 
	loadings = liger.obj@W,
	assay = DefaultAssay(seurat.obj)
  )
#   seurat.obj@reductions['X_emb'] <- inmf.obj

  return(seurat.obj)
}

## TODO: adaopt these functions to Seurat v5
runFastMNN <- function(seurat.obj, batch) {  
  # Seurat v4
  expr <- GetAssayData(seurat.obj, slot = "data")
  sce <- fastMNN(expr, batch = seurat.obj@meta.data[[batch]])
  corrected_data <- assay(sce, "reconstructed")
  # Seurat v4
  seurat.obj@assays$originalexp@data <- as.matrix(corrected_data)
	#   <- SetAssayData(
	#     seurat.obj, slot = "data", new.data = as.matrix(corrected_data)
	#   )
  seurat.obj@reductions['X_emb'] <- CreateDimReducObject(
    reducedDim(sce, "corrected"), 
    key = 'fastmnn_'
  )

  return(seurat.obj)
}
