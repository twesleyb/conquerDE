#' run_zinbwaveDESeq2
#'
#' NOTE: this workflow is aimed at scRNA data 
#'
#' @importFrom zinbwave
#' @importFrom DESeq2


root <- "~/projects/conquerDE"
devtools::load_all(root)
data(L)

counts <- L$counts
meta <- L$meta
contrast <- L$contrast
model_matrix <- L$design

# method=zingeR
# adapted from:
# https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/usoskin/benchmark_usoskin.Rmd

run_zinbwaveDESeq2 <- function(L, alpha=0.05) {

		model_matrix <- model.matrix(~ condition, data=meta)

  # compute zinbwave weights
  # BPPARAM=BiocParallel::MulticoreParam(2))
  nthreads <- BiocParallel::multicoreWorkers()-1
  zinb <- zinbwave::zinbFit(counts, X = model_matrix, epsilon = 1e12, 
							BPPARAM=BiocParallel::MulticoreParam(nthreads))

  # contrast computeObservationalWeights with zeroWeights function?
  weights <- zinbwave::computeObservationalWeights(zinb, counts)
  dimnames(weights) <- NULL

  # DESeq2 flow
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = meta, design = model_matrix)
  #dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = meta, design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")

  # add weights to DESeq2 object
  SummarizedExperiment::assays(dds, withDimnames=FALSE)[["weights"]] <- weights

  # resume DESeq2 flow
  dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6, quiet=TRUE)

  dds <- DESeq2::nbinomWaldTest(dds, betaPrior = TRUE, useT = TRUE, 
								df = rowSums(weights) - 2, minmu = 1e-6, quiet=TRUE)

  # collect results
  DESeq2::resultsNames(dds)

  res <- DESeq2::results(dds, contrast=contrast, name="Intercept")

  # which result to collect?

  return(res)
}
