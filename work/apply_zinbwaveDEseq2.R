#' run_zinbwaveDESeq2
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

run_zinbwaveDESeq2 <- function(e) {

  # compute zinbwave weights
  # this is expensive!!!!! BPPARAM=BioCParallel::bpparam()
  zinb <- zinbwave::zinbFit(counts, X = model_matrix, epsilon = 1e12)

  # conrast with zeroWeights function?
  weights <- zinbwave::computeObservationalWeights(zinb, counts)
  dimnames(weights) <- NULL

  # DESeq2 flow
  dse <- DESeq2::DESeqDataSetFromMatrix(counts, colData = meta, design = model_matrix)
  dse <- DESeq2::estimateSizeFactors(dse, type = "poscounts")

  # add weights to DESeq2 object
  SummarizedExperiment::assays(dse, withDimnames=FALSE)[["weights"]]<-weights

  # resume DESeq2 flow
  dse <- DESeq2::estimateDispersions(dse, minmu = 1e-6, quiet=TRUE)
  dse <- DESeq2::nbinomWaldTest(dse, betaPrior = TRUE, useT = TRUE, df = rowSums(weights) - 2, minmu = 1e-6, quiet=TRUE)

  # collect results
  DESeq2::resultsNames(dse)

  res <- DESeq2::results(dse, contrast=contrast, name="conditionCTRL")

  # which result to collect?

  return(res)
}
