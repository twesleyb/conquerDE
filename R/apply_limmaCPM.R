#' apply_limmaCPM
#'
#' @export run_limmatrend
#'
#' @importFrom limma EList
#' @importFrom edgeR DGEList calcNormFactors

run_limmaCPM <- function(L, alpha = 0.5) {

  # library(conquerDE)
  # root <- "~/projects/conquerDE"
  # devtools::load_all(root)
  # data(L)

  # colllect input
  model_matrix <- L$design
  contrast <- L$contrast
  counts <- L$counts
  meta <- L$meta

  # cpm normalization
  dge <- edgeR::DGEList(counts, sample = meta)
  dge <- edgeR::calcNormFactors(dge)
  counts_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 3)

  # init limma object
  y <- new("EList")
  y$E <- counts_cpm

  # limma flow
  fit <- limma::lmFit(y, design = model_matrix)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  res <- limma::topTable(fit, n = Inf, adjust.method = "BH")

  # collect results
  res$candidate <- res$adj.P.Val < alpha
  return(res)
}
