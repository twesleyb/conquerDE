#' apply_limma
#'
#' run limma DE analysis
#'
#' @export run_limma
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom limma lmFit contrasts.fit eBayes topTable

run_limma <- function(L, trend = FALSE, robust = FALSE, alpha = 0.05) {

  # collect inputs
  meta <- L$meta
  counts <- L$counts
  contrast <- L$contrast
  model_matrix <- L$design

  # perform library and cpm normalization with edgeR
  dge <- edgeR::DGEList(counts, sample = meta)
  dge <- edgeR::calcNormFactors(dge)
  counts_cpm <- edgeR::cpm(dge)

  # create EList object and add log(CPM counts)
  y <- new("EList")
  y$E <- log(counts_cpm + 1)

  # limma workflow
  fit <- limma::lmFit(y, design = model_matrix)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")

  # collect results
  tt$candidate <- tt$"adj.P.Val" < alpha
  return(tt)
}
