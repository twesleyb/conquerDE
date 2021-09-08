#' apply_edgeRqlf
#'
#' run edgeR QLF analysis
#'
#' @export run_edgeRQLF
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @importFrom edgeR glmQLFit glmQLFTest topTags

run_edgeRQLF <- function(L, alpha = 0.05) {

  # collect inputs
  counts <- L$counts
  model_matrix <- L$design
  contrast <- L$contrast
  meta <- L$meta

  # edgeR analysis
  dge <- edgeR::DGEList(counts, samples = meta)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design = model_matrix)
  fit <- edgeR::glmQLFit(dge, design = model_matrix)
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  tt <- edgeR::topTags(qlf, n = Inf)

  # collect results
  res <- as.data.frame(tt$table)
  res$candidate <- res$FDR < alpha
  return(res)
}
