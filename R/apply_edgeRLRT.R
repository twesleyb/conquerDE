#' apply_edgeRLRT
#'
#' run edgeR LRT analysis
#'
#' @export run_edgeRLRT
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @importFrom edgeR glmFit glmLRT topTags


run_edgeRLRT <- function(L, alpha = 0.05) {

  # collect inputs
  counts <- L$counts
  model_matrix <- L$design
  contrast <- L$contrast
  meta <- L$meta

  # edgeR analysis
  dge <- edgeR::DGEList(counts, samples = meta)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design = model_matrix)
  fit <- edgeR::glmFit(dge, design = model_matrix)
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  tt <- edgeR::topTags(lrt, n = Inf)

  # collect results
  res <- as.data.frame(tt$table)
  res$candidate <- res$FDR < alpha
  return(res)
}
