#' apply_edgeRzinbwave 
#'
#' run edgeR with zinbwave weights and LRT for DE (aka zinger method)
#'
#' see the zinger paper \insertRef{VandenBerge2018}{conquerDE}
#'
#' @export run_edgeRzinbwave
#'
#' @importFrom zinbwave zinbFit computeObservationalWeights glmWeightedF
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit

run_edgeRzinbwave <- function(L, alpha = 0.05) {

  # collect inputs
  counts <- L$counts
  model_matrix <- L$design
  meta <- L$meta
  contrast <- L$contrast

  # rm rows (genes) with zero counts
  drop <- rowSums(counts) == 0
  sum(drop)
  counts <- counts[!drop, ]

  # fit zinbwave model
  nthreads <- BiocParallel::multicoreWorkers() - 1
  zinb <- zinbwave::zinbFit(counts,
    X = model_matrix, epsilon = 1e12,
    BPPARAM = BiocParallel::MulticoreParam(nthreads)
  )

  # collect weights
  weights <- zinbwave::computeObservationalWeights(zinb, counts)

  # edgeR flow
  dge <- edgeR::DGEList(counts, samples=meta)
  dge <- edgeR::calcNormFactors(dge)

  # add weights to edgeR object
  dge$weights <- weights

  # resume edgeR flow
  dge <- edgeR::estimateDisp(dge, design = model_matrix)
  fit <- edgeR::glmFit(dge, design = model_matrix)

  # zinbwave LRT with weights 
  # NOTE: why not edgeR glmLRT here?
  lrt <- zinbwave::glmWeightedF(fit, contrast = contrast, independentFiltering = TRUE)

  # collect results
  res <- lrt$table
  res$candidate <- res$padjFilter < alpha

  return(res)
}
