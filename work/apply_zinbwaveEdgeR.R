runZinbwaveEdgeR <- function(e){
		# input is object islam
  library(zinbwave)
  library(edgeR)

  condition = pData(e)$condition
  design <- model.matrix(~ condition)

  # compute zinbwave weights
  zinb <- zinbwave::zinbFit(exprs(e), X = design, epsilon = 1e12)
  weights <- zinbwave::computeObservationalWeights(zinb, exprs(e))

  d <- edgeR::DGEList(exprs(e))
  d <- edgeR::calcNormFactors(d)
  d$weights <- weights
  d=edgeR::estimateDisp(d, design)
  fit=edgeR::glmFit(d,design)
  lrt=zinbwave::glmWeightedF(fit,coef=2, independentFiltering = TRUE)

  pvals = lrt$table$PValue

  list(pvals = pvals, padj = lrt$table$padjFilter,
       logfc = lrt$table$logFC)
}

