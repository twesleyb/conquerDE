#' apply_wilcoxon
#'
#' @export run_wilcoxon

run_Wilcoxon <- function(L) {

  # root <-  "~/projects/conquerDE"
  # devtools::load_all(root)
  # data(L)

  # collect intputs
  counts <- L$counts
  model_matrix <- L$design
  contrast <- L$contrast
  meta <- L$meta

  # edgeR library and cpm normalization
  dge <- edgeR::DGEList(counts)
  dge <- edgeR::calcNormFactors(dge)
  dm <- edgeR::cpm(dge)

  idx <- 1:nrow(dm)
  names(idx) <- rownames(dm)

  wilcox_p <- sapply(idx, function(i) {
    c(
      foldChange = sapply(split(dm[i, ], meta$condition), mean) %*% contrast,
      pvalue = wilcox.test(dm[i, ] ~ meta$condition)$p.value
    )
  })

  # collect results
  res <- as.data.frame(t(wilcox_p))
  res$FDR <- p.adjust(res$pvalue, method = "BH")
  res$candidate <- res$candidate <- res$FDR < 0.05

  return(res)
}
