#' apply_DESeq2
#'
#' apply DESeq2 method
#'
#' @export run_DESeq2
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results

run_DESeq2 <- function(L, alpha = 0.05) {

  # collect input data from list
  counts <- L$counts
  model_matrix <- L$design
  contrast <- L$contrast
  meta <- L$meta

  # DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = model_matrix
  )

  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # NOTE: padj method defaults to "BH"
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)

  # other methods:
  # DESeq2::nbinomLRT
  # DESeq2::nbionomWaldTest

  # colnames(res$table)
  # [1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"
  # [5] "pvalue"         "padj"

  # collect results
  res <- as.data.frame(res)
  res$candidate <- res$padj < alpha

  return(res)
}
