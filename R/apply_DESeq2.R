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

  # FIXME: expand to explicitly perform normalization and model fitting steps
  # for more control over the workflow!
  #dds <- DESeq2::estimateSizeFactors(dds,type="poscounts")
  #dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
  #dds <- DESeq2::nbinomWaldTest(dds, betaPrior=TRUE, modelMatrixType="standard", minmu = 1e-6)
  #res <- DESeq2::results(dds, name="condition_stemCell_vs_fibro")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # NOTE: padj method defaults to "BH"
  # NOTE: disable outlier imputation by setting cooksCutoff=FALSE
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)

  # NOTE: single-cell UMI data: suggested nbinomLRT
  # NOTE: other (non-UMI), the Wald test in nbinomWaldTest can be used, with
  # null distribution a t-distribution with degrees of freedom corrected for
  # downweighting. In both cases, we recommend the minimum expected count to be
  # set to a small value (minmu=1e-6). The Wald test in DESEQ2 allows for testing
  #contrasts of the coefficients (VandenBerge2018).

  # NOTE: limitations of default normalization approach: The PHYLOSEQ normalization procedure can now be applied by setting the
  #argument type equal to "poscounts" in the DESEQ2 function estimateSizeFactors. (VandenBerge2018)

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
