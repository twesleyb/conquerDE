#' apply_limma
#'
#' run limma DE analysis (trended, robust)
#'
#' @param L$counts - normalized counts e.g. log2(cpm+1) or log2(fpkm)
#'
#' @export run_limma
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom limma EList lmFit contrasts.fit eBayes topTable

run_limma <- function(L, alpha = 0.05) {

  ##library(conquerDE)
  # root <- "~/projects/conquerDE"
  # devtools::load_all(root)
  # data(Satpathy2021_counts) #log2FPKM
  # data(Satpathy2021_meta)
  # counts <- Satpathy2021_counts
  # meta <- Satpathy2021_meta
  # condition <- sapply(strsplit(colnames(Satpathy2021_counts),"_"),"[",1)
  # model_matrix <- model.matrix(~ 0 + condition)
  # contrast <- limma::makeContrasts("conditiontumor-conditionnormal",
  #   							   levels=model_matrix)
  # L <- list("meta" = meta, "counts" = counts, "design"=model_matrix, "contrast"=contrast)
  # trend = TRUE; robust=TRUE; alpha = 0.05

  # collect inputs
  meta <- L$meta
  norm_counts <- L$counts
  contrast <- L$contrast
  model_matrix <- L$design

  # create EList object and add log(CPM counts)
  y <- new("EList")
  y$E <- norm_counts

  # limma workflow
  fit <- limma::lmFit(y, design = model_matrix)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")

  # collect results
  tt$candidate <- tt$"adj.P.Val" < alpha

  return(tt)
}
