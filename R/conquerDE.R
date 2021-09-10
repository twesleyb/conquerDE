#' conquerDE
#'
#' run DE analysis ala \insertRef{Soneson2018}{conquerDE}
#'
#' @export conquerDE

conquerDE <- function(L, method = c("edgeRQLF", "edgeRLRT", "DESeq2", "glmGamPoi", "limma", "limmatrend"), alpha = 0.05) {
  method <- match.arg(method)
  switch(method,
    # DESeq2-based methods
    DESeq2 = {
      run_DESeq2(L, alpha)
    },
    # edgeR-based methods
    edgeRQLF = {
      run_edgeRQLF(L, alpha)
    },
    edgeRLRT = {
      run_edgeRLRT(L, alpha)
    },
    # glmGamPoi-based methods
    glmGamPoi = {
      run_glmGamPoi(L, alpha)
    },
    # limma-based methods
    limma = {
      run_limma(L, alpha)
    },
    limmatrend = {
      run_limma(L, trend = TRUE, robust = TRUE, alpha)
    }
  )
}
