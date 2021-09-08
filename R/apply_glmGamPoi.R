#' apply_glmGamPoi
#'
#' run glmGamPoi analysis
#'
#' @export run_glmGamPoi
#'
#' @importFrom glmGamPoi glm_gp test_de

run_glmGamPoi <- function(L, alpha = 0.05) {

  # collect inputs
  counts <- L$counts
  model_matrix <- L$design
  contrast <- L$contrast
  meta <- L$meta

  # glmGamPoi analysis
  fit <- glmGamPoi::glm_gp(counts, design = model_matrix, on_disk = FALSE)
  res <- glmGamPoi::test_de(fit, contrast)

  # collect results
  res$candidate <- res$adj_pval < alpha
  return(res)
}
