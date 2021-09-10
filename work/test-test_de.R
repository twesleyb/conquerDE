#!/usr/bin/env Rscript

# Hypothesis: performance of edgeR and glmGamPoi depends upon fraction zeros
# * edgeR is faster then glmGamPoi when fraction zeros is low
# * glmGamPoi is faster then edgeR when fraction zeros is high
# * performace is about the same when fraction zeros is modest:0.3
# * Q: what is typical fraction zeros?
# more robust test
# how to properly titrate fraction zeros in nbinom distribution

library(glmGamPoi)

n_cols <- 10 # 10
n_rows <- 30 # 30
dispersion <- 40 # 4
success <- 0.99 # 0.3
Y <- matrix(rnbinom(n = n_rows * n_cols, mu = dispersion, size = success), nrow = n_rows, ncol = n_cols)
fraction_zeros <- sum(Y == 0) / length(Y) # ~0.40 - 0.50
print(fraction_zeros)
annot <- data.frame(
  group = sample(c("A", "B"), size = n_cols, replace = TRUE),
  cont1 = rnorm(n_cols), cont2 = rnorm(n_cols, mean = 30)
)
design <- model.matrix(~ group + cont1 + cont2, data = annot)
LT <- limma::makeContrasts(-groupB, levels = design)
L <- list(counts = Y, design = design, contrast = LT, meta = annot)
fit <- glm_gp(Y, design = design)
res <- test_de(fit, contrast = LT)

bench::mark(
  glmGamPoi = {
    conquerDE(L, method = "glmGamPoi")
  },
  edgeR = {
    conquerDE(L, method = "edgeRQLF")
  }, check = FALSE, min_iterations = 3
)


fZ <- function(X) { sum(X == 0) / length(X) }

fZ(filt_counts)

x=generateData(dispersion=4,success=2)[["counts"]]
fZ(x)

generateData <- function(n_cols = 10, n_rows = 30, dispersion = 4, success = 0.3) {
  Y <- matrix(rnbinom(n = n_rows * n_cols, mu = dispersion, size = success), nrow = n_rows, ncol = n_cols)
  annot <- data.frame(
    group = sample(c("A", "B"), size = n_cols, replace = TRUE),
    cont1 = rnorm(n_cols), cont2 = rnorm(n_cols, mean = 30)
  )
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  suppressWarnings({
    LT <- limma::makeContrasts(-groupB, levels = design)
  })
  L <- list(counts = Y, design = design, contrast = LT, meta = annot)
  return(L)
}


success <- seq(0.1, 0.9, by = 0.1)
results <- list()

for (i in seq_along(success)) {
  L <- generateData(n_cols = 10, n_rows = 30, dispersion = 4, success = success[i])
  bm <- bench::mark(
    glmGamPoi = {
      conquerDE(L, method = "glmGamPoi")
    },
    edgeR = {
      conquerDE(L, method = "edgeRQLF")
    },
    check = FALSE, min_iterations = 3
  )
x <- bm %>% dplyr::pull(as.numeric(median), name=expression)
results[[i]] <- c(x,fractionZeros=fZ(L[["counts"]]))
}

df = do.call(rbind,results)

df = dplyr::bind_rows(results) %>% dplyr::mutate(delta=edgeR/glmGamPoi)
df
