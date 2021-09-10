#!/usr/bin/env Rscript

# prepare the R env
root <- "~/projects/conquerDE"
devtools::load_all(root)

# load the data
data(meta)
data(filt_counts)

# create design matrix and contrast
model_matrix <- model.matrix(~ 0 + condition, data = meta)
LT <- limma::makeContrasts("conditionIgAN-conditionCTRL", levels = model_matrix)

L <- list(counts = filt_counts, design = model_matrix, contrast = LT, meta = meta)

# benchmark conquerDE
bench::mark(
  conquerDE_glmGamPoi = {
    conquerDE(L, method = "glmGamPoi")
  },
  conquerDE_edgeR = {
    conquerDE(L, method = "edgeRQLF")
  }, check = FALSE, min_iterations = 3
)


# Benchmark adapted from glmgampoi
bench::mark(
  glmGamPoi_in_memory = {
    glmGamPoi_fit <- glmGamPoi::glm_gp(filt_counts, design = model_matrix, on_disk = FALSE)
    glmGamPoi_res <- glmGamPoi::test_de(glmGamPoi_fit, contrast = LT)
  }, glmGamPoi_on_disk = {
    glmGamPoi_fit <- glmGamPoi::glm_gp(filt_counts, design = model_matrix, on_disk = TRUE)
  }, DESeq2 = suppressMessages({
    dds <- DESeq2::DESeqDataSetFromMatrix(filt_counts, colData = meta, design = model_matrix)
    dds <- DESeq2::estimateSizeFactors(dds, "poscounts")
    dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
  }), edgeR = {
    edgeR_data <- edgeR::DGEList(filt_counts)
    edgeR_data <- edgeR::calcNormFactors(edgeR_data)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
    edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
    edgeR_res <- edgeR::glmQLFTest(edgeR_fit, contrast = LT)
  }, check = FALSE, min_iterations = 3
)
