#!/usr/bin/env Rscript

# title:
# author: twab
# description: basic zinbwave flow

root <- "~/projects/conquerDE"
devtools::load_all(root)

data(filt_counts)
data(meta)
data(model_matrix)
data(contrast)

library(zinbwave)
#library(scRNAseq)
#library(matrixStats)
#library(magrittr)
#library(ggplot2)
#library(biomaRt)

#BiocManager::install("zinbwave")

# load a dataset
sce <- scRNAseq::ReprocessedFluidigmData(assays = "tophat_counts")

# create SingleCellExperiment class
#counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(list(counts=filt_counts),colData=meta)

# filtering
#keep <- rowSums(assay(sce)>5)>5
#table(keep)

#sce <- sce[keep,]

# fit zinbwave model
fit <- zinbwave::zinbFit(sce, K=2, epsilon=1e12)

# DE testing; X = sample design matrix, V = gene design matrix
res <- zinbwave::zinbwave(sce, fitted_model = fit, K = 2, epsilon=1e12, observationalWeights = TRUE)
