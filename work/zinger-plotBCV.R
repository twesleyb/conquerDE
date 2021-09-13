#!/usr/bin/env Rscript

# title: *.R
# author: twab
# description: plotting BCV before and after applying zinbwave weights

# inputs
alpha = 0.05

# load project
root <- "~/projects/conquerDE"
devtools::load_all(root)

# load data
data(meta)
data(counts) # unfiltered counts
data(contrast)
data(model_matrix)

# rm rows (genes) with zero counts
drop <- rowSums(counts) == 0
warning("Removing ", sum(drop), " rows with zero counts.")
counts <- counts[!drop, ]

# fit zinbwave model
nthreads <- BiocParallel::multicoreWorkers() - 1
zinb <- zinbwave::zinbFit(counts,
						  X = model_matrix, epsilon = 1e12,
						  BPPARAM = BiocParallel::MulticoreParam(nthreads)
						  )

# collect weights
weights <- zinbwave::computeObservationalWeights(zinb, counts)

# edgeR flow
dge <- edgeR::DGEList(counts, samples=meta)
dge <- edgeR::calcNormFactors(dge)

# KVdB does it like this
dge$weights <- weights
dge <- edgeR::estimateGLMTagwiseDisp(edgeR::estimateGLMCommonDisp(dge, model_matrix), model_matrix, prior.df=0)

edgeR::plotBCV(dge)

dge <- edgeR::estimateDisp(dge, design = model_matrix, prior.df=0)
edgeR::plotBCV(dge)

# add weights to edgeR object
dge <- dge
dge$weights <- weights

# resume edgeR flow
dge <- edgeR::estimateDisp(wdge, design = model_matrix)

#edgeR::plotBCV(dge)

# statistical testing
fit <- edgeR::glmFit(dge, design = model_matrix)

# zinbwave LRT with weights 
# NOTE: why not edgeR glmLRT here?
lrt <- zinbwave::glmWeightedF(fit, contrast = contrast, independentFiltering = TRUE)

# collect results
res <- lrt$table
res$candidate <- res$padjFilter < alpha

hist(res$PValue)
