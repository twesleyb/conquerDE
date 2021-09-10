#!/usr/bin/env Rscript

# title:
# author:
# description:

# inputs
root <- "~/projects/conquerDE"
devtools::load_all(root)

# additional imports
suppressPackageStartupMessages({
library(dplyr)
library(data.table)
})

# prepare the data
data(meta)
data(filt_counts)

# an RNAseq experiment in four +1 parts:
# * counts matrix
# * experimental meta data (colData)
# * experimental design matrix 
# * contrast vector or matrix 
# * additional row meta data (gene ids, annotations, ect)

# create a design matrix
model_matrix <- model.matrix(~ 0 + condition, data=meta)

# create a contrast vector
contrast <- limma::makeContrasts('conditionIgAN-conditionCTRL', levels=model_matrix)

# save key pieces of data
save(model_matrix, file=file.path(root,"data","model_matrix.rda"),version=2)
save(contrast, file=file.path(root,"data","contrast.rda"),version=2)

L <- list('counts'=filt_counts, 'meta'=meta, 
		  'design'=model_matrix, 'contrast'=contrast)

res <- conquerDE(L,method="edgeRQLF")

