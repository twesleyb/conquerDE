#!/usr/bin/env Rscript

# title: *.R
# author: twab
# description:
# Usage of the fgsea package.
# BiocManager::install("fgsea")

library(fgsea)

# Load example pathways.
data(examplePathways)

# pathways are stored in a list
class(examplePathways)

str(examplePathways)

# the example list has 1,457 "pathways"
length(examplePathways)

# this is the first pathway, which contains
# all the Entrez gene IDs of the Meiotic_Synapsis pathway
examplePathways[1]

# Load example genes.
data(exampleRanks)
head(exampleRanks) # names are entrez ids.

# What happens if we reverse them?
exampleRanks <- exampleRanks[rev(order(exampleRanks))]

# Result seems to be the same!

# What happens if we randomize them.
#exampleRanks <- exampleRanks[sample(length(exampleRanks))]
# Result seems to be the same!

# GSEA Analysis
results <- fgsea(
  pathways = examplePathways,
  stats = exampleRanks, minSize = 15,
  maxSize = 500,
  nperm = 100000
)
results <- results[order(results$pval), ]

head(results)

## ---- work

root <- "~/projects/conquerDE"
devtools::load_all(root)

# need human pathways!

data(zinger_res)

zinger_res$mus_entrez <- geneLists::getHomologs(zinger_res$entrez,species="human")
gene_ranks <- zinger_res$mus_entrez
gene_ranks <- setNames(zinger_res$logFC * log(zinger_res$PValue),nm=zinger_res$entrez)
gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

idx = names(gene_ranks) %in% zinger_res$mus_entrez[zinger_res$candidate]
gene_ranks = gene_ranks[idx]

head(gene_ranks)

results <- fgseaMultilevel(
  pathways = examplePathways,
  stats = gene_ranks, minSize = 15,
  maxSize = 500,
)

results <- results[order(results$pval), ]

head(results)

