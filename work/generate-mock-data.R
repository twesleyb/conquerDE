sce <- scater::mockSCE()

counts <- assay(sce)

meta <- colData(sce)

dge <- scran::convertTo(sce, type="edgeR") # "DESeq2" or "monocle"
