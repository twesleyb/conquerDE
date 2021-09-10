#' run_edgeRLRTdeconv
#' 
#' run edgeR LRT analysis using deconvolution normalization \insertRef{Lun2016}{conquerDE}
#'
#' @importFrom scran convertTo
#' @importFrom scuttle computePooledFactors
#' @importFrom edgeR estimateDisp glmFit glmLRT topTags
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export run_edgeRLRTdeconv

run_edgeRLRTdeconv <- function(L, alpha=0.5) {

		# collect inputs
		meta <- L$meta
		counts <- L$counts
		contrast <- L$contrast
		model_matrix <- L$design

		# coerce to sce object
		sce <- SingleCellExperiment::SingleCellExperiment(list("counts"=counts))

		# sizes is the number of cells per pool
		pool_size <- unique(pmin(c(20, 40, 60, 80, 100), ncol(sce)/2))

		# scaling normalization by deconvolving size factors from cell pools
		sce <- scuttle::computePooledFactors(sce, sizes=pool_size, positive=TRUE)

		# convert back to dge
		dge <- scran::convertTo(sce, type = "edgeR")

		# edgeR flow
		dge <- edgeR::estimateDisp(dge, design = model_matrix)
		fit <- edgeR::glmFit(dge, design = model_matrix)
		lrt <- edgeR::glmLRT(fit, contrast=contrast)
		tt <- edgeR::topTags(lrt, n = Inf)

		# collect results
		res <- as.data.frame(tt$table)
		res$candidate <- res$FDR < alpha
		return(res)
}
