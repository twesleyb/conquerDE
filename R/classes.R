#' conquerList
#'
#' create a S4 class
#'
#' @slot counts A counts matrix.
#' @slot meta A table with sample meta data (colData), coerced to data.frame.
#' @slot contrast A matrix specifying contrast
#' @slot model_matrix A matrix specifying experimental design
#'
#' @export conquerList

# FIXME: how to implement? 

conquerClass <- setClass("conquerList", 
		 slots = c(counts="matrix", meta="data.frame", model_matrix="matrix",contrast="matrix")
		 )

conquerList <- function(counts, meta, model_matrix, contrast) {
		L <- conquerClass(counts=counts, meta=as.data.frame(meta), model_matrix=model_matrix, contrast=contrast)
}
