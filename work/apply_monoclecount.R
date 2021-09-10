#' run_monoclecount
#'
#' This is the description.
#'
#' @param x A description of the parameter 'x'. The
#' @param y A description of the parameter 'y'.
#'
#' @export run_monoclecount

run_monoclecount <- function(L, alpha=0.05) {

		meta <- L$meta
		counts <- L$counts

      mce <- monocle::newCellDataSet(counts,  phenoData=meta,
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())

      mce <- monocle::estimateSizeFactors(mon)

      mce <- monocle::estimateDispersions(mce)

      res <- monocle::differentialGeneTest(mce, fullModelFormulaStr = " ~ condition")
    
	  return(res)
}
