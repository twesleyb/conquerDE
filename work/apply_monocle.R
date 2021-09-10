#' This is the title.
#'
#' This is the description.
#'
#' These are further details.
#'
#' @section A Custom Section:
#'
#' Text accompanying the custom section.
#'
#' @param x A description of the parameter 'x'. The
#'   description can span multiple lines.
#' @param y A description of the parameter 'y'.
#'
#' @export
#'
#' @examples
#' add_numbers(1, 2) ## returns 3

suppressPackageStartupMessages(library(edgeR))

run_monocle <- function(L) {

		L$tpm

    mon <- newCellDataSet(as.matrix(L$tpm), 
                          phenoData = new("AnnotatedDataFrame", 
                                          data = data.frame(condition = L$condt, 
                                                            row.names = colnames(L$tpm))),
                          expressionFamily = tobit())

    monres <- differentialGeneTest(mon, fullModelFormulaStr = " ~ condition")
