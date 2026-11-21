#' @param semData GOSemSimDATA object
#' @param measure One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
#' @param combine One of "max", "avg", "rcmax", "BMA" methods, used to combine multiple term scores. "BMA" is the best-match average: it averages the maxima of the rows and of the columns of the score matrix, so that each term of a set contributes the score of its best match rather than a single overall maximum. With a single term on one side this is `(max + sum) / (n + 1)` over the `n` scores in that row or column.
