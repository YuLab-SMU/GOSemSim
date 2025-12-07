#' @keywords internal
buildGoMatrix <- function(uniqueGO, semData, measure) {
    mgoSim(uniqueGO, uniqueGO, semData, measure = measure, combine = NULL)
}

#' @keywords internal
uniqueGOFrom <- function(golist) {
    unique(unlist(golist))
}

#' @keywords internal
subsetCombine <- function(go_matrix, gos1, gos2, combine) {
    combineScores(go_matrix[gos1, gos2, drop = FALSE], combine = combine)
}
