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

getOffspringIdx <- function(ont, goids) {
    .GOSemSimEnv <- get_gosemsim_env()
    key <- paste0("offspring_idx_", ont, "_", digest::digest(goids))
    if (exists(key, envir = .GOSemSimEnv)) {
        return(get(key, envir = .GOSemSimEnv))
    }
    off <- getOffsprings(ont)
    pos <- stats::setNames(seq_along(goids), goids)
    idx <- lapply(goids, function(id) {
        ids <- off[[id]]
        if (is.null(ids)) integer(0) else as.integer(stats::na.omit(pos[ids]))
    })
    names(idx) <- goids
    assign(key, idx, envir = .GOSemSimEnv)
    idx
}
