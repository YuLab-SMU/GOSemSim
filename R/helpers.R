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
    key <- paste0("offspring_idx_", ont, "_", digest::digest(goids))
    res <- yulab.utils::get_cache_element("GOSemSim_offspring_idx", key)
    if (!is.null(res)) return(res)
    off <- getOffsprings(ont)
    pos <- stats::setNames(seq_along(goids), goids)
    idx <- lapply(goids, function(id) {
        ids <- off[[id]]
        if (is.null(ids)) integer(0) else as.integer(stats::na.omit(pos[ids]))
    })
    names(idx) <- goids
    e <- list()
    e[[key]] <- idx
    yulab.utils::update_cache_item("GOSemSim_offspring_idx", e)
    idx
}
