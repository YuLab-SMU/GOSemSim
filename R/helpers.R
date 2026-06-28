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

pairwiseCombineMatrix <- function(labels, golist, go_matrix, combine, BPPARAM = NULL) {
    n <- length(golist)
    scores <- matrix(NA_real_, nrow = n, ncol = n)
    if (!is.null(labels)) {
        rownames(scores) <- labels
        colnames(scores) <- labels
    }
    pair_index <- as.list(seq_len(n))
    if (n > 1) {
        pair_index <- c(pair_index, utils::combn(seq_len(n), 2, simplify = FALSE))
    }

    compute_pair <- function(idx) {
        if (length(idx) == 1) {
            i <- idx
            j <- idx
        } else {
            i <- idx[1]
            j <- idx[2]
        }

        gos1 <- golist[[i]]
        gos2 <- golist[[j]]
        if (length(gos1) == 0 || length(gos2) == 0) {
            return(list(i = i, j = j, score = NA_real_))
        }
        list(i = i, j = j, score = subsetCombine(go_matrix, gos1, gos2, combine))
    }

    if (is.null(BPPARAM)) {
        res <- lapply(pair_index, compute_pair)
    } else {
        rlang::check_installed("BiocParallel", "for parallel pairwise similarity calculation")
        res <- BiocParallel::bplapply(pair_index, compute_pair, BPPARAM = BPPARAM)
    }

    for (x in res) {
        scores[x$i, x$j] <- x$score
        scores[x$j, x$i] <- x$score
    }

    scores
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
