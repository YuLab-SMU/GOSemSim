#' Method TCSS for semantic similarity measuring
#'
#' @param t1 term vector
#' @param t2 term vector
#' @param semData GOSemSimDATA object
#'
#' @return vector, similarity score for t1 and t2
#' @noRd
#'
#' @examples
#' library(org.Hs.eg.db)
#' semdata <- godata(org.Hs.eg.db,
#'   keytype = "ENTREZID", ont = "MF",
#'   computeIC = TRUE, processTCSS = TRUE, cutoff = NULL
#' )
#' termSim("GO:0000003", "GO:0009987", semdata, method = "TCSS")
tcssMethod <- function(t1, t2, semData) {
    matrix(mapply(tcssMethod_internal,
                  rep(t1, length(t2)),
                  rep(t2, each = length(t1)),
                  MoreArgs = list(semData = semData)),
    dimnames = list(t1, t2), ncol = length(t2)
    )
}

#' process one term with one term
#'
#' @param ID1 term
#' @param ID2 term
#' @param semData GOSemSimDATA object
#'
#' @return numeric, similarity score for ID1 and ID2
#' @noRd
#'
tcssMethod_internal <- function(ID1, ID2, semData) {
    tcssdata <- semData@tcssdata
    ont <- semData@ont

    if (length(tcssdata) == 0) {
        stop("tcssdata not found, please re-generate your `semData` with `tcssprocess = TRUE`...")
    }

    # get common ancestors
    com_anc <- ancestors_in_common(ID1 = ID1, ID2 = ID2, ont = ont)

    if (length(com_anc) == 0) {
        return(NA)
    }

    # belonged cluster-ids for each ID
    clus1_list <- tcssdata[[ID1]][["clusid"]]
    clus2_list <- tcssdata[[ID2]][["clusid"]]

    # calculate within different clusters
    sim_value <- unlist(mapply(calc_lca,
                               rep(clus1_list, length(clus2_list)),
                               rep(clus2_list, each = length(clus1_list)),
                               MoreArgs = list(
                                   ID1 = ID1, ID2 = ID2,
                                   tcssdata = tcssdata,
                                   com_anc = com_anc, ont = ont
                               )))

    if (is.null(sim_value) || length(sim_value) == 0) {
        return(NA)
    }
    # here max value means lowest common ancestor
    max(sim_value)
}

#' calculate lowest common ancestors' value
#'
#' @param clus1 character, cluster-id for ID1
#' @param clus2 character, cluster-id for ID2
#' @param ID1 term
#' @param ID2 term
#' @param tcssdata list, belonged clusters and its elements for all nodes
#' @param com_anc character, common ancestors
#' @param ont ontology
#'
#' @return numeric/NULL, similarity value for ID1 with clus1 and ID2 with clus2
#' @noRd
#'
calc_lca <- function(clus1, clus2, ID1, ID2, tcssdata, com_anc, ont) {
    # "meta" only represents meta-graph, do not have relations
    if (clus1 == "meta" || clus2 == "meta") {
        return(NULL)
    }
    # if the two clusters are the same one
    if (identical(clus1, clus2)) {
        # all cluster-nodes inside cluster
        clus_content <- tcssdata[[clus1]]
    } else {
        # all cluster-nodes inside "meta" cluster
        # common ancestors are from clus1 and clus2
        clus_content <- tcssdata[["meta"]]
        com_anc <- ancestors_in_common(ID1 = clus1, ID2 = clus2, ont = ont)
    }

    # get common ancestors' position, nomatch helps not to introduce NA
    com_anc_pos <- match(com_anc, clus_content[["GO"]], nomatch = 0)
    # common ancestors' ica value is all the possible sim value
    sim_value <- clus_content[["ica"]][com_anc_pos]

    if (is.null(sim_value) || length(sim_value) == 0) {
        return(NULL)
    }
    # here max similarity value means one of lowest common ancestors
    max(sim_value)
}

#' collect common ancestors
#'
#' @param ID1 term
#' @param ID2 term
#' @param ont ontology
#'
#' @return character, common ancestors for ID1 and ID2
#' @noRd
#'
ancestors_in_common <- function(ID1, ID2, ont) {
    ancestor1 <- ancestors_envir(ID1, ont)
    ancestor2 <- ancestors_envir(ID2, ont)

    if (ID1 == ID2 || ID1 %in% ancestor2) {
        return(ID1)
    }

    if (ID2 %in% ancestor1) {
        return(ID2)
    }

    setdiff(intersect(ancestor1, ancestor2), "all")
}

#' get ancestors from environment
#'
#' @param ID term
#' @param ont ontology
#'
#' @return ancestors for ID
#' @noRd
ancestors_envir <- function(ID, ont) {
    if (!exists(".ancCache")) .initial()
    .ancCache <- get(".ancCache", envir = .GlobalEnv)

    if (exists(ID, envir = .ancCache)) {
        return(get(ID, envir = .ancCache))
    }
    ancestors <- getAncestors(ont)[[ID]]
    assign(ID, ancestors, envir = .ancCache)
    return(ancestors)
}
