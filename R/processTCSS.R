#' Title prepare tcss data for TCSS to calculate semantic similarity
#'
#' @param ont ontology
#' @param IC information content
#' @param cutoff the topology cutoff
#'
#' @return data.frame
#' @noRd
process_tcss <- function(ont, IC, cutoff = NULL) {

    if (length(IC) == 0) {
        stop("IC data not found, please re-generate your `semData` with `computeIC = TRUE`...")
    }

    if (is.null(cutoff)) {
        message("cutoff value is not specified, default value based on human
        data will be taken, or you can call the function 'get_cutoff' with your testdata")
    }

    GO <- unique(names(IC))

    offspring <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                      BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                      CC = AnnotationDbi::as.list(GOCCOFFSPRING))

    #calculate ict
    ict <- computeICT(GO, offspring = offspring)
    #those nodes smaller than cutoff are meta-terms
    meta_terms <- get_meta(ont = ont, ict, GO = GO, cutoff = cutoff)
    #if two parent-child's ict value too close
    meta_terms <- remove_close(meta_terms, ont = ont, ict = ict)
    #get the terms of each sub-graph
    meta_graph <- get_sub_terms(meta_terms, offspring = offspring)
    #get the max IC value for each graph
    meta_maxIC <- get_maxIC(meta_graph = meta_graph, IC = IC)

    #build a data.frame
    res <- data.frame(GO = unname(unlist(meta_graph)),
                      clusid = rep(meta_terms,
                                times = vapply(meta_graph, length, integer(1))),
                      stringsAsFactors = FALSE)

    res <- rbind(res, data.frame(GO = meta_terms, clusid = "meta",
                                 stringsAsFactors = FALSE))

    #ICA means altered IC value
    res$ica <- unname(mapply(
        function(e, f) IC[e] / meta_maxIC[f], res$GO, res$clusid))

    return(res)
}

#' Title compute ICT (information content topology) for each term
#'
#' @param GO character
#' @param offspring list
#'
#' @return numeric
#' @noRd
#'
computeICT <- function(GO, offspring) {
    vapply(GO, function(e)
        -log10(length(offspring[[e]]) / length(GO)), numeric(1))
}

#' Title according to the cutoff select meta-terms
#'
#' @param ont ontology
#' @param ict numeric
#' @param GO character
#' @param cutoff numeric
#'
#' @return character
#' @noRd
#'
get_meta <- function(ont, ict, GO, cutoff) {
    if (is.null(cutoff)) {
        cutoff <- switch(ont,
                         MF = 3.0,
                         BP = 3.8,
                         CC = 3.0)
    }
    GO[which(ict <= cutoff)]
}

#' Title for each graph get their max IC value
#'
#' @param meta_graph list, all sub-graphs
#' @param IC numeric
#'
#' @return numeric
#' @noRd
#'
get_maxIC <- function(meta_graph, IC) {
    #mic : max IC value of all terms
    mic <- max(IC[IC != Inf & IC != -Inf])

    meta_maxIC <- vapply(meta_graph, function(e) {
        all <- IC[e]
        all <- na.omit(all[all != Inf & all != -Inf])
        #if value is empty, assign the mic value
        if (length(all) == 0) mic else max(all)
        }, numeric(1))
    #"meta" means the graph contains all meta-terms
    meta_maxIC["meta"] <- mic

    return(meta_maxIC)
}

#' Title take offspring node as sub-graph-nodes
#'
#' @param meta_terms character
#' @param offspring list
#'
#' @return list
#' @noRd
#'
get_sub_terms <- function(meta_terms, offspring) {
    res <- lapply(meta_terms, function(e) {
        all <- offspring[[e]]
        #other sub-root-node
        other_sub <- intersect(all, meta_terms)
        #other sub-root-node's offspring
        other_sub_offs <- unlist(offspring[other_sub])
        #remove
        terms <- setdiff(all, other_sub_offs)
        #add term itself
        terms <- c(terms, e)
    })
    names(res) <- meta_terms
    return(res)
}

#' Title remove close-relation in meta_terms
#'
#' @param meta_terms character
#' @param ont ontology
#' @param ict numeric
#'
#' @return character
#' @noRd
#'
remove_close <- function(meta_terms, ont, ict) {
    parents <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFPARENTS),
                      BP = AnnotationDbi::as.list(GOBPPARENTS),
                      CC = AnnotationDbi::as.list(GOCCPARENTS))
    #pre-reserve all nodes in advance
    all_ <- meta_terms
    for (term1 in meta_terms) {
        #parent term
        obj <- intersect(parents[[term1]], all_)
        for (term2 in obj) {
            if (ict[term2] != 0 & ict[term1] / ict[term2] <= 1.2) {
            #remove when satisfing the condition
            meta_terms <- setdiff(meta_terms, term1)
            break
            }
        }
    }
    #return the left nodes
    return(meta_terms)
}
