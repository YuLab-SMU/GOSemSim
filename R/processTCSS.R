#' prepare tcss data for TCSS to calculate semantic similarity
#'
#' @param ont ontology
#' @param IC information content
#' @param cutoff the topology cutoff
#'
#' @return data.frame, tcssdata, the cluster-id and cluster-value of nodes
#' @noRd
process_tcss <- function(ont, IC, cutoff = NULL) {

    if (length(IC) == 0) {
        stop("IC data not found, please re-generate your `semData` with `computeIC = TRUE`...")
    }

    if (is.null(cutoff)) {
        message("cutoff value is not specified, default value based on human
        data will be taken, or you can call the function 'tcss_cutoff' with your ppidata")
    }

    GO <- names(IC[!is.infinite(IC)])

    offspring <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                      BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                      CC = AnnotationDbi::as.list(GOCCOFFSPRING))

    #calculate ICT
    ICT <- computeICT(GO, offspring = offspring)
    #those nodes smaller than cutoff are meta-terms
    meta_terms <- create_meta_terms(ont, ICT = ICT, GO = GO, cutoff = cutoff)
    #if two parent-child nodes' ICT value too close
    meta_terms <- remove_close(meta_terms, ont = ont, ICT = ICT)
    #get the terms of each sub-graph
    meta_graph <- create_sub_terms(meta_terms, offspring = offspring)
    #get the max IC value for each graph
    meta_maxIC <- calc_maxIC(meta_graph = meta_graph, IC = IC)

    #build a data.frame
    res <- data.frame(GO = unname(unlist(meta_graph)),
                      clusid = rep(meta_terms,
                                times = vapply(meta_graph, length, integer(1))),
                      stringsAsFactors = FALSE)

    res <- rbind(res, data.frame(GO = meta_terms, clusid = "meta",
                                 stringsAsFactors = FALSE))

    #ica means altered IC value
    res$ica <- unname(mapply(
        function(e, f) IC[e] / meta_maxIC[f], res$GO, res$clusid))

    return(res)
}

#' compute ICT (information content topology) for each term
#'
#' @param GO character, all go terms, species specific
#' @param offspring list, offspring nodes
#'
#' @return numeric, ICT value
#' @noRd
#'
computeICT <- function(GO, offspring) {
    all <- length(names(offspring))
    vapply(GO, function(e)
        -log10(length(offspring[[e]]) / all), numeric(1))
}

#' all nodes with ICT value under cutoff are meta_terms
#'
#' @param ont ontology
#' @param ICT numeric, ICT value
#' @param GO character, all go terms, species specific
#' @param cutoff numeric, topological cutoff
#'
#' @return character, sub-graph-root nodes
#' @noRd
#'
create_meta_terms <- function(ont, ICT, GO, cutoff) {
    if (is.null(cutoff)) {
        cutoff <- switch(ont,
                         MF = 3.0,
                         BP = 3.8,
                         CC = 3.0)
    }
    GO[which(ICT <= cutoff)]
}

#' for each graph calculate their max IC value
#'
#' @param meta_graph list, all sub-graphs
#' @param IC numeric, ICT value
#'
#' @return numeric, max IC in different graphs
#' @noRd
#'
calc_maxIC <- function(meta_graph, IC) {
    #mic : max IC value of all terms
    mic <- max(IC[!is.infinite(IC)])

    meta_maxIC <- vapply(meta_graph, function(e) {
        all <- IC[e]
        all <- all[!is.na(all) & !is.infinite(all)]
        #if value is empty, assign the mic value
        if (length(all) == 0) mic else max(all)
        }, numeric(1))
    #"meta" means the graph that contains all meta-terms
    meta_maxIC["meta"] <- mic

    return(meta_maxIC)
}

#' take offspring nodes as sub-graph-nodes
#'
#' @param meta_terms character
#' @param offspring list
#'
#' @return list, all nodes in sub-graphs
#' @noRd
#'
create_sub_terms <- function(meta_terms, offspring) {
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

#' remove close relation in meta_terms
#'
#' @param meta_terms character
#' @param ont ontology
#' @param ICT numeric, topological information content value
#'
#' @return character, meta_terms with less nodes
#' @noRd
#'
remove_close <- function(meta_terms, ont, ICT) {
    parents <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFPARENTS),
                      BP = AnnotationDbi::as.list(GOBPPARENTS),
                      CC = AnnotationDbi::as.list(GOCCPARENTS))
    #reserve all nodes in advance
    all_ <- meta_terms
    for (term1 in meta_terms) {
        #parent term
        obj <- intersect(parents[[term1]], all_)
        for (term2 in obj) {
            if (ICT[term2] != 0 && ICT[term1] / ICT[term2] <= 1.2) {
            #remove when satisfing the condition
            meta_terms <- setdiff(meta_terms, term1)
            break
            }
        }
    }
    #return the left nodes
    return(meta_terms)
}
