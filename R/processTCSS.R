#' prepare tcss data for TCSS to calculate semantic similarity
#'
#' @param ont ontology
#' @param IC information content
#' @param cutoff the topology cutoff
#'
#' @return list, belonged clusters and its elements for all nodes
#' @noRd
process_tcss <- function(ont, IC, cutoff = NULL) {
    if (length(IC) == 0) {
        stop("IC data not found, please re-generate your `semData` with `computeIC = TRUE`...")
    }

    if (is.null(cutoff)) {
        message("cutoff value is not specified, default value based on human
        data will be taken, or you can call the function 'tcss_cutoff' with your ppidata")
    } else if (cutoff <= 0) {
        stop("cutoff value must be greater than 0")
    }

    GO <- names(IC[!is.infinite(IC)])

    offspring <- switch(ont,
                        MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                        BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                        CC = AnnotationDbi::as.list(GOCCOFFSPRING)
    )

    # calculate ICT
    ICT <- computeICT(GO, offspring = offspring)
    # nodes smaller than cutoff are meta-terms
    meta_terms <- create_meta_terms(ont, ICT = ICT, GO = GO, cutoff = cutoff)
    # if two parent-child nodes' ICT value too close
    meta_terms <- remove_close(meta_terms, ont = ont, ICT = ICT)
    # relationship between cluster-id and its elements
    meta_graph <- create_sub_terms(meta_terms, offspring = offspring)
    #for all GO terms, get the contained elements
    GO_element <- lapply(GO, function(t) {
        meta_graph[meta_graph["cluster"] == t, "element"]
    })
    names(GO_element) <- GO

    # get the max IC value for each graph
    meta_maxIC <- calc_maxIC(meta_terms, GO_element = GO_element, IC = IC)
    
    #return a list, represents term's elements and their value, term's clusters
    res <- lapply(GO, function(t) {
        list("GO" = GO_element[[t]],
             "ica" = unname(IC[GO_element[[t]]] / meta_maxIC[t]),
             "clusid" = meta_graph[meta_graph[, "element"] == t, "cluster"])
    })
    names(res) <- GO
    
    #add "meta" cluster to collect meta_terms
    res$"meta" <- list("GO" = meta_terms,
                       "ica" = IC[meta_terms] / max(meta_maxIC),
                       "clusid" = NULL)

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
    filtered_offspring <- offspring[GO]
    all <- length(GO)
    num <- -log10(1 / all)

    vapply(filtered_offspring, function(off) {
        if (any(is.na(off))) {
            # only term itself
            num
        } else {
            # add term itself
            -log10((sum(off %in% GO) + 1) / all)
        }
    }, numeric(1))
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
                         MF = 3.5,
                         BP = 3.5,
                         CC = 3.2
        )
    }
    GO[ICT <= cutoff]
}

#' calculate every graph's max IC value
#'
#' @param meta_terms character, all cluster ids but "meta"
#' @param GO_element the contained elements
#' @param IC numeric, ICT value
#'
#' @return numeric, max IC in different graphs
#' @noRd
#'
calc_maxIC <- function(meta_terms, GO_element, IC) {
    # mic : max IC value of all terms
    mic <- max(IC[!is.infinite(IC)])

    meta_maxIC <- vapply(meta_terms,
                         function(t) {
                             #all the IC value of elements
                             all <- IC[GO_element[[t]]]
                             all <- all[!is.infinite(all)]
                             # if value is empty, assign the mic value
                             if (length(all) == 0) mic else max(all)
                             }, numeric(1))
    return(meta_maxIC)
}

#' get the relationship between clusters and elements
#'
#' @param meta_terms character, all cluster ids but "meta"
#' @param offspring list
#'
#' @return data.frame, relationship between clusters and elements
#' @noRd
#'
create_sub_terms <- function(meta_terms, offspring) {
    #all terms within cluster
    element <- lapply(meta_terms, function(term) {
        #when no offspring, all shows NA
        all <- offspring[[term]]
        # other sub-root-node
        other_sub <- intersect(all, meta_terms)
        # other sub-root-node's offspring
        other_sub_offs <- unlist(offspring[other_sub])
        # remove
        terms <- setdiff(all, c(other_sub_offs, other_sub))
        # add term itself
        terms <- c(terms, term)
    })
    #add "meta" cluster
    element$"meta" <- meta_terms
    #to unfold into a data.frame
    len <- lapply(element, length)

    res <- data.frame(cluster = rep(c(meta_terms, "meta"), times = len),
                      element = unlist(element))
    return(na.omit(res))
}

#' remove close relation in meta_terms
#'
#' @param meta_terms character, all cluster ids but "meta"
#' @param ont ontology
#' @param ICT numeric, topological information content value
#'
#' @return character, meta_terms with fewer nodes
#' @noRd
#'
remove_close <- function(meta_terms, ont, ICT) {
    parents <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFPARENTS),
                      BP = AnnotationDbi::as.list(GOBPPARENTS),
                      CC = AnnotationDbi::as.list(GOCCPARENTS)
    )
    # reserve all nodes in advance
    all_ <- meta_terms
    for (term1 in meta_terms) {
        # parent term
        obj <- intersect(parents[[term1]], all_)
        for (term2 in obj) {
            if (ICT[term2] != 0 && ICT[term1] / ICT[term2] <= 1.2) {
                # remove when satisfying the condition
                meta_terms <- setdiff(meta_terms, term1)
                break
            }
        }
    }
    # return the left nodes
    return(meta_terms)
}
