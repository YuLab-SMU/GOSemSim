#' prepare tcss data for TCSS to calculate semantic similarity
#'
#' @param ont ontology
#' @param IC information content
#' @param cutoff the topology cutoff, users can use tcss_cutoff() function to calculate cutoff value  
#' @param ICT optional pre-computed ICT vector (to avoid recomputation in tcss_cutoff)
#' @param offspring optional pre-fetched offspring list
#'
#' @return list, belonged clusters and its elements for all nodes
#' @noRd
process_tcss <- function(ont, IC, cutoff = NULL, ICT = NULL, offspring = NULL) {
    if (length(IC) == 0 || !any(is.finite(IC))) {
        stop("IC data not found, please re-generate your `semData` with `computeIC = TRUE`...")
    }

    if (is.null(cutoff)) {
        message("As cutoff value is not provided, default value based on human will be used")
        cutoff <- switch(ont,
                         MF = 3.5,
                         BP = 3.5,
                         CC = 3.2,
                         DO = 3.5,
                         MPO = 3.5,
                         HPO = 3.5
                         )
    } else if (cutoff <= 0) {
        stop("cutoff value must be positive")
    }

    GO <- names(IC[!is.infinite(IC)])

    if (is.null(offspring)) {
        offspring <- getOffsprings(ont) 
    }
    if (is.null(ICT)) {
        # calculate ICT
        ICT <- computeICT(GO, offspring = offspring)
    }
    # nodes smaller than cutoff are meta-terms
    meta_terms <- create_meta_terms(ICT = ICT, cutoff = cutoff)
    # if two parent-child nodes' ICT value too close
    meta_terms <- remove_close(meta_terms, ont = ont, ICT = ICT)
    # sort by ICT descending for deterministic cluster assignment
    meta_terms <- meta_terms[order(ICT[meta_terms], decreasing = TRUE)]
    # relationship between cluster-id and its elements
    meta_graph <- create_sub_terms(meta_terms, offspring = offspring)

    # the value each graph is normalised by
    meta_maxIC <- calc_maxIC(meta_graph, IC = IC)
    ica <- lapply(seq_along(meta_graph), function(i) {
        IC[meta_graph[[i]]] / meta_maxIC[i]
    })
    names(ica) <- meta_terms

    aa <- utils::stack(meta_graph)
    bb <- split(as.character(aa$ind), aa$values)
    clusid <- bb[GO]

    res <- list(
        ## meta_terms == names(meta_graph)
        meta_graph = meta_graph,
        ica = ica,
        clusid = clusid
    )

    return(res)
}

#' compute ICT (Topology Information Content) for each term
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

    res <- vapply(filtered_offspring, function(off) {
        if (any(is.na(off))) {
            ## only term itself
            num
        } else {
            ## add term itself
            -log10((sum(off %in% GO) + 1) / all)
        }
    }, numeric(1))
    names(res) <- GO
    return(res)
}


#' all nodes with ICT value under cutoff are meta_terms
#'
#' @param ICT numeric, ICT value with corresponding GO term as name attributes
#' @param cutoff numeric, topological cutoff
#'
#' @return character, sub-graph-root nodes
#' @noRd
create_meta_terms <- function(ICT, cutoff) {
    res <- ICT[ICT <= cutoff]
    names(res)
}

#' the IC value each graph is normalised by
#'
#' @param meta_graph list, the elements of each cluster
#' @param IC numeric, information content
#'
#' @return numeric, the IC each graph is divided by
#' @noRd
#'
calc_maxIC <- function(meta_graph, IC) {
    # mic : max IC value of all terms
    mic <- max(IC[is.finite(IC)])

    meta_maxIC <- vapply(meta_graph,
                         function(t) {
                                        #all the IC value of elements
                             all <- IC[t]
                                        #all <- all[!is.infinite(all)]
                                        # all <- all[!is.infinite(all) & !is.na(all)]
                             all <- all[is.finite(all)]
                                        # A cluster with a single member has no
                                        # internal spread, so `IC / max(IC)` is
                                        # identically 1 no matter how general or
                                        # how specific that term is. Dividing by
                                        # `mic` instead scores such a term on the
                                        # same scale as every other term.
                                        #
                                        # This is not a rare corner case: at the
                                        # default BP cutoff (3.5) GO:0008150, the
                                        # root of the ontology, is the only
                                        # member of its own cluster. Every pair
                                        # of BP terms has the root as a common
                                        # ancestor, so an ICA of 1 for the root
                                        # made `max(sim_value)` in
                                        # tcssMethod_internal() return 1 for
                                        # essentially every pair. Other cutoffs
                                        # and other ontologies produce further
                                        # singletons (GO:0050896 and GO:0046337
                                        # for BP, GO:0052745 for MF), which the
                                        # same collapse would inflate to 1.
                             if (length(all) <= 1) mic else max(all)
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
    ## #add "meta" cluster
    ## element$"meta" <- meta_terms
    ## #to unfold into a data.frame
    ## len <- lapply(element, length)

    ## res <- data.frame(cluster = rep(c(meta_terms, "meta"), times = len),
    ##                   element = unlist(element))
    ## return(na.omit(res))

    names(element) <- meta_terms
    return(element)
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
    parents <- getParents(ont) 

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

