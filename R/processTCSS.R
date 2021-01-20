#prepare tcss data for TCSS to calculate semantic similarity
#' Title
#'
#' @param ont ontology
#' @param geneAnno the annotation for genes or gene products
#' @param IC information content
#' @param cutoff the topology cutoff
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#'     ## 例子必须要能跑的通，这里goAnno 和 IC都还没赋值，不能直接用。
#'     tcssdata <- process_tcss("BP", goAnno, IC, 3.8)
#' }
process_tcss <- function(ont, geneAnno, IC, cutoff = NULL) {

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

    #mic : max IC value of all terms
    mic <- get_maxIC(meta_terms, IC = IC, mic = NULL)
    #get the corrsponding max ic for each meta-term
    meta_maxIC <- vapply(meta_graph, get_maxIC, IC = IC, mic = mic, numeric(1))

    meta_maxIC["meta"] <- mic

    #build a data.frame

    res <- data.frame(GO = unname(unlist(meta_graph)),
                      clusid = rep(meta_terms,
                                times = vapply(meta_graph, length, integer(1))),
                      stringsAsFactors = FALSE)

    res <- rbind(res, data.frame(GO = meta_terms, clusid = "meta",
                                 stringsAsFactors = FALSE))

    #ICA means altered IC value
    ica <- unname(mapply(
        function(e, f) IC[e] / meta_maxIC[f], res$GO, res$clusid))

    res$ica <- ica

    return(res)
}

#compute ICT (information content topology) for each term
computeICT <- function(GO, offspring) {
    vapply(GO, function(e)
        -log10(length(offspring[[e]]) / length(GO)), numeric(1))
}

#according to the cutoff select meta-terms
get_meta <- function(ont, ict, GO, cutoff) {
    if (is.null(cutoff)) {
        cutoff <- switch(ont,
                         MF = 3.0, 
                         BP = 3.8, 
                         CC = 3.0)
    }
    GO[which(ict <= cutoff)]
}

#for each graph get their max IC value
get_maxIC <- function(terms, IC, mic) {
    all <- IC[terms]
    all <- all[all != Inf & all != -Inf]
    all <- na.omit(all)
    if (length(all) == 0) {
        #assign the mic value
        mic
    }else {
        max(all)
    }
}

#Each meta-term represent a sub-graph-root-node,
#its offspring nodes as the sub-graph-nodes,
#of course, redundancy is needed.
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

#in the meta-terms, if two terms's ict are close(± 20%)
#and satisfied with parent-child relation, then
#the child term is removed.
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
