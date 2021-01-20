#' Title Method TCSS for semantic similarity measuring
#'
#' @param t1 term vector
#' @param t2 term vector
#' @param semData GOSemSimDATA object
#'
#' @return score vector
#' @noRd
#' @importFrom stats na.omit
#'
#' @examples
#' library(org.Hs.eg.db)
#' semdata <- godata(org.Hs.eg.db, "ENTREZID", "MF", computeIC = TRUE,
#' processTCSS = TRUE, cutoff = NULL)
#' termSim("GO:0000003", "GO:0009987", semdata, method = "TCSS")
#process two term vectors
tcssMethod <- function(t1, t2, semData) {
    matrix( mapply( tcssMethod_internal,
                    rep( t1, length(t2) ),
                    rep( t2, each = length(t1) ),
                    MoreArgs = list( semData = semData ) ),
            dimnames = list( t1, t2 ), ncol = length(t2) ) 
}

#' Title process one term with one term
#'
#' @param ID1 term
#' @param ID2 term
#' @param semData GOSemSimDATA object
#'
#' @return numeric
#' @noRd
#'
tcssMethod_internal <- function(ID1, ID2, semData) {

    tcssdata <- semData@tcssdata
    ont <- semData@ont

    if (length(tcssdata) == 0) {
        stop("tcss data not found, please re-generate your `semData` with `tcssprocess = TRUE`...")
    }

    #get common ancestors
    com_anc <- get_common_anc(ID1 = ID1, ID2 = ID2, ont = ont)

    if (length(com_anc) == 0) return(NA)

    #get cluster-ids for each ID
    clus1_list <- tcssdata[tcssdata[, "GO"] == ID1, "clusid"]
    clus2_list <- tcssdata[tcssdata[, "GO"] == ID2, "clusid"]
    #calculate within different clusters
    value <- lapply(clus1_list, function(e) {
                lapply(clus2_list, get_lca, e,
                       ID1 = ID1, ID2 = ID2,
                       tcssdata = tcssdata, com_anc = com_anc, ont = ont)
    })
    
    value <- na.omit(unlist(value))
    value <- value[value != Inf & value != -Inf]

    if (is.null(value) || length(value) == 0) {
        NULL
    }else {
        #here max value means lowest common ancestor
        return(max(value))
    }
}

#' Title get lowest common ancestors's value 
#'
#' @param clus1 character, cluster-id for ID1
#' @param clus2 character, cluster-id for ID2
#' @param ID1 term
#' @param ID2 term
#' @param tcssdata data.frame
#' @param com_anc character
#' @param ont ontology
#'
#' @return numeric or NULL
#' @noRd
#'
get_lca <- function(clus1, clus2, ID1, ID2, tcssdata, com_anc, ont) {

    #
    if (clus1 == "meta" || clus2 == "meta") {
        return(NULL)
    }
    #if the two clusters are the same one
    if (identical(clus1, clus2)) {
        #get common ancestors inside the cluster
        clus_content <- tcssdata[tcssdata[, "clusid"] == clus1, ]
        com_anc_loc <- match(com_anc, clus_content[, "GO"])
        value <- clus_content[com_anc_loc, "ica"]

    }else {
        #get common ancestors in the meta-cluster
        #clus1 replace ID1, clus2 replace ID2
        clus_content <- tcssdata[tcssdata[, "clusid"] == "meta", ]
        com_anc <- get_common_anc(ID1 = clus1, ID2 = clus2, ont = ont)
        com_anc_loc <- match(com_anc, clus_content[, "GO"])
        value <- clus_content[com_anc_loc, "ica"]
    }

    #value <- na.omit(value[value != Inf & value != -Inf])
    value <- na.omit(value)

    if (is.null(value) || length(value) == 0) return(NULL)
    
    #here max value means one of lowest common ancestors
    max(value)
    
}

#' Title get common ancestors
#'
#' @param ID1 term
#' @param ID2 term
#' @param ont ontology
#'
#' @return character
#' @noRd
#'
get_common_anc <- function(ID1, ID2, ont) {
    ancestor1 <- getAncestors(ont)[[ID1]]
    ancestor2 <- getAncestors(ont)[[ID2]]
    
    if (ID1 == ID2 || ID1 %in% ancestor2) return(ID1)
    
    if (ID2 %in% ancestor1) return(ID2)
    
    setdiff(intersect(ancestor1, ancestor2), "all")
    
}
