#Method TCSS for semantic similarity measuring
#' Title
#'
#' @param ID1 term vector
#' @param ID2 term vector
#' @param semData GOSemSimDATA object
#'
#' @return score vector
#' @export
#'
#' @examples
#' tcssMethod("GO:0000003", "GO:0009987", semdata)
tcssMethod <- function(ID1, ID2, semData) {

    tcssdata <- semData@tcssdata
    ont <- semData@ont

    if (length(tcssdata) == 0) {
        stop("tcss data not found, please re-generate your `semData` with `tcssprocess = TRUE`...")
    }

    #get common ancestors
    com_anc <- get_common_anc(ID1, ID2, ont)

    if (length(com_anc) == 0) return(NA)

    #get cluster-ids for each ID
    clus1_list <- tcssdata[tcssdata$GO == ID1, ]$clusid
    clus2_list <- tcssdata[tcssdata$GO == ID2, ]$clusid
    #calculate within different clusters
    value <- sapply(clus1_list, function(e) {
                sapply(clus2_list, get_lca, e, ID1, ID2,
                       tcssdata, com_anc, ont)
    })

    #value <- value[value != Inf & value != -Inf]

    value <- na.omit(unlist(value))

    if (is.null(value) || length(value) == 0) {
        NULL
    }else {
        #here max value means lowest common ancestor
        return(max(value))
    }
}


#lca : lowest common ancestors
#when term1 belongs to clus1, term2 belongs to clus2
#calculate the semantic similarity value for term1 and term2
get_lca <- function(clus1, clus2, ID1, ID2, tcssdata, com_anc, ont) {

    #
    if (clus1 == "meta" || clus2 == "meta") {
        return(NULL)
    }
    #if the two clusters are the same one
    if (identical(clus1, clus2)) {
        #get common ancestors inside the cluster
        clus_content <- tcssdata[tcssdata$clusid == clus1, ]
        com_anc_loc <- match(com_anc, clus_content$GO)
        value <- clus_content[com_anc_loc, ]$ica

    }else {
        #get common ancestors in the meta-cluster
        #clus1 replace ID1, clus2 replace ID2
        clus_content <- tcssdata[tcssdata$clusid == "meta", ]
        com_anc <- get_common_anc(clus1, clus2, ont)
        com_anc_loc <- match(com_anc, clus_content$GO)
        value <- clus_content[com_anc_loc, ]$ica
    }

    #value <- na.omit(value[value != Inf & value != -Inf])
    value <- na.omit(value)

    if (is.null(value) || length(value) == 0) {
        NULL
    }else {
        #here max value means one of lowest common ancestors
        return(max(value))
    }
}

#get common ancestors
get_common_anc <- function(ID1, ID2, ont) {
    ancestor1 <- getAncestors(ont)[[ID1]]
    ancestor2 <- getAncestors(ont)[[ID2]]
    if (ID1 == ID2) {
        ID1
    } else if (ID1 %in% ancestor2) {
        ID1
    } else if (ID2 %in% ancestor1) {
        ID2
    } else {
        setdiff(intersect(ancestor1, ancestor2), "all")
    }
}

# get_common_anc <- function(ID1, ID2, ont) {
#     #copy from ICMethods.R 
#     if (ont %in% c("MF", "BP", "CC", "DO")) {
#         .anc <- AnnotationDbi::as.list(getAncestors(ont))
#         allid <- union(ID1, ID2)
#         .anc <- .anc[allid]
#         .anc <- .anc[!vapply(.anc, is.empty, logical(1))]
#     } else {
#         mesh_getAnc <- eval(parse(text = "meshes:::getAncestors"))
#         .anc <- lapply(union(ID1, ID2), mesh_getAnc)
#         names(.anc) <- union(ID1, ID2)
#     }
#     #
#     if (ID1 == ID2) {
#         ID1
#     } else if (ID1 %in% .anc[[1]]) {
#         ID1
#     } else if (ID2 %in% .anc[[2]]) {
#         ID2
#     } else {
#         setdiff(intersect(.anc[[1]], .anc[[2]]), "all")
#     }
# }
