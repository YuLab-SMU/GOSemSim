wangMethod <- function(t1, t2, ont) {
    n1 <- length(t1)
    n2 <- length(t2)
    if (n2 > 256) {
        res <- matrix(NA_real_, nrow = n1, ncol = n2, dimnames = list(t1, t2))
        step <- 256L
        for (start in seq.int(1L, n2, by = step)) {
            end <- min(start + step - 1L, n2)
            block <- t2[start:end]
            cols <- lapply(block, function(b) vapply(t1, function(a) wangMethod_internal(a, b, ont = ont), numeric(1)))
            res[, start:end] <- do.call(cbind, cols)
        }
        return(res)
    }
    matrix(mapply(wangMethod_internal,
                  rep(t1, n2),
                  rep(t2, each = n1),
                  MoreArgs = list(ont = ont)),
           dimnames = list(t1, t2), ncol = n2)
}


#' Method Wang for semantic similarity measuring
#'
#' @title wangMethod
#' @param ID1 Ontology Term
#' @param ID2 Ontology Term
#' @param ont Ontology
#' @return semantic similarity score
#' @export
#' @author Guangchuang Yu <https://yulab-smu.top>
wangMethod_internal <- function(ID1, ID2, ont="BP") {
    if (ID1 == ID2)
        return (sim=1)

    if (ont %in% c("BP", "CC", "MF")) {
        .GOSemSimEnv <- get_gosemsim_env()
        rel_df <- get("gotbl", envir=.GOSemSimEnv)
    } else if (ont %in% c("HDO", "HPO", "MPO")) {
        rel_df <- get_rel_df(ont)
    } else {
        .meshesEnv <- get(".meshesEnv", envir=.GlobalEnv)
        rel_df <- get("meshtbl", envir=.meshesEnv)
    } 
    
    
    sv.a <- getSV(ID1, ont, rel_df)
    sv.b <- getSV(ID2, ont, rel_df)

    if(all(is.na(sv.a)) || all(is.na(sv.b)))
        return (NA)

    idx         <- intersect(names(sv.a), names(sv.b))
    inter.sva   <- sv.a[idx]
    inter.svb   <- sv.b[idx]
    if (is.null(inter.sva) ||
        is.null(inter.svb) ||
        length(inter.sva) == 0 ||
        length(inter.svb) ==0) {
        return (NA)
    } 
    
    sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
    return(sim)
}

get_rel_df <- function(ont) {
    ontbl <- sprintf("%stbl", ont)
    .GOSemSimEnv <- get_gosemsim_env()

    if (exists(ontbl, envir=.GOSemSimEnv)) {
        res <- get(ontbl, envir=.GOSemSimEnv)
        return(res)
    }

    ont_db <- load_onto(ont)
    gtb <- toTable(ont_db)
    gtb <- gtb[,1, drop=FALSE]
    gtb <- unique(gtb)

    id <- gtb$id
    parent <- getParents(ont)
    pid <- parent[id]
    cid <- rep(names(pid), times=sapply(pid, length))

    ptb <- data.frame(id=cid,
                      relationship = 'other',
                      parent = unlist(pid),
                      Ontology = ont,
                      stringsAsFactors = FALSE)

    rel_df <- merge(gtb, ptb, by="id")
    rel_df <- rel_df[!is.na(rel_df$id), ]
    rel_df <- rel_df[!is.na(rel_df$parent), ]

    assign(ontbl, rel_df, envir = .GOSemSimEnv)
    return(rel_df)
}


getSV <- function(ID, ont, rel_df, weight=NULL) {
    sv <- yulab.utils::get_cache_element("GOSemSim_SemSimCache", ID)
    if (!is.null(sv)) return(sv)

    if (ont == "HDO") {
        topNode <- "DOID:4"
    } else if (ont == "MPO") {
       topNode <- "MP:0000001"
    } else {
        topNode <- "all"
    }
    
    if (ID == topNode) {
        sv <- 1
        names(sv) <- topNode
        return (sv)
    }
    
    if (is.null(weight)) {
        weight <- c(0.8, 0.6, 0.7)
        names(weight) <- c("is_a", "part_of", "other")
    }

    rel_df <- rel_df[rel_df$Ontology == ont,]
    if (! 'relationship' %in% colnames(rel_df))
        rel_df$relationship <- "other"
    
    rel_df$relationship[!rel_df$relationship %in% c("is_a", "part_of")] <- "other"
    

    sv <- 1
    names(sv) <- ID
    allid <- ID

    idx <- which(rel_df[,1] %in% ID)
    while (length(idx) != 0) {
        p <- rel_df[idx,]
        pid <- p$parent
        allid <- c(allid, pid)
        
        sv <- c(sv, weight[p$relationship]*sv[p[,1]])
        names(sv) <- allid
        idx <- which(rel_df[,1] %in% pid)
    }

    sv <- sv[!is.na(names(sv))]
    sv <- sv[!duplicated(names(sv))]

    if (!(ont %in% c("DO", "MPO")))
        sv[topNode] <- 0

    e <- list()
    e[[ID]] <- sv
    yulab.utils::update_cache_item("GOSemSim_SemSimCache", e)
    
    return(sv)
}

