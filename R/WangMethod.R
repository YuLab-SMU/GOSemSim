wangMethod <- function(t1, t2, ont) {
    n1 <- length(t1)
    n2 <- length(t2)

    if (is_supported_go(ont)) {
        rel_df <- ensure_gotbl_cached()
    } else if (is_supported_do(ont)) {
        rel_df <- get_rel_df(ont)
    } else {
        .meshesEnv <- get(".meshesEnv", envir=.GlobalEnv)
        rel_df <- get("meshtbl", envir=.meshesEnv)
    }

    terms <- unique(c(t1, t2))
    sv_list <- lapply(terms, getSV, ont = ont, rel_df = rel_df)
    names(sv_list) <- terms

    calc_sim <- function(ID1, ID2) {
        if (ID1 == ID2) {
            return(1)
        }
        sv.a <- sv_list[[ID1]]
        sv.b <- sv_list[[ID2]]

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

        sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
    }

    res <- matrix(NA_real_, nrow = n1, ncol = n2, dimnames = list(t1, t2))
    if (identical(t1, t2)) {
        for (i in seq_len(n1)) {
            for (j in seq_len(i)) {
                res[i, j] <- calc_sim(t1[i], t2[j])
                res[j, i] <- res[i, j]
            }
        }
        return(res)
    }

    for (j in seq_len(n2)) {
        for (i in seq_len(n1)) {
            res[i, j] <- calc_sim(t1[i], t2[j])
        }
    }
    res
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

    if (is_supported_go(ont)) {
        rel_df <- ensure_gotbl_cached()
    } else if (is_supported_do(ont)) {
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
    get_gosemsim_env()

    res <- yulab.utils::get_cache_element(".GOSemSimEnv", ontbl)
    if (!is.null(res)) return(res)

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

    e <- list()
    e[[ontbl]] <- rel_df
    yulab.utils::update_cache_item(".GOSemSimEnv", e)
    return(rel_df)
}


getSV <- function(ID, ont, rel_df, weight=NULL) {
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

    weight_key <- paste(names(weight), weight, sep = "=", collapse = ";")
    cache_key <- paste(ID, ont, weight_key, sep = "|")
    sv <- yulab.utils::get_cache_element("GOSemSim_SemSimCache", cache_key)
    if (!is.null(sv)) return(sv)

    rel_df <- rel_df[rel_df$Ontology == ont,]
    if (! 'relationship' %in% colnames(rel_df))
        rel_df$relationship <- "other"

    ## GO.db emits relationship names as "isa" / "part of" (see GOBPPARENTS),
    ## but the weight table above is keyed on "is_a" / "part_of". Without
    ## normalising these spellings, every edge fails the match below and is
    ## remapped to "other", so Wang similarity silently uses a uniform 0.7
    ## weight for all edges instead of is_a = 0.8 / part_of = 0.6.
    rel_df$relationship[rel_df$relationship %in% c("isa", "is a")] <- "is_a"
    rel_df$relationship[rel_df$relationship %in% c("part of")]     <- "part_of"
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
    e[[cache_key]] <- sv
    yulab.utils::update_cache_item("GOSemSim_SemSimCache", e)
    
    return(sv)
}

