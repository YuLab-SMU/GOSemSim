.initial <- function() {
    e <- new.env(parent = emptyenv())
    gotbl <- tryCatch(utils::data("gotbl", package="GOSemSim", envir=e),
                      error = function(e) NULL)
    if (is.null(gotbl)) return(NULL)
    
    yulab.utils::update_cache_item(".GOSemSimEnv", list(gotbl = e$gotbl))
}

ensure_gotbl_cached <- function() {
    gotbl <- yulab.utils::get_cache_element(".GOSemSimEnv", "gotbl")
    if (is.null(gotbl) || !is.data.frame(gotbl)) {
        .initial()
        gotbl <- yulab.utils::get_cache_element(".GOSemSimEnv", "gotbl")
    }
    return(gotbl)
}

get_gosemsim_env <- function() {
    ensure_gotbl_cached()
    invisible(NULL)
}

supported_GO <- function() c("BP", "CC", "MF")
supported_DO <- function() c("DO", "HDO", "HPO", "MPO")

is_supported_go <- function(ont) ont %in% supported_GO()
is_supported_do <- function(ont) ont %in% supported_DO()


## getAncestors(), getParents() and getOffsprings() each materialise a whole
## ontology mapping before the caller picks out a single term with `[[ID]]`.
## For GO that means AnnotationDbi::as.list(GOBPANCESTOR) and friends, which
## costs about half a second; for the downloadable ontologies it is a full
## dbReadTable() of the relation table. The callers typically ask for only one
## or two IDs, so without caching the very same mapping is rebuilt over and
## over again -- a single geneSim() call can rebuild it hundreds of times, and
## it dominates the running time of TCSS and of the IC-based measures.
##
## The mappings are derived from an annotation source that is itself cached
## (GO.db is a static annotation package and load_onto() caches the OntDb), so
## they cannot change within a session and are cached here on the same basis.
##
## Note that the cache key is the ontology as a whole, not the individual term.
## The per-term cache in tcssMethod_internal() does not help here: it saves the
## lookup but every newly seen term still pays for building the whole mapping.
onto_relation <- function(ont, kind, fun) {
    cache_key <- paste(kind, ont, sep = "|")
    res <- yulab.utils::get_cache_element("GOSemSim_ontoRelation", cache_key)
    if (!is.null(res)) {
        return(res)
    }

    res <- fun()

    yulab.utils::update_cache_item("GOSemSim_ontoRelation",
                                   setNames(list(res), cache_key))
    return(res)
}

##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
getAncestors <- function(ont) {
    if (is_supported_go(ont)) {
        return(onto_relation(ont, 'ancestor', function() {
            Ancestors <- switch(ont,
                                MF = GOMFANCESTOR,
                                BP = GOBPANCESTOR,
                                CC = GOCCANCESTOR
                                )
            AnnotationDbi::as.list(Ancestors)
        }))
    }

    onto_relation(ont, 'ancestor',
                  function() get_onto_data(ont, output = 'list', 'ancestor'))
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
getParents <- function(ont) {
    if (is_supported_go(ont)) {
        return(onto_relation(ont, 'parent', function() {
            Parents <- switch(ont,
                            MF = GOMFPARENTS,
                            BP = GOBPPARENTS,
                            CC = GOCCPARENTS
                            )
            AnnotationDbi::as.list(Parents)
        }))
    }

    onto_relation(ont, 'parent',
                  function() get_onto_data(ont, output = 'list', 'parent'))
}

##' @importFrom GO.db GOMFOFFSPRING
##' @importFrom GO.db GOBPOFFSPRING
##' @importFrom GO.db GOCCOFFSPRING
getOffsprings <- function(ont) {
    if (is_supported_go(ont)) {
        return(onto_relation(ont, 'offspring', function() {
            Offsprings <- switch(ont,
                            MF = GOMFOFFSPRING,
                            BP = GOBPOFFSPRING,
                            CC = GOCCOFFSPRING
                            )
            AnnotationDbi::as.list(Offsprings)
        }))
    }

    onto_relation(ont, 'offspring',
                  function() get_onto_data(ont, output = 'list', 'offspring'))
}

##' @importFrom GO.db GOTERM
##' @importFrom AnnotationDbi toTable
prepare_relation_df <- function() {
    gtb <- toTable(GOTERM)
    gtb <- gtb[,c(2:4)]
    gtb <- unique(gtb)
    
    ptb <- lapply(c("BP", "MF", "CC"), function(ont) {
        id <- with(gtb, go_id[Ontology == ont])
        parentMap <- getParents(ont)
        # pid <- AnnotationDbi::mget(id, parentMap)
        pid <- parentMap[id]

        n <- sapply(pid, length)
        cid <- rep(names(pid), times=n)
        relationship <- unlist(lapply(pid, names))
        
        data.frame(id=cid,
                   relationship=relationship,
                   parent=unlist(pid),
                   stringsAsFactors = FALSE)
    }) 
    ptb <- do.call('rbind', ptb)

    gotbl <- merge(gtb, ptb, by.x="go_id", by.y="id")
    save(gotbl, file="gotbl.rda", compress="xz")
    invisible(gotbl)
}

#' @title Get organism name from OrgDb object
#' @param object OrgDb object or OrgDb package name
#' @return Organism name
#' @importFrom yulab.utils load_OrgDb
#' @export
#' @author Guangchuang Yu
get_organism <- function(object) {
    OrgDb <- load_OrgDb(object)
    AnnotationDbi::species(OrgDb)
}

#' @importFrom yulab.utils load_OrgDb
#' @export
yulab.utils::load_OrgDb