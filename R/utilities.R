.initial <- function() {
    pos <- 1
    envir <- as.environment(pos) 
    assign(".GOSemSimEnv", new.env(), envir = envir)
    assign(".SemSimCache", new.env(), envir = envir)
    assign(".ancCache", new.env(), envir = envir)
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    
    tryCatch(utils::data(list="gotbl",
                         package="GOSemSim"))
    gotbl <- get("gotbl")
    assign("gotbl", gotbl, envir = .GOSemSimEnv)
    rm(gotbl, envir = .GlobalEnv)
}

get_gosemsim_env <- function() {
    if (!exists(".GOSemSimEnv")) {
        .initial()
    }
    get(".GOSemSimEnv")
}



##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
getAncestors <- function(ont) {
    if (ont %in% c("MF", "BP", "CC")) {
        Ancestors <- switch(ont,
                            MF = GOMFANCESTOR,
                            BP = GOBPANCESTOR,
                            CC = GOCCANCESTOR
                            )
        anc <- AnnotationDbi::as.list(Ancestors)
        return(anc)
    }

    get_onto_data(ont, output = 'list', 'ancestor') 
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
getParents <- function(ont) {
    if (ont %in% c("MF", "BP", "CC")) {
        Parents <- switch(ont,
                        MF = GOMFPARENTS,
                        BP = GOBPPARENTS,
                        CC = GOCCPARENTS
                        )
        parent <- AnnotationDbi::as.list(Parents)
        return(parent)
    }

    get_onto_data(ont, output = 'list', 'parent') 
}

##' @importFrom GO.db GOMFOFFSPRING
##' @importFrom GO.db GOBPOFFSPRING
##' @importFrom GO.db GOCCOFFSPRING
getOffsprings <- function(ont) {
    if (ont %in% c("MF", "BP", "CC")) {
        Offsprings <- switch(ont,
                        MF = GOMFOFFSPRING,
                        BP = GOBPOFFSPRING,
                        CC = GOCCOFFSPRING
                        )
        offspring <- AnnotationDbi::as.list(Offsprings)
        return(offspring)
    }

    get_onto_data(ont, output = 'list', 'offspring') 
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
