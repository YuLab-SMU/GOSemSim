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

##' load OrgDb
##'
##' 
##' @title load_OrgDb
##' @param OrgDb OrgDb object or OrgDb name
##' @return OrgDb object
##' @importFrom methods is
##' @importFrom utils getFromNamespace 
##' @export
##' @author Guangchuang Yu
load_OrgDb <- function(OrgDb) {
    #if (is(OrgDb, "character")) {
    #    require(OrgDb, character.only = TRUE)
    #    OrgDb <- eval(parse(text=OrgDb))
    #}
    if (is(OrgDb, "character")) {
        OrgDb <- utils::getFromNamespace(OrgDb, OrgDb)
    } 
    
    return(OrgDb)
}

##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
getAncestors <- function(ont) {
    Ancestors <- switch(ont,
                        MF = "GOMFANCESTOR",
                        BP = "GOBPANCESTOR",
                        CC = "GOCCANCESTOR",
                        DO = "HDO.db::HDOANCESTOR",
                        MPO = "MPO.db::MPOANCESTOR"
                        )
    if (ont == "DO") {
        db <- "HDO.db"
        ## require(db, character.only=TRUE)
        requireNamespace(db)
    }

    if (ont == "MPO") {
        db <- "MPO.db"
        requireNamespace(db)
    }
    return (eval(parse(text=Ancestors)))
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
getParents <- function(ont) {
    Parents <- switch(ont,
                      MF = "GOMFPARENTS",
                      BP = "GOBPPARENTS",
                      CC = "GOCCPARENTS",
                      DO = "HDO.db::HDOPARENTS"
                      )
    if (ont == "DO") {
        db <- "HDO.db"
        requireNamespace(db)
    }
    Parents <- eval(parse(text=Parents))
    return(Parents)
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
        pid <- AnnotationDbi::mget(id, parentMap)
        
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
