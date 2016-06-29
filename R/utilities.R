.initial <- function() {
    assign(".GOSemSimEnv", new.env(),.GlobalEnv)
    assign(".SemSimCache", new.env(), .GlobalEnv)
    assign(".ICEnv", new.env(), .GlobalEnv)
    tryCatch(utils::data(list="godata",
                         package="GOSemSim"))
    godata <- get("godata")
    assign("godata", godata, envir = .GOSemSimEnv)

    assign("SupportedSpecies", c("anopheles",
                                  "arabidopsis",
                                  "bovine",
                                  "canine",
                                  "chicken",
                                  "chimp",
                                  "ecolik12",
                                  "ecsakai",
                                  "fly",
                                  "gondii",
                                  "human",
                                  "malaria",
                                  "mouse",
                                  "pig",
                                  "rat",
                                  "rhesus",
                                  "coelicolor",
                                  "celegans",
                                 "xenopus",
                                 "yeast",
                                 "zebrafish"),
           envir=.GOSemSimEnv)
    ## remove "coelicolor" as it is not supported by Bioconductor
}

##' get supported organisms
##'
##'
##' @title getSupported_Org
##' @return supported organisms
##' @export
##' @author Yu Guangchuang
getSupported_Org <- function() {
    if (!exists("GOSemSimEnv")) .initial()
    supported_Org <- get("SupportedSpecies", envir=GOSemSimEnv)
    return(supported_Org)
}


##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
.getParents <- function(ont) {
    Parents <- switch(ont,
                      MF = "GOMFPARENTS",
                      BP = "GOBPPARENTS",
                      CC = "GOCCPARENTS",
                      DO = "DO.db::DOPARENTS"
                      )
    if (ont == "DO") {
        db <- "DO.db"
        requireNamespace(db)
    }
    Parents <- eval(parse(text=Parents))
    return(Parents)
}

prepare_relation_df <- function() {
    gtb <- toTable(GOTERM)
    gtb <- gtb[,c(2:4)]
    gtb <- unique(gtb)

    ptb <- lapply(c("BP", "MF", "CC"), function(ont) {
        id <- with(gtb, go_id[Ontology == ont])
        pid <- mget(id, .getParents(ont))
        
        n <- sapply(pid, length)
        cid <- rep(names(pid), times=n)
        relationship <- lapply(pid, names) %>% unlist
        
        data.frame(id=cid,
                   relationship=relationship,
                   parent=unlist(pid),
                   stringsAsFactors = FALSE)
    }) %>% do.call('rbind', .)
    
    godf <- merge(gtb, ptb, by.x="go_id", by.y="id")
    return(godf)
}
