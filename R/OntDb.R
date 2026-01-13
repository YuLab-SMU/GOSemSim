##' @importClassesFrom AnnotationDbi AnnotationDb
##' @importFrom methods setRefClass
setRefClass("OntDb", contains="AnnotationDb")

#' @importMethodsFrom AnnotationDbi keys
#' @importMethodsFrom AnnotationDbi toTable
setMethod("keys", "OntDb",
    function(x, keytype, ...){
        if(missing(keytype)) keytype <- "id"
        term <- toTable(x)
        term[, keytype]
    }
)

#' @importMethodsFrom AnnotationDbi keytypes
setMethod("keytypes", "OntDb",
    function(x) {
        c("id", "term")
    }

)


#' @importFrom DBI dbReadTable
setMethod("toTable", "OntDb",
    function(x) {
        setNames(dbReadTable(dbconn(x), 'term'), c("id", "term"))
    }
)



#' @importMethodsFrom AnnotationDbi select
#' @importMethodsFrom AnnotationDbi dbconn
setMethod("select", "OntDb",
    function(x, keys, columns, keytype, ...){
        if (missing(keytype)) keytype <- "id"
        keytype <- match.arg(keytype, c("id","term"))
        strKeys <- paste0("\"", keys, "\"", collapse = ",")
        if (keytype == "term") {
            sql_key <- paste("SELECT doid FROM do_term WHERE term in (",
                strKeys, ")")
            doids <- dbQuery(dbconn(x), sql_key)[, 1]
            strKeys <- paste0("\"", doids, "\"", collapse = ",")
        }
        columns <- unique(c("id", columns))

        sqls <- paste("SELECT ", paste(columns, collapse = ","),
            " FROM term")
        columns2 <- setdiff(columns, c("id", "term"))
        for (col in columns2) {
            leftJoin <- paste0("LEFT JOIN  ", col, " USING (id)")
            sqls <- c(sqls, leftJoin)
        }
        sqls <- c(sqls, paste0("WHERE term.id in (", strKeys, ")"))
        sqls <- paste(sqls, collapse = " ")
        res <- dbQuery(dbconn(x), sqls)
        res
    }
)

dbQuery <- getFromNamespace("dbQuery", "AnnotationDbi")

#' @importMethodsFrom AnnotationDbi columns
setMethod("columns", "OntDb",
    function(x) {
        c("id","term", "alias", "synonym", "parent", "children",
            "ancestor", "offspring")
    }
)


get_onto_data <- function(ont = "HDO", output='list', table="offspring") {
    x <- load_onto(ont)
    output <- match.arg(output, c("data.frame", "list"))
    res <- dbReadTable(dbconn(x), table)
    if (output == 'data.frame') return(res)

    # column 1 is ID, column 2 is the related term
    split(res[,2], res[,1]) 
}

#' Load Ontology Database
#' 
#' @param onto character. The ontology to load (e.g., "HDO").
#' @return An `AnnotationDb` object.
#' @importFrom digest digest
#' @importFrom AnnotationDbi loadDb
#' @importFrom yulab.utils download_yulab_file
#' @importFrom yulab.utils user_dir
#' @keywords internal
load_onto <- function(onto = "HDO") {
    .onto <- sprintf(".onto_%s", onto)
    
    db <- yulab.utils::get_cache_element(".GOSemSimEnv", .onto)
    if (!is.null(db)) return(db)

    dbfile <- sprintf("%s.sqlite", onto)
    urls <- c("https://yulab-smu.top/DOSE",
              "https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages")

    # use download_yulab_file from yulab.utils to handle multiple mirrors
    dbfile <- download_yulab_file(dbfile, urls, gzfile = TRUE, appname = "GOSemSim")

    db <- loadDb(dbfile)
    
    yulab.utils::update_cache_item(".GOSemSimEnv", setNames(list(db), .onto))
    return(db)
}

