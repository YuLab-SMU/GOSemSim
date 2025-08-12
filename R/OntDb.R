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
        dbReadTable(dbconn(x), 'term') |>
        setNames(c("id", "term"))
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

#' @importFrom digest digest
#' @importFrom AnnotationDbi loadDb
#' @importFrom R.utils gunzip
load_onto <- function(onto = "HDO") {
    .env <- get_gosemsim_env()
    .onto <- sprintf(".onto_%s", onto)
    if (exists(.onto, envir=.env)) {
        db <- get(.onto, envir=.env)
        return(db)
    }

    dir <- rappdirs::user_data_dir("GOSemSim", appauthor=NULL)

    if (!dir.exists(dir)) dir.create(dir)

    dbfile0 <- sprintf("%s.sqlite", onto)
    dbfile <- file.path(dir, dbfile0)

    if (file.exists(dbfile)) {
        base_url <- 'https://yulab-smu.top/DOSE'
        md5_url <- sprintf("%s/md5.txt", base_url)
        md5 <- tryCatch(read.delim(md5_url, header=FALSE), 
            error = function(e) NULL)
        if (is.null(md5)) {
            base_url <- 'https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages'
            md5_url <- sprintf("%s/md5.txt", base_url)
            md5 <- read.delim(md5_url, header=FALSE)
        }
        md5_remote <- md5[md5[,1] == dbfile0, 2]
        md5_local <- digest::digest(dbfile, algo='md5', file=TRUE)
        if (md5_remote != md5_local) {
            msg <- sprintf("%s is outdated, download the latest version...\n", dbfile0)
            cat(msg)
            need_dl <- TRUE
        } else {
            need_dl <- FALSE
        }
    } else {
        msg <- sprintf("%s is not found, download it online...\n", dbfile0)
        cat(msg)
        need_dl <- TRUE
    }

    if (need_dl) {
        url <- sprintf('%s/%s.gz', base_url, dbfile0)
        gzdbfile <- sprintf("%s.gz", dbfile)
        yulab.utils:::mydownload(url, gzdbfile)
        R.utils::gunzip(gzdbfile, overwrite = TRUE)
    } 

    db <- loadDb(dbfile)
    assign(.onto, db, envir = .env)
    return(db)
}

