##' parse GAF files
##'
##' given a GAF file, this function extracts the information from it
##' @title read.gaf
##' @rdname read-gaf
##' @param file GAF file
##' @param asis logical, whether output the original contains of the file and only works if 'add_indirect_GO = FALSE'
##' @param add_indirect_GO whether to add indirect GO annotation 
##' @return A data.frame. Original table if 'asis' works, otherwise contains 3 conlumns of 'GENE', 'GO' and 'ONTOLOGY'
##' @export
read.gaf <- function(file, asis = FALSE, add_indirect_GO = FALSE) {
  GafFile <- read.gaf2(file)
  if (!add_indirect_GO && asis) return(GafFile)

  new.data.frame <- GafFile[, c("DB_Object_ID", "GOID", "Aspect")]
  names(new.data.frame) <- c("GENE", "GO", "ONTOLOGY")
  ont <- setNames(c("MF", "CC", "BP"), c("F", "C", "P"))
  new.data.frame$ONTOLOGY <- ont[new.data.frame$ONTOLOGY]

  if (!add_indirect_GO) return(new.data.frame)

  ## use buildGOmap function to append indirect annotation
  buildGOmap(new.data.frame)
}

##' @importFrom GO.db GO.db
##' @importFrom AnnotationDbi columns
goid2term <- function(simplify = TRUE) {
  go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
  if (simplify) {
    go.ALL <- go.ALL[, c("GOID", "TERM")]
  }

  return(go.ALL)
}

##' @rdname read-gaf
##' @export
parse_gff <- read.gaf

# only read the file with selected columns
##' @importFrom utils read.delim
read.gaf2 <- function(GafFile, nrows = -1) {
  cat("Reading ", GafFile, ": ", sep = "")
  GafFile <-
    read.delim(
      GafFile,
      sep = "\t",
      as.is = TRUE,
      quote = "\"",
      fill = TRUE,
      header = FALSE,
      nrows = nrows,
      comment.char = "!"
    )
  GafFile <- GafFile[, c(2, 3, 5, 7, 9, 10)]
  names(GafFile) <- c(
    "DB_Object_ID",
    "DB_Object_Symbol",
    "GOID",
    "Evidence_Code",
    "Aspect",
    "DB_Object_Name"
  )
  cat("found", nrow(GafFile), "rows in this GAF file")
  return(GafFile)
}
