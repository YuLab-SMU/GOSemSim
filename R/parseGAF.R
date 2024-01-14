##' parse GAF files
##'
##' given a GAF file, this function extracts the information from it and add indirect GO annotation
##' @title read.gaf
##' @rdname read-gaf
##' @param GafFile GAF file
##' @return a list with two dataframes
##' @export
##' @importFrom GO.db GO.db
##' @importFrom utils read.delim
##' @importFrom AnnotationDbi columns
read.gaf <- function(GafFile) {
  GafFile <- read.gaf2(GafFile)
  new.data.frame <- GafFile[, c("GOID", "DB_Object_ID")]
  
  ##use buildGOmap function to get information related
  build.df <- buildGOmap(new.data.frame)
  
  ##rename this colnames to be same with the buildGOmap function to facilatate bind information.
  names(new.data.frame) <- c("GO", "Gene")
  
  ##bind the information needed
  bind.info <- rbind(build.df, new.data.frame)
  
  ##process  data
  bind.info <- bind.info[order(bind.info$GO, bind.info$Gene), ]
  bind.info <- bind.info[!duplicated(bind.info), ]
  bind.info[, "Gene"] <- as.character(bind.info$Gene)
  
  go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
  need.anno <- go.ALL[go.ALL$GOID %in% unique(bind.info$GO), ]
  
  # get Ontology information
  if ("Ontology" %in% colnames(go.ALL)) {
      bind.info$Ontology <- need.anno$Ontology
      bind.info$Ontology <- gsub("molecular_function", "MF", bind.info$Ontology)
      bind.info$Ontology <- gsub("cellular_component", "CC", bind.info$Ontology)
      bind.info$Ontology <- gsub("biological_process", "BP", bind.info$Ontology)
  } 

  # detach extra output
  extra.otp <- need.anno[, c("GOID", "TERM")]

  # changed output
  #list(TERM2GENE = bind.info[, c("GO", "Gene", "Ontology")], TERM2NAME = extra.otp)
  return(list(TERM2GENE = bind.info[, c("GO", "Gene", "Ontology")], TERM2NAME = extra.otp))
}


##' @rdname read-gaf
##' @export
parse_gff <- read.gaf

# only read the file with selected columns
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
