gene2GO <- function(gene, godata, dropCodes) {
    goAnno <- godata@geneAnno
    if (! "EVIDENCE" %in% colnames(goAnno)) {
        warning("Evidence codes not found, 'drop' parameter will be ignored...")
    } else {
        goAnno <- goAnno[!goAnno$EVIDENCE %in% dropCodes,]
    }
    go <- as.character(unique(goAnno[goAnno[,1] == gene, "GO"]))
    go[!is.na(go)]
}


