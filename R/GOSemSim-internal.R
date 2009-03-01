WangMethod <- function(GOID1, GOID2, ont="MF") {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	weight.isa = 0.8
	weight.partof = 0.6

	if (GOID1 == GOID2)
		return (gosim=1)
	
	Parents <- switch(wh_ont,
		MF = AnnotationDbi::as.list(GOMFPARENTS) ,
		BP = AnnotationDbi::as.list(GOBPPARENTS) , 
		CC = AnnotationDbi::as.list(GOCCPARENTS)	)
		
	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- ygcSemVal(GOID1, Parents, sv.a, sw, weight.isa, weight.partof)
	sv.b <- ygcSemVal(GOID2, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- uniqsv(sv.a)
	sv.b <- uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- sv.a[idx]
	inter.svb <- sv.b[idx]
	sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	return(sim)
}

uniqsv <- function(sv) {
	una <- unique(names(sv))
	sv <- sapply(una, function(x) {max(sv[names(sv)==x])})
}

ygcSemVal <- function(goid, Parents, sv, w, weight.isa=0.8, weight.partof=0.6) {	
	p <- Parents[goid]
	p <- unlist(p[[1]])
	relations <- names(p)
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "isa") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- ygcSemVal(p[i], Parents, sv, w)
		}
	}
	return (sv)
}

`ygcGetOnt` <-
  function(gene, ontology, dropCodes) {
    allGO <- org.Hs.egGO[[gene]]
    if (is.null(allGO)) {
    	return (NA)
    }
    if (sum(!is.na(allGO)) == 0) {
    	return (NA)
    }
    if(!is.null(dropCodes)) { 
      evidence<-sapply(allGO, function(x) x$Evidence) 
      drop<-evidence %in% dropCodes
      allGO<-allGO[!drop]
    }
    
    category<-sapply(allGO, function(x) x$Ontology)
    return (unlist(unique(names(allGO[category %in% ontology]))))
}


`InfoContentMethod` <- function(GOID1, GOID2, ont, measure) {
	cnt <- "Count"
	GOCount <- NULL
	rm(GOCount)
	data("GOCount", package="GOSemSim")
	rootCount <- switch(ont, 
                     MF=GOCount["GO:0003674", cnt],
                     BP=GOCount["GO:0008150", cnt],
                     CC=GOCount["GO:0005575", cnt] )
	p1 <- GOCount[GOID1, cnt]/rootCount
	p2 <- GOCount[GOID2, cnt]/rootCount        
	if (p1 == 0 || p2 == 0) return (NA)
	
	Ancestor <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
		BP = AnnotationDbi::as.list(GOBPANCESTOR) , 
		CC = AnnotationDbi::as.list(GOCCANCESTOR)	)
		
	ancestor1 <- unlist(Ancestor[GOID1])
	ancestor2 <- unlist(Ancestor[GOID2])
	commonAncestor <- intersect(ancestor1, ancestor2)
	if (length(commonAncestor) == 0) return (NA)
	pms <- min(GOCount[commonAncestor,cnt], na.rm=TRUE)/rootCount
	
	sim<-switch(measure,
   	    Resnik=-log(pms),
   	    Lin=2*log(pms)/(log(p1)+log(p2)),
   	    Jiang=1/(1-2*log(pms)+log(p1)+log(p2)), 
   	    Rel=(2*log(pms)/(log(p1)+log(p2)))*(1-pms) )	
	return (sim)
}
