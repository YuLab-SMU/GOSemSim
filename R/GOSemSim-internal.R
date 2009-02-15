`ygcQueryAncestor` <-
function (goid, ont="MF") {
	wh_ont<-match.arg(ont, c("MF", "BP", "CC"))
	wh_ont <- tolower(wh_ont)	
	
	qid <- paste('\'', goid, '\'', sep='')
	qid <- paste("go_term.go_id=", qid, sep='')
	go.parents <- paste("go_", wh_ont, "_parents", sep='')
	go.parents.id <- paste(go.parents, "._id", sep='')
	go.parents.parents.id <- paste(go.parents, "._parent_id", sep='')
	
	SQL <- paste("SELECT go_term2.go_id,", go.parents,".relationship_type", " FROM go_term,", go.parents, ",go_term as go_term2 WHERE ",  qid, " AND go_term._id=", go.parents.id," AND ", go.parents.parents.id, "=go_term2._id;", sep='')
	
	ancestor <- dbGetQuery(GO_dbconn(), SQL)
	if( length(ancestor) ){
		ancestor1 <-  unlist(lapply (ancestor, ygcQueryAncestor))
		ancestor <- c(ancestor, ancestor1)
	} 
	return (ancestor)
}



`ygcSemVal` <-
function(query.go, weight.isa = 0.8, weight.partof = 0.6)
{
	sv <- 1
	sw <- 1
	n <- length(query.go)
	if (n > 3) {
		sequence <- seq(1, n-1, by=2)	
		w <- 1
		for(i in sequence) {
			old.w <- w
			for (j in query.go[[i+1]]) {
				if (j == "isa") {
					w <- old.w * weight.isa
				} else {
					w <- old.w*weight.partof
				}
				sw <- c(sw, w)
				sv <- sv + w
			}
		}
	}
	return (list(sv=sv, sw=sw))
}


`ygcStoreAncestor` <- function(GO, ont) {
	all_id <- unique(unlist(GO))
	Ancestor <- list()
	for (m in 1:length(all_id)) {
		Ancestor[[m]] <- ygcQueryAncestor(all_id[m], ont)
	}
	names(Ancestor) <- all_id
	return (Ancestor)
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

`WangMethod` <- 
function(GOID1, GOID2, ont="MF"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	weight.isa = 0.8
	weight.partof = 0.6
	if (GOID1 == GOID2)
		return (gosim=1)
	if (!exists("GOSemSimenvironment")) 
		.initial()
	if (exists("Ancestor", envir=GOSemSimenvironment)) {
		Ancestor <- get("Ancestor", envir=GOSemSimenvironment)	
		go1 <- Ancestor[[which(names(Ancestor)==GOID1)]]
		go2 <- Ancestor[[which(names(Ancestor)==GOID2)]]
	} else {
		go1 <- ygcQueryAncestor(GOID1, wh_ont)
		go2 <- ygcQueryAncestor(GOID2, wh_ont)
	}	
	
	size1 <- length(go1)
	size2 <- length(go2)
	if (size1 == 0 || size2 == 0)	return (NA)
	
	sv1 <- ygcSemVal(go1, weight.isa, weight.partof)
	sv2 <- ygcSemVal(go2, weight.isa, weight.partof)
	
	ancestor1 <- unlist(go1[seq(1, length(go1)-3, by=2)])
	ancestor2 <- unlist(go2[seq(1, length(go2)-3, by=2)])
	ancestor1 <- c(GOID1, ancestor1)
	ancestor2 <- c(GOID2, ancestor2)
	
	commonAncestor<-intersect(ancestor1, ancestor2)
	ca1.idx <- unlist(sapply(commonAncestor, function(x) which(ancestor1==x)))
	ca2.idx <- unlist(sapply(commonAncestor, function(x) which(ancestor2==x)))
	
	gosim <- sum(sv1$sw[ca1.idx], sv2$sw[ca2.idx]) / (sv1$sv + sv2$sv)
	gosim <- round(gosim, digits=3)
	
	return (gosim)
}

`InfoContentMethod` <- function(GOID1, GOID2, ont, measure) {
	cnt <- "Count"
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
   	    Jiang=2*log(pms)-log(p1)-log(p2), 
   	    Rel=(2*log(pms)/(log(p1)+log(p2)))*(1-pms) )	
	return(round(sim, digits=3))
}