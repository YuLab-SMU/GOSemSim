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
	sequence <- seq(1, length(query.go)-3, by=2)	##drop the 'all' node
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
    
    if(!is.null(dropCodes)) { 
      evidence<-sapply(allGO, function(x) x$Evidence) 
      drop<-evidence %in% dropCodes
      allGO<-allGO[!drop]
    }
    
    category<-sapply(allGO, function(x) x$Ontology)
    unlist(unique(names(allGO[category %in% ontology])))
}
