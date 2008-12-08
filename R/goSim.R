`goSim` <-
function(GOID1, GOID2, ont="MF", weight.isa = 0.8, weight.partof = 0.6){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))

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
	
	return (goSim=gosim)
}

