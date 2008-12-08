`mgoSim` <- 
function(GO1, GO2, ont="MF", weight.isa = 0.8, weight.partof = 0.6){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	
	GO1 <- unlist(GO1)
	GO2 <- unlist(GO2)
	m <- length(GO1)
	n <- length(GO2)
	allGO <- unlist(c(GO1,GO2))
	flag <- 0
	if (!exists("GOSemSimenvironment")) 
		.initial()
	if (!exists("Ancestor", envir=GOSemSimenvironment)) {
		Ancestor <- assign("Ancestor", GOSemSimenvironment)
		Ancestor <- ygcStoreAncestor(allGO, wh_ont)
		flag <- 1
	} 
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- GO1
	colnames(scores) <- GO2
	for( i in 1:m) {
		for (j in 1:n) {
			if(is.na(scores[i,j])) {
				scores[i,j] <- goSim(GO1[i], GO2[j])
			} 
		}
	}
	if (n ==1 || m == 1) {
		return (max(scores))
	}
	
	mgosim <- (sum(sapply(1:m,function(x) {max(scores[x,])})) + sum(sapply(1:n, function(x) {max(scores[,x])})))/(m+n) 
	if (flag==1 && exists("Ancestor", envir=GOSemSimenvironment)) {
		rm(Ancestor)
	}
	return (mgosim)
}