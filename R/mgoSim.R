`mgoSim` <- 
function(GO1, GO2, ont="MF", organism="human", measure="Wang"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast") )
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))

	GO1 <- unlist(GO1)
	GO2 <- unlist(GO2)
	m <- length(GO1)
	n <- length(GO2)
	 
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- GO1
	colnames(scores) <- GO2
	for( i in 1:m) {
		for (j in 1:n) {
			scores[i,j] <- goSim(GO1[i], GO2[j], wh_ont, wh_organism, wh_measure)
		}
	}

	if (!sum(!is.na(scores))) return (NA)	
	if (n ==1 || m == 1) {
		return (max(scores))
	}
	
	sim <- (sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)})) + sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)})))/(m+n)	
			
	return (round(sim,digits=3))
}
