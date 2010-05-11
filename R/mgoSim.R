mgoSim <- 
function(GO1, GO2, ont="MF", organism="human", measure="Wang", combine="rcmax.avg"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	if(!exists("GOSemSimEnv")) .initial()
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	wh_combine <- match.arg(combine, c("max", "average", "rcmax", "rcmax.avg"))
	
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
		if (wh_combine == "average") {
			return (round(mean(scores),digits=3))
		}
		return (round(max(scores),digits=3))
	}
	
	sim <- ygcCombine(scores, wh_combine)
	return (round(sim,digits=3))
}
