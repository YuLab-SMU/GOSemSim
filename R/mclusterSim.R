`mclusterSim` <- 
function(clusters, ont="MF", organism="human", measure="Wang", drop="IEA") {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast"))
	
	size <- length(clusters)
	
	simmat <- matrix(NA, nrow=size, ncol=size)
	rownames(simmat) <- names(clusters)
	colnames(simmat) <- names(clusters)
	
	for (i in 1:size) {
		for (j in 1:i) {
			simmat[i,j] <- clusterSim(clusters[[i]], clusters[[j]], wh_ont, wh_organism, wh_measure, drop)
			if ( i != j) {
				simmat[j, i] <- simmat[i,j]
			} 
		}
	}	
	
	removeNA <- apply(!is.na(simmat), 1, sum) > 0
	return(simmat[removeNA, removeNA])
}