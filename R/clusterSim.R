`clusterSim` <- 
function(cluster1, cluster2, ont="MF", organism="human", measure="Wang", drop="IEA"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	
	size1 <- length(cluster1)
	size2 <- length(cluster2)
	if (size1 == 0 || size2 == 0) {
		return (NA)
	}
	
	allSim <- matrix(data=NA, nrow=size1, ncol=size2)
	for (i in 1:size1) {
		for (j in 1:size2){
			allSim[i,j] <- geneSim(cluster1[i], cluster2[j], wh_ont, wh_organism, wh_measure, drop)$geneSim
		}
	}	
	
	if (!sum(!is.na(allSim))) {
		return (NA)
	}
	
	result <- sum(allSim, na.rm=TRUE)/sum(!is.na(allSim))
	return (round(result, digits=3))
}