`mclusterSim` <- 
function(clusters, ont="MF", organism="human", measure="Wang", drop="IEA") {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	
	size <- length(clusters)
	cluster_gos=list()
	for(i in 1:size){
		cluster_gos[[i]]=sapply(clusters[[i]], function(x) ygcGetOnt(x, organism= wh_organism, ontology= wh_ont, dropCodes=drop))
	}
	simmat <- matrix(NA, nrow=size, ncol=size)
	rownames(simmat) <- names(clusters)
	colnames(simmat) <- names(clusters)
	assign("GOSemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)
	for (i in 1:size) {
		for (j in 1:i) {
			#simmat[i,j] <- clusterSim(clusters[[i]], clusters[[j]], wh_ont, wh_organism, wh_measure, drop)
			size1 <- length(clusters[[i]])
			size2 <- length(clusters[[j]])
			if (size1 == 0 || size2 == 0) {
				return (NA)
			}

			gos1 <- cluster_gos[[i]]
			gos2 <- cluster_gos[[j]]
			allSim <- matrix(data=NA, nrow=size1, ncol=size2)
			for (m in 1:size1) {
				for (n in 1:size2){
					if(any(!is.na(gos1[m])) &&  any(!is.na(gos2[n]))){
						sim <- mgoSim(gos1[m],gos2[n], wh_ont, wh_organism, wh_measure)
						#sim <- round(sim, digits=3)
						allSim[m,n] <- sim
					}
				}
			}
			if (!sum(!is.na(allSim))) {
				return (NA)
			}
			result <- sum(allSim, na.rm=TRUE)/sum(!is.na(allSim))
			simmat[i,j] <- round(result, digits=3)
			if ( i != j) {
				simmat[j, i] <- simmat[i,j]
			} 
		}
	}	
	remove("GOSemSimCache", envir=.GlobalEnv)
	removeNA <- apply(!is.na(simmat), 1, sum) > 0
	return(simmat[removeNA, removeNA])

}