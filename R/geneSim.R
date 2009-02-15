`geneSim` <-
function(gene1, gene2, ont="MF", drop="IEA", measure="Resnik"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	
	go1 <- ygcGetOnt(gene1, ontology= wh_ont, dropCodes=drop)
	go2 <- ygcGetOnt(gene2, ontology= wh_ont, dropCodes=drop)
	if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
		return (list(geneSim=NA, GO1=go1, GO2=go2))
	}
	sim <- mgoSim(go1,go2, wh_ont, wh_measure)
	sim <- round(sim, digits=3)
	return (list(geneSim=sim, GO1=go1, GO2=go2))
}