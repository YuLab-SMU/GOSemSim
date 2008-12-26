`geneSim` <-
function(gene1, gene2, ont="MF", drop="IEA", weight.isa = 0.8, weight.partof = 0.6){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	go1 <- ygcGetOnt(gene1, ontology= ont, dropCodes=drop)
	go2 <- ygcGetOnt(gene2, ontology= ont, dropCodes=drop)
	if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
		return (list(geneSim=NA, GO1=go1, GO2=go2))
	}
	geneSim <- mgoSim(go1,go2, wh_ont, weight.isa, weight.partof)
	geneSim <- round(geneSim, digits=3)
	return (list(geneSim=geneSim, GO1=go1, GO2=go2))
}