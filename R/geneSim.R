geneSim <-
function(gene1, gene2, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg"){
	gene1 <- as.character(gene1)
	gene2 <- as.character(gene2)
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	if(!exists("GOSemSimEnv")) .initial()
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	wh_combine <- match.arg(combine, c("max", "average", "rcmax", "rcmax.avg"))
	
	go1 <- ygcGetOnt(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
	go2 <- ygcGetOnt(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
	if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
		return (list(geneSim=NA, GO1=go1, GO2=go2))
	}
	
	sim <- mgoSim(go1,go2, wh_ont, wh_organism, wh_measure, wh_combine)
	sim <- round(sim, digits=3)
	return (list(geneSim=sim, GO1=go1, GO2=go2))
}
