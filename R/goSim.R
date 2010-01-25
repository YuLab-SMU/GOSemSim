`goSim` <-
function(GOID1, GOID2, ont="MF", organism="human", measure="Wang"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	if (wh_measure == "Wang") {
		sim <- ygcWangMethod(GOID1, GOID2, ont=wh_ont, wh_organism)
	} else {
		sim <- ygcInfoContentMethod(GOID1, GOID2, ont=wh_ont, measure=wh_measure, wh_organism)
	}
	sim <- unname(sim, force=TRUE)
	return(round(sim, digits=3))
	
}

