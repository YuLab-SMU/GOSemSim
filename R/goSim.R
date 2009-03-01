`goSim` <-
function(GOID1, GOID2, ont="MF", measure="Resnik"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	
	if (wh_measure == "Wang") {
		sim <- WangMethod(GOID1, GOID2, ont=wh_ont)
	} else {
		sim <- InfoContentMethod(GOID1, GOID2, ont=wh_ont, measure=wh_measure)
	}
	return(round(sim, digits=3))
}

