.initial <- function() {
	assign("GOSemSimEnv", new.env(),.GlobalEnv)
	assign("GOSemSimCache", new.env(), .GlobalEnv)
	assign("SupportedSpecies", c("anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "zebrafish"), envir=GOSemSimEnv)
}
