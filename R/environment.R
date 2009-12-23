.initial <- function() {
	assign("GOSemSimEnv", new.env(),.GlobalEnv)
	assign("GOSemSimCache", new.env(), .GlobalEnv)
}
