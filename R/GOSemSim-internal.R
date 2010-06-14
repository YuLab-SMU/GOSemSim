ygcWangMethod <- function(GOID1, GOID2, ont="MF", organism="human") {
	if(!exists("GOSemSimEnv")) .initial()
	weight.isa = 0.8
	weight.partof = 0.6

	if (GOID1 == GOID2)
		return (gosim=1)		

	Parents.name <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"
	)
	if (!exists(Parents.name, envir=GOSemSimEnv)) {
		ygcGetParents(ont)
	}
	Parents <- get(Parents.name, envir=GOSemSimEnv)
	
	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- ygcSemVal(GOID1, ont, Parents, sv.a, sw, weight.isa, weight.partof)
	sv.b <- ygcSemVal(GOID2, ont, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- uniqsv(sv.a)
	sv.b <- uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	return(sim)
}



uniqsv <- function(sv) {
	sv <- unlist(sv)
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}

ygcSemVal_internal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	p <- Parents[goid]
	p <- unlist(p[[1]])
	if (length(p) == 0) {
		#warning(goid, " may not belong to Ontology ", ont)
		return(0)
	}
	relations <- names(p)
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "is_a") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- ygcSemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof)
		}
	}
	return (sv)
}

ygcSemVal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	if(!exists("GOSemSimCache")) return(ygcSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof))
	goid.ont <- paste(goid, ont, sep=".")
	if (!exists(goid.ont, envir=GOSemSimCache)) {
	  	value <- ygcSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof)
	  	assign(goid.ont, value, envir=GOSemSimCache)
		#cat("recompute ", goid, value, "\n")
	}
	else{
		#cat("cache ", goid, get(goid, envir=GOSemSimCache), "\n")
	}
	return(get(goid.ont, envir=GOSemSimCache))
}

ygcInfoContentMethod <- function(GOID1, GOID2, ont, measure, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	fname <- paste("Info_Contents", ont, organism, sep="_")
	tryCatch(utils::data(list=fname, package="GOSemSim", envir=GOSemSimEnv))
	IC <- get("IC", envir=GOSemSimEnv)

	p1 <- IC[GOID1]
	p2 <- IC[GOID2]
	
	if (p1 == 0 || p2 == 0) return (NA)
	Ancestor.name <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)	
	if (!exists(Ancestor.name, envir=GOSemSimEnv)) {
		ygcGetAncestors(ont)
	}
	
	Ancestor <- get(Ancestor.name, envir=GOSemSimEnv)					
	ancestor1 <- unlist(Ancestor[GOID1])
	ancestor2 <- unlist(Ancestor[GOID2])
	if (GOID1 == GOID2) {
		commonAncestor <- GOID1
	} else if (GOID1 %in% ancestor2) {
		commonAncestor <- GOID1
	} else if (GOID2 %in% ancestor1) {
		commonAncestor <- GOID2
	} else { 
		commonAncestor <- intersect(ancestor1, ancestor2)
	}
	if (length(commonAncestor) == 0) return (NA)
	
	#Information Content of the most informative common ancestor (MICA)
	pms <- max(IC[commonAncestor])  
	
	## IC is biased
	## because the IC of a term is dependent of its children but not on its parents.
	sim<-switch(measure,
   	    Resnik = pms, ## Resnik does not consider how distant the terms are from their common ancestor.
   	    ## Lin and Jiang take that distance into account.
   	    Lin = 2*pms/(p1+p2),
   	    Jiang = 1 - min(1, -2*pms + p1 + p2), 
   	    Rel = 2*pms/(p1+p2)*(1-exp(-pms))
	)   	
	return (sim)
}

ygcCombine <- function(SimMatrix, combine="average") {
	wh_combine <- match.arg(combine, c("max", "average", "rcmax", "rcmax.avg"))
        if(!sum(!is.na(SimMatrix))) return (NA)
        if (is.vector(SimMatrix)) {
                if (wh_combine == "average") {
                  return(round(mean(SimMatrix), digits=3))
                } else {
                  return (round(max(SimMatrix), digits=3)) 
                }
        }
        m <- nrow(SimMatrix)
        n <- ncol(SimMatrix)
        if (n==1 || m==1) {
                if (wh_combine == "average") {
                  return(round(mean(SimMatrix), digits=3))
                } else {
                  return (round(max(SimMatrix), digits=3)) 
                }
        }
        if (wh_combine == "average") {
		result <- mean(SimMatrix, na.rm=TRUE)
	} else if (wh_combine == "max") {
		result <- max(SimMatrix, na.rm=TRUE)
	} else if (wh_combine == "rcmax") {
		rowScore <- mean(apply(SimMatrix, 1, max, na.rm=TRUE))
		colScore <- mean(apply(SimMatrix, 2, max, na.rm=TRUE))
		result <- max(rowScore, colScore)
	} else if (wh_combine == "rcmax.avg") {
		result <- sum( apply(SimMatrix, 1, max, na.rm=TRUE), apply(SimMatrix, 2, max, na.rm=TRUE) ) / sum(dim(SimMatrix))
	}
	
	return (round(result, digits=3))
}

