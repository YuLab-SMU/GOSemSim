setClass(Class="Params",
         representation(
                        ontology="character",
                        organism="character",
                        method="character",
                        combine="character",
                        dropCodes="character"),
         prototype=prototype (dropCodes="NULL")
         )

setValidity("Params",
            function(object) {
		if(!exists("GOSemSimEnv")) .initial()
		organs <- get("SupportedSpecies",envir=GOSemSimEnv)
		onts <- c("MF", "BP", "CC")
		mets <- c("Resnik", "Jiang", "Lin", "Rel", "Wang")
		combines <- c("max", "avg", "rcmax", "rcmax.avg")
		if (object@ontology %in% onts) {

		} else {
                    stop("*ontology* must be one of ", paste(onts, collapse=","))
		}
		if(length(object@organism != 0)) {
                    if (object@organism %in% organs) {
                                        #loadAnnoPkg(object) ##load annotation pkg.
                    } else {
                        stop("*organism* must be one of ", paste(organs, collapse=","))
                    }
		} else {
		}
		if (object@method %in% mets) {

		} else {
                    stop("*method* must be one of ", paste(mets, collapse=","))
		}
		if(length(object@combine) != 0) {
                    if (object@combine %in% combines) {

                    } else {
                        stop("*combine* must be one of ", paste(combines, collapse=","))
                    }
		} else {

		}
		return (TRUE)
            }
            )


setMethod(
          f= "loadGOMap",
          signature= "Params",
          definition=function(params){
              annoDb <- switch(params@organism,
                              anopheles	=	"org.Ag.eg.db",
                              arabidopsis = "org.At.tair.db",
                              bovine	= "org.Bt.eg.db",
                              canine	= "org.Cf.eg.db",
                              chicken	=	"org.Gg.eg.db",
                              chimp	=	"org.Pt.eg.db",
                              ##coelicolor = "org.Sco.eg.db", ## this package is no longer supports.
                              ecolik12 = "org.EcK12.eg.db",
                              ecsakai	=	"org.EcSakai.eg.db",
                              fly = "org.Dm.eg.db",
                              human = "org.Hs.eg.db",
                              malaria	=	"org.Pf.plasmo.db",
                              mouse = "org.Mm.eg.db",
                              pig	= "org.Ss.eg.db",
                              rat = "org.Rn.eg.db",
                              rhesus	=	"org.Mmu.eg.db",
                              worm = "org.Ce.eg.db",
                              xenopus	=	"org.Xl.eg.db",
                              yeast = "org.Sc.sgd.db",
                              zebrafish = "org.Dr.eg.db",
                              )

                                        # loading annotation pakcage
              require(annoDb, character.only = TRUE)

              gomap <- switch(params@organism,
                              anopheles	=	"org.Ag.egGO",
                              arabidopsis = "org.At.tairGO",
                              bovine	= "org.Bt.egGO",
                              canine	= "org.Cf.egGO",
                              chicken	=	"org.Gg.egGO",
                              chimp	=	"org.Pt.egGO",
                              ##coelicolor = "org.Sco.egGO", ## no longer supports.
                              ecolik12 = "org.EcK12.egGO",
                              ecsakai	=	"org.EcSakai.egGO",
                              fly = "org.Dm.egGO",
                              human = "org.Hs.egGO",
                              malaria	=	"org.Pf.plasmoGO",
                              mouse = "org.Mm.egGO",
                              pig	= "org.Ss.egGO",
                              rat = "org.Rn.egGO",
                              rhesus	=	"org.Mmu.egGO",
                              worm = "org.Ce.egGO",
                              xenopus	=	"org.Xl.egGO",
                              yeast = "org.Sc.sgdGO",
                              zebrafish = "org.Dr.egGO",
                              )
              gomap <- eval(parse(text=gomap))

              print("loading GOMap...")

              mapped_genes <- mappedkeys(gomap)
              gomap = AnnotationDbi::as.list(gomap[mapped_genes])
              if (!is.null(params@dropCodes)){
                  gomap <- sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% params@dropCodes, y$Ontology)))
                  gomap <- sapply(gomap, function(x) x[2,x[1,]=="FALSE"]) ## filt out dropCodes Evidence.
                  gomap <- gomap[sapply(gomap,length) >0]
              }else {
                  gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology))
              }
              assign("gomap", gomap, envir=GOSemSimEnv)

              print("Done...")
          }
          )

setMethod(
          f= "computeIC",
          signature= "Params",
          definition=function(params){
              gomap <- get("gomap", envir=GOSemSimEnv)

######### all GO terms appearing in an given ontology ###########
              goterms <- unlist(sapply(gomap, function(x) names(x[x == params@ontology])), use.names=FALSE)

              require(GO.db)
              if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
                  assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
              }
              goids <- get("ALLGOID", envir=GOSemSimEnv)
                                        #goids <- toTable(GOTERM)

                                        # all go terms which belong to the corresponding ontology..
              goids <- unique(goids[goids[,"Ontology"] == params@ontology, "go_id"])
              gocount <- table(goterms)
              goname <- names(gocount) #goid of specific organism and selected category.

              ## ensure goterms not appearing in the specific annotation have 0 frequency..
              go.diff <- setdiff(goids, goname)
              m <- double(length(go.diff))
              names(m) <- go.diff
              gocount <- as.vector(gocount)
              names(gocount) <- goname
              gocount <- c(gocount, m)

              Offsprings.name <- switch(params@ontology,
                                        MF = "MFOffsprings",
                                        BP = "BPOffsprings",
                                        CC = "CCOffsprings"
                                        )
              if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
                  .getOffsprings(params@ontology)
              }
              Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)
              cnt <- sapply(goids,function(x){ n=gocount[ Offsprings[[x]] ]; gocount[x]+sum(n[!is.na(n)])})
              names(cnt) <- goids
                                        # the probabilities of occurrence of GO terms in a specific corpus.
              p <- cnt/sum(gocount)
                                        # IC of GO terms was quantified as the negative log likelihood.
              IC <- -log(p)

              save(IC, file=paste(paste("Info_Contents", params@organism, params@ontology, sep="_"), ".rda", sep=""), compress="xz")
          }
          )


#########################################################
## setter functions for modifing parameters.
setReplaceMethod(
                 f="setOntology",
                 signature="Params",
                 definition=function(object, value){
                     object@ontology <- value
                     return(object)
                 }
                 )

setReplaceMethod(
                 f="setOrganism",
                 signature="Params",
                 definition=function(object, value){
                     object@organism <- value
                                        #loadAnnoPkg(object)
                     loadGOMap(object)  ## update GOMap to the corresponding species...
                     return(object)
                 }
                 )

setReplaceMethod(
                 f="setMethod",
                 signature="Params",
                 definition=function(object, value){
                     object@method <- value
                     return(object)
                 }
                 )

setReplaceMethod(
                 f="setCombineMethod",
                 signature="Params",
                 definition=function(object, value){
                     object@combine <- value
                     return(object)
                 }
                 )

############

setMethod(
          f="[",
          signature=signature(x="Params",i="character"),
          definition=function(x,i,j, ..., drop=TRUE) {
              if(i=="ontology")
                  return(x@ontology)
              if(i=="organism")
                  return(x@organism)
              if(i=="method")
                  return(x@method)
              if(i=="combine")
                  return(x@combine)
          }
          )

