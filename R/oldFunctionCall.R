goSim <- function(GOID1, GOID2, ont="MF", organism="human", measure="Wang"){
    params <- new("Params", ontology=ont, organism=organism, method=measure)
    gos <- new("GOSet", GOSet1=GOID1, GOSet2=GOID2)
    result <- sim(gos, params)
    result <- as.vector(result)
    return(round(result, digits=3))
}

                                        #goSim("GO:0043121", "GO:0019838", ont="MF", organism="human", measure="Wang")

mgoSim <- function(GO1, GO2, ont="MF", organism="human", measure="Wang", combine="rcmax.avg"){
    params <- new("Params", ontology=ont, organism=organism, method=measure, combine=combine)
    gos <- new("GOSet", GOSet1=GO1, GOSet2=GO2)
    result <- sim(gos, params)
    return(round(result, digits=3))
}

                                        #go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
                                        #go2 <- c("GO:0009055", "GO:0020037")
                                        #mgoSim("GO:0003824", go2, measure="Wang")

geneSim <- function(gene1, gene2, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg"){
    params <- new("Params", ontology=ont, organism=organism, method=measure, combine=combine, dropCodes=drop)
    gs <- new("GeneSet", GeneSet1 = gene1, GeneSet2=gene2)
    result <- sim(gs, params)
    result <- round(result, digits=3)
    go1 <- gene2GO(gene1, params)
    go2 <- gene2GO(gene2, params)
    return (list(geneSim=result, GO1=go1, GO2=go2))
}

                                        #geneSim("241", "251", ont="MF", organism="human", measure="Wang")

mgeneSim <- function (genes, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg") {
    params <- new("Params", ontology=ont, organism=organism, method=measure, combine=combine, dropCodes=drop)
    gs <- new("GeneSet", GeneSet1 = genes, GeneSet2=genes)
    result <- sim(gs, params)
    return(round(result, digits=3))
}

                                        #mgeneSim(genes=c("835", "5261","241", "994"), ont="MF", organism="human", measure="Wang")

clusterSim <- function(cluster1, cluster2, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg"){
    params <- new("Params", ontology=ont, organism=organism, method=measure, combine=combine, dropCodes=drop)
    geneClusters <- new("GeneClusterSet", GeneClusters=list(cluster1=cluster1, cluster2=cluster2))
    result <- sim(geneClusters, params)
    result <- result[1,2]
    return(round(result, digits=3))
}

mclusterSim <- function(clusters, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg") {
    params <- new("Params", ontology=ont, organism=organism, method=measure, combine=combine, dropCodes=drop)
    geneClusters <- new("GeneClusterSet", GeneClusters=clusters)
    result <- sim(geneClusters, params)
    return(round(result, digits=3))
}
