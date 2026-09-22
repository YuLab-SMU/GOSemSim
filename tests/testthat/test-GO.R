library(GOSemSim)
library(GO.db)

context("GO")

test_that("parent node", {
    goid <- 'GO:0004022'
    x <- as.list(GOSemSim:::getParents('MF')[goid])
    expect_equal(names(x),goid)

    pid <- x[[1]]
    expect_true( goid %in% as.list(GOMFCHILDREN[pid])[[1]] )
})

test_that("ontology mappings are cached as a whole, and match GO.db", {
    ## Each getter used to rebuild its entire mapping on every call, which cost
    ## about half a second, while callers only ever ask for one or two terms.
    ## The cache key is therefore the ontology, not the individual term.
    yulab.utils::initial_cache()

    expect_null(yulab.utils::get_cache_element("GOSemSim_ontoRelation", "ancestor|BP"))

    expect_identical(GOSemSim:::getAncestors("BP"),
                     AnnotationDbi::as.list(GOBPANCESTOR))
    expect_identical(GOSemSim:::getParents("BP"),
                     AnnotationDbi::as.list(GOBPPARENTS))
    expect_identical(GOSemSim:::getOffsprings("BP"),
                     AnnotationDbi::as.list(GOBPOFFSPRING))

    ## the mapping is now held under the ontology, so asking again -- for the
    ## same or for a different term -- returns it without rebuilding
    cached <- yulab.utils::get_cache_element("GOSemSim_ontoRelation", "ancestor|BP")
    expect_false(is.null(cached))
    expect_identical(cached, GOSemSim:::getAncestors("BP"))

    ## and the three ontologies are cached independently
    expect_false(identical(GOSemSim:::getAncestors("BP"),
                           GOSemSim:::getAncestors("MF")))
})

