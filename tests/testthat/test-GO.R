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

