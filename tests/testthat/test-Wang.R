library(GOSemSim)

context("Wang")

test_that("Wang's method", {
    hsGO <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
    x <- goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
    expect_true(x >= 0 && x<=1)
})
