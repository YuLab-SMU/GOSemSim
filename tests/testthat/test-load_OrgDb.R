library(GOSemSim)

context("load_OrgDb")

test_that("load OrgDb with package name", {
    db <- load_OrgDb("org.Hs.eg.db")
    expect_true(is(db, "OrgDb"))
})

