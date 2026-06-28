library(GOSemSim)

context("Wang")

test_that("Wang's method", {
    hsGO <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
    x <- goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
    expect_true(x >= 0 && x<=1)
})

test_that("Wang semantic value cache is ontology-specific", {
    yulab.utils::initial_cache()
    rel_df <- data.frame(
        go_id = c("GO:test", "GO:test"),
        relationship = c("is_a", "is_a"),
        parent = c("GO:bp_parent", "GO:mf_parent"),
        Ontology = c("BP", "MF"),
        stringsAsFactors = FALSE
    )

    bp_sv <- GOSemSim:::getSV("GO:test", "BP", rel_df)
    mf_sv <- GOSemSim:::getSV("GO:test", "MF", rel_df)

    expect_equal(unname(bp_sv["GO:bp_parent"]), 0.8)
    expect_true(is.na(bp_sv["GO:mf_parent"]))
    expect_equal(unname(mf_sv["GO:mf_parent"]), 0.8)
    expect_true(is.na(mf_sv["GO:bp_parent"]))
})
