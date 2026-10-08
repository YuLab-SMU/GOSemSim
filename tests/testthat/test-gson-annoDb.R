library(GOSemSim)

context("godata with a GSON annotation")

make_gson <- function(gsid, gene) {
    gson::gson(
        gsid2gene = data.frame(gsid = gsid, gene = gene),
        gsid2name = data.frame(gsid = gsid, name = paste("set", gsid)),
        species = "Homo sapiens",
        gsname = "Gene Ontology",
        version = "test",
        accessed_date = "2026-10-08",
        keytype = "ENTREZID"
    )
}

test_that("GO annotations can be passed as a GSON object", {
    skip_if_not_installed("gson")

    g <- make_gson(
        c("GO:0008150", "GO:0006954", "GO:0003824", "GO:0005634"),
        c("10", "10", "20", "30")
    )

    d <- godata(annoDb = g, ont = "BP", computeIC = FALSE)
    expect_true(is(d, "GOSemSimDATA"))
    expect_identical(d@ont, "BP")
    expect_identical(sort(unique(d@geneAnno$GO)), c("GO:0006954", "GO:0008150"))
    expect_identical(unique(d@geneAnno$ONTOLOGY), "BP")
    expect_identical(d@geneAnno$GENE, c("10", "10"))

    ## the GSON route and the data.frame route describe the same annotation
    df <- data.frame(
        GENE = c("10", "10"),
        GO = c("GO:0008150", "GO:0006954"),
        ONTOLOGY = "BP"
    )
    expect_identical(d@geneAnno, godata(annoDb = df, ont = "BP", computeIC = FALSE)@geneAnno)
})

test_that("the whole ontology is available, not only the requested one", {
    skip_if_not_installed("gson")

    g <- make_gson(
        c("GO:0008150", "GO:0003824", "GO:0005634"),
        c("10", "20", "30")
    )

    for (ont in c("BP", "MF", "CC")) {
        d <- godata(annoDb = g, ont = ont, computeIC = FALSE)
        expect_identical(unique(d@geneAnno$ONTOLOGY), ont)
        expect_identical(d@geneAnno$GENE, c(BP = "10", MF = "20", CC = "30")[[ont]])
    }

    ## an ontology with no annotation yields an empty table rather than an error
    mf_only <- make_gson("GO:0003824", "20")
    expect_equal(nrow(godata(annoDb = mf_only, ont = "CC", computeIC = FALSE)@geneAnno), 0)
})

test_that("GO IDs stored without the 'GO:' prefix are recognized", {
    skip_if_not_installed("gson")

    g <- make_gson(c("0008150", "0006954"), c("10", "10"))
    d <- godata(annoDb = g, ont = "BP", computeIC = FALSE)
    expect_identical(sort(unique(d@geneAnno$GO)), c("GO:0006954", "GO:0008150"))
})

test_that("gene sets that are not GO are rejected with a usable message", {
    skip_if_not_installed("gson")

    kegg <- make_gson(c("hsa00010", "hsa00020"), c("10", "20"))
    expect_error(godata(annoDb = kegg, ont = "BP"), "no GO ID")

    empty <- gson::gson(gsid2gene = data.frame(gsid = character(), gene = character()))
    expect_error(godata(annoDb = empty, ont = "BP"), "empty")

    ## unsupported input still reaches the data.frame branch untouched
    expect_error(godata(annoDb = list(a = 1), ont = "BP"), "GSON")
})

test_that("GO IDs missing from GO.db are dropped with a warning", {
    skip_if_not_installed("gson")

    g <- make_gson(c("GO:0008150", "GO:9999999"), c("10", "10"))
    expect_warning(
        d <- godata(annoDb = g, ont = "BP", computeIC = FALSE),
        "not found in"
    )
    expect_identical(unique(d@geneAnno$GO), "GO:0008150")
})

test_that("GSON provenance is recorded in the metadata slot", {
    skip_if_not_installed("gson")

    g <- make_gson("GO:0008150", "10")
    d <- godata(annoDb = g, ont = "BP", computeIC = FALSE)

    expect_identical(names(d@metadata), c("name", "value"))
    expect_identical(d@metadata$value[d@metadata$name == "species"], "Homo sapiens")
    expect_identical(d@metadata$value[d@metadata$name == "keytype"], "ENTREZID")
    ## slots the object does not carry are left out instead of stored as NA
    expect_false("urlpattern" %in% d@metadata$name)

    multi <- gson::gson(
        gsid2gene = data.frame(gsid = "GO:0008150", gene = "10"),
        species = c("Homo sapiens", "Mus musculus")
    )
    d2 <- godata(annoDb = multi, ont = "BP", computeIC = FALSE)
    expect_identical(
        d2@metadata$value[d2@metadata$name == "species"],
        "Homo sapiens, Mus musculus"
    )
})

test_that("IC is computed for a GSON annotation", {
    skip_if_not_installed("gson")

    g <- make_gson(c("GO:0008150", "GO:0006954"), c("10", "20"))
    d <- godata(annoDb = g, ont = "BP", computeIC = TRUE)
    expect_gt(length(d@IC), 0)
    expect_true(all(c("GO:0008150", "GO:0006954") %in% names(d@IC)))
})
