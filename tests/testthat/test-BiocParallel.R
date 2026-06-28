library(GOSemSim)

context("BiocParallel")

test_that("mgeneSim supports optional BiocParallel backend", {
    skip_if_not_installed("BiocParallel")

    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "MF", computeIC = FALSE)
    genes <- c("835", "5261", "241")

    serial <- mgeneSim(genes, semData = hsGO, measure = "Wang", verbose = FALSE)
    parallel <- mgeneSim(
        genes, semData = hsGO, measure = "Wang", verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_equal(parallel, serial)
})

test_that("mclusterSim supports optional BiocParallel backend", {
    skip_if_not_installed("BiocParallel")

    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "MF", computeIC = FALSE)
    clusters <- list(
        a = c("835", "5261", "241"),
        b = c("578", "582"),
        c = c("307", "308", "317")
    )

    serial <- mclusterSim(clusters, semData = hsGO, measure = "Wang")
    parallel <- mclusterSim(
        clusters, semData = hsGO, measure = "Wang",
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_equal(parallel, serial)
})
