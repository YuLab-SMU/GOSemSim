library(GOSemSim)

context("TCSS")

test_that("TCSS self similarity is one", {
    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE, processTCSS = TRUE)
    x <- termSim("GO:0009987", "GO:0009987", hsGO, method = "TCSS")
    expect_equal(unname(x[1, 1]), 1)
})

test_that("TCSS handles invalid terms and missing tcssdata", {
    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE, processTCSS = TRUE)
    x <- termSim("BAD", "GO:0009987", hsGO, method = "TCSS")
    expect_true(is.na(x[1, 1]))

    hsGO_no_tcss <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)
    expect_error(
        termSim("GO:0009987", "GO:0009987", hsGO_no_tcss, method = "TCSS"),
        "tcssdata not found"
    )
})

test_that("TCSS ancestor cache is ontology-specific", {
    yulab.utils::initial_cache()
    bp_anc <- GOSemSim:::ancestors_envir("GO:0009987", "BP")
    mf_anc <- GOSemSim:::ancestors_envir("GO:0009987", "MF")

    expect_false(identical(bp_anc, mf_anc))
    expect_false(is.null(bp_anc))
    expect_true(is.null(mf_anc) || all(is.na(mf_anc)))
})

test_that("TCSS requires finite IC values", {
    expect_error(
        GOSemSim:::process_tcss("BP", numeric()),
        "IC data not found"
    )
    expect_error(
        GOSemSim:::process_tcss("BP", c("GO:test" = Inf)),
        "IC data not found"
    )
})

test_that("TCSS cutoff helpers validate and return numeric predictions", {
    expect_error(
        tcss_cutoff(ont = "BP", combine_method = "bad",
                    ppidata = data.frame(a = "1", b = "2", label = TRUE)),
        "'arg' should be one of"
    )

    expect_error(
        GOSemSim:::create_filtered_ppidata(
            all_pro = c("1", "2"),
            ppidata = data.frame(a = 1, b = "2", label = TRUE)
        ),
        "ppidata must be"
    )

    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)
    genes <- unique(hsGO@geneAnno[, 1])
    ppidata <- data.frame(
        a = as.character(genes[1:2]),
        b = as.character(genes[2:3]),
        label = c(TRUE, FALSE),
        stringsAsFactors = FALSE
    )
    filtered <- GOSemSim:::create_filtered_ppidata(genes, ppidata)
    pred <- GOSemSim:::computePre(
        cutoff = 3.5,
        filtered_ppidata = filtered,
        semdata = hsGO,
        combine_method = "max"
    )

    expect_type(pred, "double")
    expect_length(pred, nrow(filtered))
})

test_that("TCSS cutoff AUC/F1 helper consumes numeric predictions", {
    testthat::skip_if_not_installed("ROCR")
    filtered <- data.frame(
        a = c("1", "2", "3", "4"),
        b = c("2", "3", "4", "5"),
        label = c(TRUE, FALSE, TRUE, FALSE),
        stringsAsFactors = FALSE
    )
    res <- GOSemSim:::calc_auc_F1_score(
        predict_result = list(c(0.9, 0.2, 0.8, NA),
                              c(0.7, 0.3, 0.6, 0.1)),
        filtered_ppidata = filtered
    )

    expect_s3_class(res, "data.frame")
    expect_named(res, c("auc", "F1_score"))
    expect_equal(nrow(res), 2)
})
