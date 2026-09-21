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

test_that("TCSS handles terms that belong to more than one cluster", {
    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE, processTCSS = TRUE)
    tcssdata <- hsGO@tcssdata

    ## GO is a DAG, so a term can have several meta-term ancestors and
    ## `clusid[[term]]` is then a vector rather than a scalar. Indexing `ica`
    ## with `[[` used to raise "subscript out of bounds" / "recursive indexing
    ## failed at level 2".
    n_clusters <- vapply(tcssdata$clusid, length, integer(1))
    expect_gt(sum(n_clusters > 1), 0)

    ## GO:0000018 is a common ancestor of GO:0000019 and sits in two clusters
    skip_if(!"GO:0000018" %in% names(tcssdata$clusid),
            "GO:0000018 is not part of this ontology build")

    expect_no_error(termSim("GO:0000018", "GO:0000019", hsGO, method = "TCSS"))
    res <- termSim("GO:0000018", "GO:0000019", hsGO, method = "TCSS")
    expect_true(is.finite(unname(res[1, 1])))
})

test_that("TCSS does not collapse to 1 for every term pair", {
    hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE, processTCSS = TRUE)
    tcssdata <- hsGO@tcssdata

    ## A cluster with a single member has no internal spread, so normalising it
    ## by its own max IC makes ICA identically 1 no matter how general the term
    ## is. GO:0008150 (the BP root) is such a cluster, and because it is a
    ## common ancestor of every BP pair, `max(sim_value)` used to return 1 for
    ## essentially every pair -- the whole method degenerated to a constant.
    sz <- vapply(tcssdata$meta_graph, length, numeric(1))
    singleton <- names(sz)[sz == 1]
    expect_gt(length(singleton), 0)      # the degenerate case does occur

    for (s in singleton) {
        expect_lt(tcssdata$ica[[s]], 1)
    }

    ## the root must score like the general term it is, not like a leaf
    if ("GO:0008150" %in% singleton) {
        expect_lt(tcssdata$ica[["GO:0008150"]], 1e-3)
    }

    ## and pairwise scores must actually vary instead of being constant
    set.seed(1)
    go <- sample(names(tcssdata$clusid), 30)
    m <- termSim(go, go, hsGO, method = "TCSS")
    off_diagonal <- m[upper.tri(m)]

    expect_false(any(is.na(off_diagonal)))
    expect_gt(length(unique(off_diagonal)), 1)
    expect_lt(mean(off_diagonal == 1), 0.5)
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
