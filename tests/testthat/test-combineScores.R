library(GOSemSim)

context("combineScores")

## The reference is written from the definitions in
## https://doi.org/10.1093/bioinformatics/btq406 and deliberately avoids the
## early returns of combineScores(), so that a wrong shortcut there cannot be
## mirrored here. A bare vector is treated as a single row, which is what the
## combining rules imply for one set holding one term.
combineScores_reference <- function(SimScores, combine) {
    if (length(combine) == 0) {
        return(round(SimScores, digits = 3))
    }

    if (!sum(!is.na(SimScores))) return(NA)

    if (is.null(dim(SimScores))) SimScores <- matrix(SimScores, nrow = 1)

    row.na.idx <- apply(SimScores, 1, function(i) all(is.na(i)))
    if (any(row.na.idx)) {
        SimScores <- SimScores[-which(row.na.idx), , drop = FALSE]
    }

    col.na.idx <- apply(SimScores, 2, function(i) all(is.na(i)))
    if (any(col.na.idx)) {
        SimScores <- SimScores[, -which(col.na.idx), drop = FALSE]
    }

    if (combine == "avg") {
        result <- mean(SimScores, na.rm = TRUE)
    } else if (combine == "max") {
        result <- max(SimScores, na.rm = TRUE)
    } else if (combine == "rcmax") {
        rowScore <- mean(apply(SimScores, 1, max, na.rm = TRUE))
        colScore <- mean(apply(SimScores, 2, max, na.rm = TRUE))
        result <- max(rowScore, colScore)
    } else if (combine == "rcmax.avg" || combine == "BMA") {
        result <- sum(
            apply(SimScores, 1, max, na.rm = TRUE),
            apply(SimScores, 2, max, na.rm = TRUE)
        ) / sum(dim(SimScores))
    }

    round(result, digits = 3)
}

test_that("combineScores optimized path matches previous behavior", {
    sim <- matrix(
        c(
            0.3, NA, 0.8, NA,
            0.5, 0.7, NA, NA,
            NA, 0.4, 0.9, NA,
            NA, NA, NA, NA
        ),
        nrow = 4,
        byrow = TRUE
    )

    for (combine in c("avg", "max", "rcmax", "rcmax.avg", "BMA")) {
        expect_equal(
            combineScores(sim, combine),
            combineScores_reference(sim, combine)
        )
    }
})

test_that("combineScores agrees with the reference on degenerate shapes", {
    set.seed(1)
    shapes <- list(
        matrix(c(1, 0.4, 0.2), nrow = 1),
        matrix(c(1, 0.4, 0.2), ncol = 1),
        matrix(c(1), 1, 1),
        matrix(round(runif(6), 3), nrow = 1),
        matrix(round(runif(6), 3), ncol = 1),
        c(1, 0.4, 0.2),
        c(0.9, 0.1, 0.5, 0.3)
    )
    for (sim in shapes) {
        for (combine in c("avg", "max", "rcmax", "rcmax.avg", "BMA")) {
            expect_equal(
                combineScores(sim, combine),
                combineScores_reference(sim, combine)
            )
        }
    }
})

## hand-computed values, so that a regression is caught even if the reference
## above is ever changed to match a wrong implementation again
test_that("BMA on a single row or column is (max + sum) / (n + 1), not the max", {
    S <- matrix(c(1, 0.4, 0.2), nrow = 1)

    expect_equal(combineScores(S, "BMA"), 0.65)        # (1 + 1.6) / 4
    expect_equal(combineScores(S, "rcmax.avg"), 0.65)
    expect_equal(combineScores(t(S), "BMA"), 0.65)     # the same as one column
    expect_equal(combineScores(c(1, 0.4, 0.2), "BMA"), 0.65)  # and as a vector
    expect_equal(combineScores(matrix(0.7, 1, 1), "BMA"), 0.7) # (0.7 + 0.7) / 2

    ## 'max' and 'rcmax' do equal the max in these shapes, so they stay put
    expect_equal(combineScores(S, "max"), 1)
    expect_equal(combineScores(S, "rcmax"), 1)
    expect_equal(combineScores(S, "avg"), round(mean(c(1, 0.4, 0.2)), 3))

    ## a single row must not score 1 against a set it merely overlaps
    expect_lt(combineScores(matrix(c(1, 0.1411), nrow = 1), "BMA"), 1)
})

test_that("2-D BMA is unchanged by the degenerate fix", {
    m <- matrix(c(0.9, 0.2, 0.4, 0.1, 0.8, 0.3, 0.5, 0.6, 0.7), nrow = 3, byrow = TRUE)
    expect_equal(combineScores(m, "BMA"), 0.8)   # (2.4 + 2.4) / 6
    expect_equal(combineScores(m, "rcmax.avg"), 0.8)
})

test_that("NA cells are dropped before combining, as in the 2-D path", {
    S <- matrix(c(0.3, NA, 0.8), nrow = 1)
    expect_equal(combineScores(S, "BMA"), round((0.8 + 1.1) / 3, 3))
})

test_that("subsetCombine keeps a single-term set as a one-row matrix", {
    go_matrix <- matrix(c(1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2,
                        dimnames = list(c("a", "b"), c("a", "b")))

    ## row 'a' against both terms is 1 x 2, and must get the BMA
    expect_equal(GOSemSim:::subsetCombine(go_matrix, "a", c("a", "b"), "BMA"),
                 round((1 + 1.3) / 3, 3))
    expect_equal(GOSemSim:::subsetCombine(go_matrix, "a", "a", "BMA"), 1)
    expect_equal(GOSemSim:::subsetCombine(go_matrix, "a", c("a", "b"), "max"), 1)
})
