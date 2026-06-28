library(GOSemSim)

context("combineScores")

combineScores_reference <- function(SimScores, combine) {
    if (length(combine) == 0) {
        return(round(SimScores, digits = 3))
    }

    if (!sum(!is.na(SimScores))) return(NA)

    if (is.vector(SimScores) || nrow(SimScores) == 1 || ncol(SimScores) == 1) {
        if (combine == "avg") {
            return(round(mean(SimScores, na.rm = TRUE), digits = 3))
        } else {
            return(round(max(SimScores, na.rm = TRUE), digits = 3))
        }
    }

    row.na.idx <- apply(SimScores, 1, function(i) all(is.na(i)))
    if (any(row.na.idx)) {
        SimScores <- SimScores[-which(row.na.idx), ]
    }

    if (!is.null(dim(SimScores))) {
        col.na.idx <- apply(SimScores, 2, function(i) all(is.na(i)))
        if (any(col.na.idx)) {
            SimScores <- SimScores[, -which(col.na.idx)]
        }
    }
    if (is.vector(SimScores) || nrow(SimScores) == 1 || ncol(SimScores) == 1) {
        if (combine == "avg") {
            return(round(mean(SimScores, na.rm = TRUE), digits = 3))
        } else {
            return(round(max(SimScores, na.rm = TRUE), digits = 3))
        }
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
