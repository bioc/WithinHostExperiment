# tests/testthat/test-strand-bias.R

library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

# ---- Helper: build WHE with strand assays ----

.make_strand_whe <- function(fwd_alt, rev_alt, fwd_total, rev_total) {
    n <- length(fwd_alt)
    gr <- GRanges("seg", IRanges::IRanges(start = seq_len(n), width = 1))
    assays <- list(
        altFreq       = matrix(0.1, nrow = n, ncol = 1),
        totalDepth    = matrix(1000L, nrow = n, ncol = 1),
        qcPass        = matrix(TRUE, nrow = n, ncol = 1),
        fwdAltCount   = matrix(fwd_alt, ncol = 1),
        revAltCount   = matrix(rev_alt, ncol = 1),
        fwdTotalCount = matrix(fwd_total, ncol = 1),
        revTotalCount = matrix(rev_total, ncol = 1)
    )
    WithinHostExperiment(
        rowRanges = gr, assays = assays,
        colData = DataFrame(sample_id = "s1"))
}

# ============================================================
# sciStrandBias tests
# ============================================================

test_that("sciStrandBias returns zero for balanced strands", {
    ## When per-site ratio matches the global ratio, SCI ~ 0
    sci <- sciStrandBias(
        fwd_alt   = c(50, 50, 50),
        rev_alt   = c(50, 50, 50),
        fwd_total = c(500, 500, 500),
        rev_total = c(500, 500, 500)
    )
    expect_length(sci, 3L)
    expect_true(all(abs(sci) < 0.1))
})

test_that("sciStrandBias detects strong forward bias", {
    ## Site with nearly all reads on fwd strand
    sci <- sciStrandBias(
        fwd_alt   = c(100, 50),
        rev_alt   = c(1,   50),
        fwd_total = c(500, 500),
        rev_total = c(500, 500)
    )
    ## First site should have large positive SCI
    expect_true(sci[1] > 1.5)
    ## Second site should be near zero
    expect_true(abs(sci[2]) < 0.5)
})

test_that("sciStrandBias() refuses negative counts (was: -Inf)", {
    expect_error(sciStrandBias(-1, 2, 3, 4),
                 "'fwd_alt' must hold non-negative, finite counts")
    expect_error(sciStrandBias(1, 2, 3, c(-4)),
                 "'rev_total' must hold non-negative, finite counts")
})

test_that("sciStrandBias() refuses non-finite and non-numeric counts", {
    expect_error(sciStrandBias(Inf, 2, 3, 4), "'fwd_alt'.*finite")
    expect_error(sciStrandBias(1, NaN, 3, 4), "'rev_alt'.*finite")
    expect_error(sciStrandBias(1, 2, "3", 4),
                 "'fwd_total' must be numeric")
    expect_error(sciStrandBias(TRUE, 2, 3, 4),
                 "'fwd_alt' must be numeric")
})

# ============================================================
# flagStrandBiasSCI tests
# ============================================================

test_that("flagStrandBiasSCI flags biased sites and updates qcLog", {
    whe <- .make_strand_whe(
        fwd_alt   = c(100L, 50L, 2L),
        rev_alt   = c(1L,   50L, 40L),
        fwd_total = c(500L, 500L, 500L),
        rev_total = c(500L, 500L, 500L)
    )
    result <- flagStrandBiasSCI(whe, maxSCI = 1.5)

    expect_s4_class(result, "WithinHostExperiment")
    qc <- assay(result, "qcPass")[, 1]
    ## Site 1 (strong fwd bias) and site 3 (strong rev bias) flagged
    expect_false(qc[1])
    expect_true(qc[2])
    expect_false(qc[3])
    ## qcLog records the step
    log <- qcLog(result)
    expect_true("strand_bias_sci_filter" %in% log$step)
})

test_that("sciStrandBias() takes the global ratio over sites with all four counts (was: all NA)", {
    fa <- c(100, 50, NA, 30)
    ra <- c(1, 50, 20, 28)
    ft <- c(500, 500, NA, 900)
    rt <- c(500, 500, 400, 100)
    sci <- sciStrandBias(fa, ra, ft, rt)
    ok <- c(1L, 2L, 4L)
    expect_true(is.na(sci[3]))
    expect_equal(sci[ok], sciStrandBias(fa[ok], ra[ok], ft[ok], rt[ok]))
    expect_true(all(is.na(sciStrandBias(NA_real_, 1, 2, 3))))
})

test_that("flagStrandBiasSCI() checks samples that have absent calls (was: flagged nothing)", {
    n <- 6L
    gr <- GRanges("seg", IRanges::IRanges(start = seq_len(n), width = 1))
    fa <- matrix(c(200L, 50L, 48L, 52L, 47L, 51L,
                   200L, 50L, NA, 52L, NA, 51L), ncol = 2)
    ra <- matrix(c(2L, 50L, 50L, 49L, 50L, 50L,
                   2L, 50L, NA, 49L, NA, 50L), ncol = 2)
    ft <- matrix(500L, nrow = n, ncol = 2)
    ft[c(3, 5), 2] <- NA
    rt <- ft
    qc <- matrix(TRUE, nrow = n, ncol = 2)
    qc[c(3, 5), 2] <- NA
    whe <- WithinHostExperiment(
        rowRanges = gr,
        assays = list(
            altFreq = matrix(0.1, nrow = n, ncol = 2),
            qcPass = qc,
            fwdAltCount = fa, revAltCount = ra,
            fwdTotalCount = ft, revTotalCount = rt),
        colData = DataFrame(sample_id = c("s1", "s2")))
    res <- flagStrandBiasSCI(whe, maxSCI = 1.5)
    q <- assay(res, "qcPass")
    expect_identical(unname(q[1, ]), c(FALSE, FALSE))
    expect_true(all(q[-1, 1]))
    expect_identical(unname(q[c(2, 4, 6), 2]), c(TRUE, TRUE, TRUE))
    expect_true(all(is.na(q[c(3, 5), 2])))
    log <- qcLog(res)
    expect_identical(log$n_flagged[log$step == "strand_bias_sci_filter"], 2L)
})

test_that("flagStrandBiasSCI errors on missing assays", {
    gr <- GRanges("seg", IRanges::IRanges(start = 1:2, width = 1))
    assays <- list(
        altFreq    = matrix(0.1, nrow = 2, ncol = 1),
        totalDepth = matrix(1000L, nrow = 2, ncol = 1),
        qcPass     = matrix(TRUE, nrow = 2, ncol = 1)
    )
    whe <- WithinHostExperiment(
        rowRanges = gr, assays = assays,
        colData = DataFrame(sample_id = "s1"))
    expect_error(flagStrandBiasSCI(whe), "Missing required assays")
})
