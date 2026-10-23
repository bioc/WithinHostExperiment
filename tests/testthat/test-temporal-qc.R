# tests/testthat/test-temporal-qc.R

library(GenomicRanges)
library(S4Vectors)

test_that("flagTemporalInconsistency flags transient iSNVs", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200, 300), width = 1))
    mcols(gr)$ref <- c("A", "C", "G")
    mcols(gr)$alt <- c("T", "G", "A")
    freq <- matrix(c(0.1, NA,  0.2,
                     0.05, NA,  NA,
                     0.15, 0.2, 0.1),
                   nrow = 3, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))

    whe2 <- flagTemporalInconsistency(whe, minTimepoints = 2L)
    tc <- mcols(SummarizedExperiment::rowRanges(whe2))$temporal_class
    expect_equal(tc[1], "persistent")   # 2 timepoints
    expect_equal(tc[2], "transient")    # 1 timepoint
    expect_equal(tc[3], "persistent")   # 3 timepoints
    # transient should be flagged
    qc <- SummarizedExperiment::assay(whe2, "qcPass")
    expect_true(all(!qc[2, ]))  # row 2 all FALSE
})

test_that("flagTemporalInconsistency with minConsecutive", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200), width = 1))
    mcols(gr)$ref <- c("A", "C")
    mcols(gr)$alt <- c("T", "G")
    # var1: days 1,3 (not consecutive), var2: days 1,2,3
    freq <- matrix(c(0.1,  0.15,
                     NA,   0.2,
                     0.12, 0.1),
                   nrow = 2, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))

    whe2 <- flagTemporalInconsistency(whe, minTimepoints = 2L,
                                       minConsecutive = 2L)
    tc <- mcols(SummarizedExperiment::rowRanges(whe2))$temporal_class
    expect_equal(tc[1], "transient")    # 2 detections but not consecutive
    expect_equal(tc[2], "persistent")   # 3 consecutive
})

.make_rep_whe <- function(freq, qc = NULL,
                          times = c(1, 1, 2),
                          hosts = rep("h1", length(times))) {
    n <- nrow(freq)
    gr <- GRanges("seg", IRanges::IRanges(seq_len(n) * 100, width = 1))
    mcols(gr)$ref <- rep("A", n)
    mcols(gr)$alt <- rep("T", n)
    assays <- list(altFreq = freq)
    if (!is.null(qc)) assays$qcPass <- qc
    WithinHostExperiment(
        assays = assays, rowRanges = gr,
        colData = DataFrame(
            sample_id = paste0("s", seq_along(times)),
            host_id = hosts, timepoint = times))
}

test_that("flagTemporalInconsistency() counts timepoints, not replicate samples (was: persistent)", {
    ## Samples: two replicates at timepoint 1, one sample at timepoint 2
    freq <- rbind(
        c(0.10, 0.12, NA),     # timepoint 1 only, in both replicates
        c(0.10, NA,   0.20),   # timepoints 1 and 2
        c(NA,   0.08, 0.05),   # timepoints 1 and 2
        c(NA,   NA,   0.30))   # timepoint 2 only
    whe <- flagTemporalInconsistency(.make_rep_whe(freq))
    tc <- mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
    expect_identical(tc, c("transient", "persistent", "persistent",
                           "transient"))
    qc <- SummarizedExperiment::assay(whe, "qcPass")
    expect_true(all(!qc[c(1, 4), ]))
    expect_true(all(qc[c(2, 3), ]))
    log <- qcLog(whe)
    expect_identical(
        as.integer(log$n_flagged[log$step == "temporal_consistency"]), 6L)
})

test_that("flagTemporalInconsistency() leaves a host sampled at one timepoint unclassified", {
    ## Two replicates of one timepoint are not a longitudinal series
    freq <- rbind(c(0.10, 0.12), c(0.10, NA))
    whe <- flagTemporalInconsistency(.make_rep_whe(freq, times = c(1, 1)))
    tc <- mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
    expect_identical(tc, c(NA_character_, NA_character_))
    expect_true(all(SummarizedExperiment::assay(whe, "qcPass")))
})

test_that("flagTemporalInconsistency() treats NA qcPass as not passed (was: error)", {
    freq <- rbind(c(0.10, 0.20, 0.15),
                  c(0.10, 0.20, 0.15))
    qc <- rbind(c(TRUE, NA, TRUE),
                c(TRUE, NA, NA))
    whe <- flagTemporalInconsistency(
        .make_rep_whe(freq, qc = qc, times = c(1, 2, 3)))
    tc <- mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
    ## Site 1 passed at timepoints 1 and 3; site 2 at timepoint 1 only
    expect_identical(tc, c("persistent", "transient"))
    out <- SummarizedExperiment::assay(whe, "qcPass")
    expect_identical(unname(out[1, ]), c(TRUE, NA, TRUE))
    ## The passed cell is flagged; the NA cells stay NA
    expect_identical(unname(out[2, ]), c(FALSE, NA, NA))
    log <- qcLog(whe)
    expect_identical(
        as.integer(log$n_flagged[log$step == "temporal_consistency"]), 1L)
})

test_that("flagTemporalInconsistency() applies minConsecutive to distinct timepoints", {
    ## Timepoints 1, 1, 2, 3: site 1 is seen in both replicates of
    ## timepoint 1 and at timepoint 3, which are not consecutive
    freq <- rbind(c(0.10, 0.10, NA,   0.20),
                  c(0.10, NA,   0.10, NA))
    whe <- flagTemporalInconsistency(
        .make_rep_whe(freq, times = c(1, 1, 2, 3)),
        minTimepoints = 2L, minConsecutive = 2L)
    tc <- mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
    expect_identical(tc, c("transient", "persistent"))
})

test_that("flagTemporalInconsistency() warns about samples with an NA timepoint", {
    freq <- rbind(c(0.10, NA, 0.20), c(0.10, 0.20, NA))
    expect_warning(
        whe <- flagTemporalInconsistency(
            .make_rep_whe(freq, times = c(1, 2, NA))),
        "1 sample\\(s\\) with NA in 'timepoint'")
    tc <- mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
    expect_identical(tc, c("transient", "persistent"))
    ## The flag covers every sample of the host, the undated one too
    qc <- SummarizedExperiment::assay(whe, "qcPass")
    expect_identical(unname(qc[1, ]), c(FALSE, FALSE, FALSE))
})

test_that("flagTemporalInconsistency logs the step", {
    gr <- GRanges("seg", IRanges::IRanges(100, width = 1))
    freq <- matrix(c(0.1, NA, 0.1), nrow = 1, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))
    whe2 <- flagTemporalInconsistency(whe)
    log <- qcLog(whe2)
    expect_true("temporal_consistency" %in% log$step)
})
