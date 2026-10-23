## Helper: build a longitudinal WHE
.make_longitudinal_whe <- function() {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
    S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
    ## 3 variants x 3 time points for host H1
    freq <- matrix(c(
        0.02, 0.10, 0.05,  # t1
        0.08, 0.15, 0.04,  # t2
        0.20, 0.20, 0.03   # t3
    ), nrow = 3, ncol = 3)
    WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = c("H1_t1", "H1_t2", "H1_t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))
}

test_that("trackFrequency returns correct structure with qc_passed", {
    whe <- .make_longitudinal_whe()
    traj <- trackFrequency(whe)
    expect_s4_class(traj, "DataFrame")
    expect_true(all(c("host_id", "variant_key", "timepoint",
                      "frequency", "qc_passed") %in% colnames(traj)))
    expect_true(nrow(traj) > 0L)
    ## No qcPass assay => all qc_passed should be TRUE
    expect_true(all(traj$qc_passed))
})

test_that("trackFrequency returns all non-NA entries", {
    whe <- .make_longitudinal_whe()
    traj <- trackFrequency(whe)
    ## 3 variants x 3 time points = 9 entries (all non-NA)
    expect_equal(nrow(traj), 9L)
})

test_that("trackFrequency errors with missing hostCol", {
    whe <- .make_longitudinal_whe()
    expect_error(trackFrequency(whe, hostCol = "patient"),
                 "patient")
})

test_that("trackFrequency errors with missing timeCol", {
    whe <- .make_longitudinal_whe()
    expect_error(trackFrequency(whe, timeCol = "day"),
                 "day")
})

test_that("detectEmergingVariants finds crossing variants", {
    whe <- .make_longitudinal_whe()
    emerging <- detectEmergingVariants(whe, freqThreshold = 0.05,
                                      minIncrease = 0.05)
    expect_s4_class(emerging, "DataFrame")
    expect_true(nrow(emerging) > 0L)
    ## Variant at pos 100: 0.02 -> 0.08 (increase 0.06, from below 0.05)
    expect_true(any(grepl("100", emerging$variant_key)))
    ## All freq_from should be below threshold
    expect_true(all(emerging$freq_from < 0.05))
    ## All increases should be >= minIncrease
    expect_true(all(emerging$delta_freq >= 0.05))
})

test_that("detectEmergingVariants returns empty for no crossings", {
    whe <- .make_longitudinal_whe()
    ## Very high threshold: nothing emerges
    emerging <- detectEmergingVariants(whe, freqThreshold = 0.01,
                                      minIncrease = 0.50)
    expect_equal(nrow(emerging), 0L)
})

.make_two_host_whe <- function(hosts = c("B", "A", "B", "A", "A")) {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("T", "G")
    ## Columns interleave the hosts; s4 and s5 are replicates of host A
    ## at timepoint 2
    freq <- matrix(c(0.20, 0.30,
                     0.01, 0.04,
                     0.02, NA,
                     0.10, 0.12,
                     0.03, 0.20), nrow = 2)
    WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = paste0("s", 1:5),
            host_id = hosts,
            timepoint = c(2, 1, 1, 2, 2)))
}

test_that("detectEmergingVariants() orders by host, variant and time, keeping tied replicates in column order", {
    em <- detectEmergingVariants(.make_two_host_whe())
    expect_s4_class(em, "DataFrame")
    expect_identical(colnames(em),
                     c("host_id", "variant_key", "time_from", "time_to",
                       "freq_from", "freq_to", "delta_freq"))
    expect_identical(em$host_id, c("B", "A", "A"))
    expect_identical(em$variant_key,
                     c("seg1:100:T", "seg1:100:T", "seg1:200:G"))
    expect_identical(em$time_from, c(1, 1, 1))
    expect_identical(em$time_to, c(2, 2, 2))
    expect_identical(em$freq_from, c(0.02, 0.01, 0.04))
    expect_identical(em$freq_to, c(0.20, 0.10, 0.12))
    expect_equal(em$delta_freq, c(0.18, 0.09, 0.08))
    expect_null(rownames(em))
})

test_that("detectEmergingVariants() leaves out samples with no host (was: error)", {
    whe <- .make_two_host_whe(hosts = c("B", "A", "B", NA, "A"))
    expect_warning(em <- detectEmergingVariants(whe),
                   "1 sample\\(s\\) with NA in 'host_id' are left out")
    expect_identical(em$host_id, c("B", "A"))
    expect_identical(em$variant_key, c("seg1:100:T", "seg1:200:G"))
    expect_identical(em$freq_from, c(0.02, 0.04))
})

test_that("detectEmergingVariants validates freqThreshold", {
    whe <- .make_longitudinal_whe()
    expect_error(detectEmergingVariants(whe, freqThreshold = -1),
                 "freqThreshold")
    expect_error(detectEmergingVariants(whe, freqThreshold = 2),
                 "freqThreshold")
})

test_that("plotFrequencyTrajectory returns ggplot", {
    skip_if_not_installed("ggplot2")
    whe <- .make_longitudinal_whe()
    p <- plotFrequencyTrajectory(whe, hostId = "H1")
    expect_s3_class(p, "gg")
})

test_that("plotFrequencyTrajectory handles missing host", {
    skip_if_not_installed("ggplot2")
    whe <- .make_longitudinal_whe()
    expect_message(
        plotFrequencyTrajectory(whe, hostId = "NONEXISTENT"),
        "No trajectory data")
})
