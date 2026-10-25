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

test_that("qc_passed reports each call's QC state with usePassedOnly = FALSE (was: always TRUE)", {
    whe <- .make_longitudinal_whe()
    qc <- matrix(TRUE, 3, 3)
    qc[1, 2] <- FALSE
    qc[3, 3] <- NA
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <- qc
    all_calls <- trackFrequency(whe, usePassedOnly = FALSE)
    expect_equal(nrow(all_calls), 9L)
    ## rows run site by site within each sample, as the matrix does
    expect_identical(all_calls$qc_passed, as.vector(!is.na(qc) & qc))
    passed <- trackFrequency(whe)
    expect_equal(nrow(passed), 7L)
    expect_true(all(passed$qc_passed))
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

test_that("detectEmergingVariants() orders by host, variant and time, combining replicates", {
    em <- detectEmergingVariants(.make_two_host_whe())
    expect_s4_class(em, "DataFrame")
    expect_identical(colnames(em),
                     c("host_id", "variant_key", "time_from", "time_to",
                       "freq_from", "freq_to", "delta_freq", "qc_from"))
    expect_identical(em$host_id, c("B", "A", "A"))
    expect_identical(em$variant_key,
                     c("seg1:100:T", "seg1:100:T", "seg1:200:G"))
    expect_identical(em$time_from, c(1, 1, 1))
    expect_identical(em$time_to, c(2, 2, 2))
    expect_identical(em$freq_from, c(0.02, 0.01, 0.04))
    ## host A's replicates s4 and s5 at timepoint 2 are combined by their
    ## mean, 0.10 and 0.03, then 0.12 and 0.20
    expect_equal(em$freq_to, c(0.20, 0.065, 0.16))
    expect_equal(em$delta_freq, c(0.18, 0.055, 0.12))
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

## One site of host H1 sampled at the given timepoints
.one_site_whe <- function(time, freq) {
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(100, width = 1))
    S4Vectors::mcols(gr)$ref <- "A"
    S4Vectors::mcols(gr)$alt <- "T"
    WithinHostExperiment(
        assays = list(altFreq = matrix(freq, nrow = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = paste0("s", seq_along(time)),
            host_id = "H1", timepoint = time))
}

test_that("timepoints are put in time order and an NA timepoint is left out with a warning (was: text sorted as text, NA placed last)", {
    ## numbers given as text: "10" comes after "2"
    em <- detectEmergingVariants(
        .one_site_whe(c("1", "2", "10"), c(0.01, 0.30, 0.02)))
    expect_identical(c(em$time_from, em$time_to), c("1", "2"))
    expect_equal(c(em$freq_from, em$freq_to), c(0.01, 0.30))
    ## text without a time order is refused, as by flagTemporalInconsistency()
    day <- .one_site_whe(c("day1", "day2", "day10"), c(0.01, 0.30, 0.02))
    expect_error(detectEmergingVariants(day), "no time order")
    expect_error(trackFrequency(day), "no time order")
    ## a sample without a timepoint cannot be placed in time
    na_time <- .one_site_whe(c(1, 2, NA), c(0.20, 0.01, 0.40))
    expect_warning(em <- detectEmergingVariants(na_time),
                   "1 sample\\(s\\) with NA in 'timepoint' are left out")
    expect_equal(nrow(em), 0L)
    expect_warning(tr <- trackFrequency(na_time), "NA in 'timepoint'")
    expect_equal(tr$frequency, c(0.20, 0.01))
    ## a factor's levels give its order, and dates keep their class
    fac <- .one_site_whe(factor(c("late", "early"), levels = c("early", "late")),
                         c(0.30, 0.01))
    expect_identical(as.character(detectEmergingVariants(fac)$time_from),
                     "early")
    dates <- as.Date(c("2024-03-01", "2024-01-15"))
    em <- detectEmergingVariants(.one_site_whe(dates, c(0.30, 0.01)))
    expect_identical(em$time_from, dates[2])
    ## the plot's axis follows time order
    skip_if_not_installed("ggplot2")
    p <- plotFrequencyTrajectory(
        .one_site_whe(c("1", "2", "10"), c(0.05, 0.30, 0.15)), "H1")
    expect_identical(levels(p$data$timepoint), c("1", "2", "10"))
})

test_that("replicates at one timepoint are combined before intervals are formed (was: a zero-length interval, and an answer that depended on column order)", {
    ## 0 at t1, then two replicates at t2
    em <- detectEmergingVariants(.one_site_whe(c(1, 2, 2), c(0, 0.01, 0.12)))
    expect_identical(c(em$time_from, em$time_to), c(1, 2))
    expect_equal(c(em$freq_from, em$freq_to), c(0, 0.065))
    swapped <- .one_site_whe(c(1, 2, 2), c(0, 0.12, 0.01))
    expect_equal(detectEmergingVariants(swapped), em)
    ## replicates alone, at one timepoint, are not an interval
    expect_equal(nrow(detectEmergingVariants(
        .one_site_whe(c(5, 5), c(0.01, 0.12)))), 0L)
})

test_that("after flagISNV() a rise from below minFreq is reported and a failing timepoint is not bridged (was: starts below 3% dropped, t1 joined to t3)", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300, 400), width = 1))
    S4Vectors::mcols(gr)$ref <- rep("A", 4L)
    S4Vectors::mcols(gr)$alt <- c("T", "G", "C", "G")
    ## at t2 the third variant is at 0.60, above the default maxFreq
    freq <- matrix(c(0.00, 0.30, 0.30,
                     0.01, 0.30, 0.30,
                     0.04, 0.60, 0.45,
                     0.04, 0.30, 0.30), nrow = 4, byrow = TRUE)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = matrix(2000L, 4, 3),
                      altCount = matrix(as.integer(round(freq * 2000)), 4, 3)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("t1", "t2", "t3"),
                                       host_id = "H1", timepoint = 1:3))
    wq <- flagISNV(whe, ISNVFilter())
    em <- detectEmergingVariants(wq)
    expect_identical(em$variant_key,
                     c("seg1:100:T", "seg1:200:G", "seg1:400:G"))
    expect_identical(c(em$time_from, em$time_to), c(1L, 1L, 1L, 2L, 2L, 2L))
    expect_equal(em$freq_from, c(0, 0.01, 0.04))
    ## the two starts below minFreq failed QC, and say so
    expect_identical(em$qc_from, c(FALSE, FALSE, TRUE))
    ## without QC the third variant's rise to 0.60 is reported, from t1 to t2
    all_calls <- detectEmergingVariants(wq, usePassedOnly = FALSE)
    expect_identical(all_calls$variant_key[3], "seg1:300:C")
    expect_identical(c(all_calls$time_from[3], all_calls$time_to[3]), 1:2)
})
