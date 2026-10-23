# tests/testthat/test-contamination.R

library(GenomicRanges)
library(S4Vectors)

test_that("detectCrossContamination finds index-hopping signature", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200, 300), width = 1))
    mcols(gr)$ref <- c("A", "C", "G")
    mcols(gr)$alt <- c("T", "G", "A")
    freq <- matrix(c(0.80, 0.05, 0.30,   # sample A
                     0.01, 0.40, 0.00),   # sample B
                   nrow = 3, ncol = 2)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("A", "B"),
            run = c("run1", "run1")))

    res <- detectCrossContamination(whe, runCol = "run")
    expect_true(nrow(res) > 0)
    # pos 100: A=80%, B=1% -> should be flagged
    hit <- res[res$position == 100 & res$source_sample == "A", ]
    expect_equal(nrow(hit), 1L)
    expect_true(hit$freq_ratio < 0.05)
})

test_that("detectCrossContamination returns empty for clean data", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200), width = 1))
    mcols(gr)$ref <- c("A", "C")
    mcols(gr)$alt <- c("T", "G")
    freq <- matrix(c(0.80, 0.30,
                     0.00, 0.00),
                   nrow = 2, ncol = 2)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(sample_id = c("A", "B")))
    res <- detectCrossContamination(whe)
    expect_equal(nrow(res), 0L)
})

test_that("detectCrossContamination rejects invalid thresholds", {
    gr <- GRanges("seg", IRanges::IRanges(100, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(0.5, nrow = 1, ncol = 1)),
        rowRanges = gr,
        colData = DataFrame(sample_id = "S1"))
    expect_error(detectCrossContamination(whe, sinkMaxFreq = 0.5),
                 "sinkMaxFreq.*sourceMinFreq")
})

test_that("detectCrossContamination respects run grouping", {
    gr <- GRanges("seg", IRanges::IRanges(100, width = 1))
    mcols(gr)$ref <- "A"; mcols(gr)$alt <- "T"
    freq <- matrix(c(0.80, 0.01), nrow = 1, ncol = 2)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("A", "B"),
            run = c("run1", "run2")))  # different runs
    res <- detectCrossContamination(whe, runCol = "run")
    expect_equal(nrow(res), 0L)  # not flagged across runs
})

.make_run_whe <- function(runs) {
    n <- length(runs)
    gr <- GRanges("seg", IRanges::IRanges(100, width = 1))
    mcols(gr)$ref <- "A"; mcols(gr)$alt <- "T"
    WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.80, rep(0.01, n - 1L)),
                                       nrow = 1, ncol = n)),
        rowRanges = gr,
        colData = DataFrame(sample_id = LETTERS[seq_len(n)], run = runs))
}

test_that("detectCrossContamination() stops on a runCol missing from colData (was: compared all samples)", {
    whe <- .make_run_whe(c("run1", "run2"))
    expect_error(detectCrossContamination(whe, runCol = "Run"),
                 "Column 'Run' not found in colData.*Available columns: sample_id, run")
})

test_that("detectCrossContamination() refuses a runCol that is not one column name", {
    whe <- .make_run_whe(c("run1", "run2"))
    for (bad in list(1, c("run", "sample_id"), NA_character_, TRUE)) {
        expect_error(detectCrossContamination(whe, runCol = bad),
                     "'runCol' must be NULL or a single column name")
        expect_error(detectCrossContamination(whe, runCol = bad),
                     "Available columns: sample_id")
    }
})

test_that("detectCrossContamination() warns when runCol has NA values, which are excluded", {
    whe <- .make_run_whe(c("run1", "run1", NA))
    expect_warning(res <- detectCrossContamination(whe, runCol = "run"),
                   "1 sample\\(s\\) with NA in 'run'")
    expect_equal(nrow(res), 1L)
    expect_identical(res$sink_sample, "B")
    expect_identical(res$run, "run1")
})
