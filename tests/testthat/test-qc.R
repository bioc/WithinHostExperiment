# tests/testthat/test-qc.R

library(S4Vectors)
library(SummarizedExperiment)

# ---- Helper: build test WHE from VCF ----

.load_donor_whe <- function() {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
}

# ============================================================
# flagISNV tests
# ============================================================

test_that("flagISNV returns WithinHostExperiment", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter())
    expect_s4_class(result, "WithinHostExperiment")
})

test_that("flagISNV records steps in qcLog", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter())
    log <- qcLog(result)
    expect_true(nrow(log) >= 4L)
    expect_true("depth_filter" %in% log$step)
    expect_true("min_freq_filter" %in% log$step)
})

test_that("flagISNV stores filter in qcFilters", {
    whe <- .load_donor_whe()
    filt <- ISNVFilter(minDepth = 500L)
    result <- flagISNV(whe, filt)
    expect_equal(length(qcFilters(result)), 1L)
    expect_s4_class(qcFilters(result)[[1]], "ISNVFilter")
})

test_that("flagISNV depth filter flags low-depth sites", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter(minDepth = 200L))
    qc <- assay(result, "qcPass")[, 1]
    depth <- assay(result, "totalDepth")[, 1]
    low_depth <- which(!is.na(depth) & depth < 200L)
    expect_true(length(low_depth) >= 1L)
    expect_true(all(!qc[low_depth]))
})

test_that("flagISNV minFreq filter flags low-frequency sites", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter(minFreq = 0.10))
    qc <- assay(result, "qcPass")[, 1]
    freq <- assay(result, "altFreq")[, 1]
    low_freq <- which(!is.na(freq) & freq < 0.10)
    expect_true(length(low_freq) >= 1L)
    expect_true(all(!qc[low_freq]))
})

test_that("flagISNV maxFreq filter flags high-frequency sites", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter(maxFreq = 0.40))
    qc <- assay(result, "qcPass")[, 1]
    freq <- assay(result, "altFreq")[, 1]
    high_freq <- which(!is.na(freq) & freq > 0.40)
    expect_true(length(high_freq) >= 1L)
    expect_true(all(!qc[high_freq]))
})

test_that("flagISNV does NOT remove rows", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter(minDepth = 5000L))
    expect_equal(nrow(result), nrow(whe))
    expect_true(any(!assay(result, "qcPass")))
})

test_that("passedISNV returns fewer rows after flagISNV", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter(minDepth = 200L))
    passed <- passedISNV(result)
    expect_true(nrow(passed) < nrow(whe))
    expect_true(nrow(passed) > 0L)
})

test_that("flagISNV with default filter keeps most sites", {
    whe <- .load_donor_whe()
    result <- flagISNV(whe, ISNVFilter())
    qc <- assay(result, "qcPass")[, 1]
    expect_true(sum(qc, na.rm = TRUE) >= 10L)
})

test_that("sequential flagISNV calls accumulate qcLog", {
    whe <- .load_donor_whe()
    whe <- flagISNV(whe, ISNVFilter(minDepth = 100L))
    n1  <- nrow(qcLog(whe))
    whe <- flagISNV(whe, ISNVFilter(minFreq = 0.10))
    n2  <- nrow(qcLog(whe))
    expect_true(n2 > n1)
})

# ============================================================
# qcSummary tests
# ============================================================

test_that("qcSummary returns correct columns", {
    whe <- .load_donor_whe()
    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    smry <- qcSummary(whe)
    expect_s4_class(smry, "DataFrame")
    expect_true(all(c("sample_id", "total_sites", "passed_sites",
                      "flagged_sites", "pass_rate") %in% colnames(smry)))
})

test_that("qcSummary counts match expectations", {
    whe <- .load_donor_whe()
    whe <- flagISNV(whe, ISNVFilter(minDepth = 200L))
    s   <- qcSummary(whe)
    expect_equal(s$total_sites, 15L)
    expect_equal(s$passed_sites + s$flagged_sites, s$total_sites)
    expect_true(s$pass_rate > 0 & s$pass_rate <= 1)
})

# ============================================================
# flagReplicateDiscordance tests
# ============================================================

test_that("flagReplicateDiscordance works with replicate VCFs", {
    rep1 <- system.file("extdata", "test_donor_rep1.vcf",
                        package = "WithinHostExperiment")
    rep2 <- system.file("extdata", "test_donor_rep2.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(rep1, rep2),
        colData = DataFrame(
            sample_id       = c("rep1", "rep2"),
            replicate_group = c("donor_reps", "donor_reps"),
            replicate_id    = c("1", "2")),
        caller = "ivar")
    result <- flagReplicateDiscordance(whe, freqTolerance = 0.05)
    expect_s4_class(result, "WithinHostExperiment")
    log <- qcLog(result)
    expect_true("replicate_discordance" %in% log$step)
})

test_that("flagReplicateDiscordance flags known discordant sites", {
    rep1 <- system.file("extdata", "test_donor_rep1.vcf",
                        package = "WithinHostExperiment")
    rep2 <- system.file("extdata", "test_donor_rep2.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(rep1, rep2),
        colData = DataFrame(
            sample_id       = c("rep1", "rep2"),
            replicate_group = c("donor_reps", "donor_reps"),
            replicate_id    = c("1", "2")),
        caller = "ivar")
    result <- flagReplicateDiscordance(whe, freqTolerance = 0.05)
    qc <- assay(result, "qcPass")
    total_flagged <- sum(!qc, na.rm = TRUE)
    expect_true(total_flagged >= 2L)
})

test_that("flagReplicateDiscordance flags rep-only sites", {
    rep1 <- system.file("extdata", "test_donor_rep1.vcf",
                        package = "WithinHostExperiment")
    rep2 <- system.file("extdata", "test_donor_rep2.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(rep1, rep2),
        colData = DataFrame(
            sample_id       = c("rep1", "rep2"),
            replicate_group = c("donor_reps", "donor_reps"),
            replicate_id    = c("1", "2")),
        caller = "ivar")
    result <- flagReplicateDiscordance(whe, freqTolerance = 0.05,
                                      requireBothDetected = TRUE)
    qc    <- assay(result, "qcPass")
    freq1 <- assay(result, "altFreq")[, 1]
    freq2 <- assay(result, "altFreq")[, 2]
    rep_only <- which(!is.na(freq1) & is.na(freq2))
    if (length(rep_only) > 0L) {
        expect_true(any(!qc[rep_only, 1]))
    }
})

test_that("flagReplicateDiscordance skips when no replicate_group", {
    whe <- .load_donor_whe()
    expect_warning(
        result <- flagReplicateDiscordance(whe),
        "replicate_group"
    )
    expect_true(all(assay(result, "qcPass")))
})

test_that("flagReplicateDiscordance with high tolerance flags fewer", {
    rep1 <- system.file("extdata", "test_donor_rep1.vcf",
                        package = "WithinHostExperiment")
    rep2 <- system.file("extdata", "test_donor_rep2.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(rep1, rep2),
        colData = DataFrame(
            sample_id       = c("rep1", "rep2"),
            replicate_group = c("donor_reps", "donor_reps"),
            replicate_id    = c("1", "2")),
        caller = "ivar")
    strict  <- flagReplicateDiscordance(whe, freqTolerance = 0.03)
    lenient <- flagReplicateDiscordance(whe, freqTolerance = 0.20)
    strict_flags  <- sum(!assay(strict, "qcPass"), na.rm = TRUE)
    lenient_flags <- sum(!assay(lenient, "qcPass"), na.rm = TRUE)
    expect_true(strict_flags >= lenient_flags)
})

# ============================================================
# Combined workflow test
# ============================================================

test_that("flagISNV + flagReplicateDiscordance pipeline works", {
    rep1 <- system.file("extdata", "test_donor_rep1.vcf",
                        package = "WithinHostExperiment")
    rep2 <- system.file("extdata", "test_donor_rep2.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(rep1, rep2),
        colData = DataFrame(
            sample_id       = c("rep1", "rep2"),
            replicate_group = c("donor_reps", "donor_reps"),
            replicate_id    = c("1", "2")),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter(minDepth = 200L, minFreq = 0.03))
    whe <- flagReplicateDiscordance(whe, freqTolerance = 0.05)
    log <- qcLog(whe)
    expect_true("depth_filter" %in% log$step)
    expect_true("replicate_discordance" %in% log$step)
    passed <- passedISNV(whe)
    expect_true(nrow(passed) > 0L)
    expect_true(nrow(passed) < nrow(whe))
})

test_that("a qcPass that is not a base logical matrix is refused (was: silent on backends)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3, width = 1))
    build <- function(qc) WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1),
                      qcPass = qc),
        rowRanges = gr,
        colData = DataFrame(sample_id = "S1"))
    ## an ordinary matrix of the wrong type
    expect_error(build(matrix(c(1L, 1L, 0L), ncol = 1)),
                 "'qcPass' assay must be a logical matrix")
    ## an array-like object that is not a base matrix, as a DelayedArray is
    whe <- build(matrix(TRUE, 3, 1))
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <-
        array(c(TRUE, TRUE, FALSE), c(3, 1, 1))
    expect_error(validObject(whe),
                 "must be an ordinary R logical matrix; DelayedArray")
})

test_that("flagReplicateDiscordance() checks requireBothDetected (was: any value taken as TRUE)", {
    ## the same checks as ConcordanceRule(), so the two routes agree
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.4, 0.1, 0.1), ncol = 2,
                                       dimnames = list(NULL, c("a", "b")))),
        rowRanges = GenomicRanges::GRanges("seg",
            IRanges::IRanges(c(10, 20), width = 1)),
        colData = DataFrame(sample_id = c("a", "b"),
                            replicate_group = c("g", "g"),
                            row.names = c("a", "b")))
    expect_error(flagReplicateDiscordance(whe, requireBothDetected = 0.5),
                 "'requireBothDetected' must be TRUE or FALSE")
    expect_error(flagReplicateDiscordance(whe, requireBothDetected = NA),
                 "'requireBothDetected' must be TRUE or FALSE")
    expect_error(flagReplicateDiscordance(whe, freqTolerance = NA_real_),
                 "'freqTolerance' must be a single positive number")
    expect_error(flagReplicateDiscordance(whe, freqTolerance = 1,
                                          requireBothDetected = FALSE),
                 "flags nothing")
    ## replicate intersection is still allowed
    expect_s4_class(flagReplicateDiscordance(whe, freqTolerance = Inf),
                    "WithinHostExperiment")
})
