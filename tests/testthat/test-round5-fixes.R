# tests/testthat/test-round5-fixes.R

test_that("QC and diversity methods reject unused arguments", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    expect_error(flagISNV(whe, ISNVFilter(), minDepht = 100L),
                 "Unused argument")
    expect_error(calcDiversity(whe, usePassedonly = FALSE),
                 "Unused argument")
    expect_error(calcNeutralityTests(whe, genomeLength = 1000L,
                                     usePassedonly = FALSE),
                 "Unused argument")
})

test_that("flagTemporalInconsistency classifies and flags per host", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("T", "G")
    ## Site 1: two timepoints in H1, one in H2. Site 2: one timepoint each.
    freq <- matrix(c(0.1, 0.05,
                     0.2, NA,
                     0.3, NA,
                     NA,  NA), nrow = 2)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = c("h1t1", "h1t2", "h2t1", "h2t2"),
            host_id = c("H1", "H1", "H2", "H2"),
            timepoint = c(1, 2, 1, 2)))
    out <- flagTemporalInconsistency(whe, minTimepoints = 2L)
    qc <- SummarizedExperiment::assay(out, "qcPass")
    ## Site 1 is persistent in H1 and transient in H2
    expect_true(all(qc[1, 1:2]))
    expect_false(any(qc[1, 3:4]))
    classes <- S4Vectors::mcols(
        SummarizedExperiment::rowRanges(out))$temporal_class
    expect_identical(classes[1], "persistent")
    expect_identical(classes[2], "transient")
})

test_that("exactBottleneck matches sites across segments and warns at maxNb", {
    gr <- GenomicRanges::GRanges(c("segA", "segB"),
        IRanges::IRanges(c(100, 100), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("T", "G")
    ## Both segments carry a variant at position 100; only segB's is
    ## transmitted, so matching by position alone would use the wrong row
    whe <- WithinHostExperiment(
        assays = list(
            altFreq = matrix(c(0.8, 0.8, 0.0, 0.5), nrow = 2),
            altCount = matrix(c(800L, 800L, 0L, 500L), nrow = 2),
            totalDepth = matrix(c(1000L, 1000L, 1000L, 1000L), nrow = 2)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role = c("donor", "recipient"),
            pair_id = c("p1", "p1")))
    res <- exactBottleneck(whe, pairId = "p1", maxNb = 40L,
                           method = "exact_bb")
    expect_equal(res$n_variants, 2L)
    expect_true(is.finite(res$Nb))
    expect_warning(exactBottleneck(whe, pairId = "p1", maxNb = 1L,
                                   method = "presence_absence"),
                   "censored")
})

test_that("trackFrequency keeps its documented order and content", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("T", "G")
    freq <- matrix(c(0.05, NA, 0.15, 0.20), nrow = 2)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = c("t1", "t2"),
            host_id = c("H1", "H1"),
            timepoint = c(1, 2)))
    traj <- trackFrequency(whe)
    expect_equal(nrow(traj), 3L)
    expect_equal(as.numeric(traj$frequency), c(0.05, 0.15, 0.20))
    expect_equal(traj$variant_key,
                 c("seg1:100:T", "seg1:100:T", "seg1:200:G"))
    expect_equal(as.numeric(traj$timepoint), c(1, 2, 2))
})

test_that("flagISNV warns about thresholds it cannot apply", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    expect_warning(flagISNV(whe, ISNVFilter(replicateConc = TRUE)),
                   "flagReplicateDiscordance")
    expect_no_warning(flagISNV(whe, ISNVFilter()))
})

test_that("the constructor requires an altFreq assay", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200), width = 1))
    expect_error(WithinHostExperiment(
        assays = list(totalDepth = matrix(c(10L, 20L), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1")),
        "altFreq")
})
