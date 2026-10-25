# tests/testthat/test-visualization.R
# Comprehensive tests for all visualization functions

# -- Shared fixture builders -----------------------------------------------

.make_single_whe <- function() {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
}

.make_pair_whe <- function() {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role      = c("donor", "recipient"),
            pair_id   = c("p1", "p1")),
        caller = "ivar")
}

.make_multi_whe <- function() {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    vcf3 <- system.file("extdata", "test_independent.vcf",
                        package = "WithinHostExperiment")
    readWithinHost(c(vcf1, vcf2, vcf3),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient", "independent"),
            role      = c("donor", "recipient", "independent"),
            pair_id   = c("p1", "p1", NA)),
        caller = "ivar")
}

.make_longitudinal_whe <- function() {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
    S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
    freq <- matrix(c(0.05, 0.10, 0.30,
                     0.12, 0.15, 0.25,
                     0.20, 0.08, 0.40),
                   nrow = 3, ncol = 3)
    WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id   = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))
}

# ==========================================================================
# plotFrequencySpectrum
# ==========================================================================

test_that("plotFrequencySpectrum returns ggplot for single sample", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    p <- plotFrequencySpectrum(whe)
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum returns ggplot for multiple samples", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    p <- plotFrequencySpectrum(whe)
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum density type works", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    p <- plotFrequencySpectrum(whe, type = "density")
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum threshold line is drawn", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    p <- plotFrequencySpectrum(whe, threshold = 0.03)
    expect_s3_class(p, "ggplot")
    # Check annotation layer exists
    layers <- p$layers
    expect_true(length(layers) >= 2)
})

test_that("plotFrequencySpectrum respects freqRange", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    p <- plotFrequencySpectrum(whe, freqRange = c(0, 0.5))
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum usePassedOnly filters data", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    whe <- flagISNV(whe, ISNVFilter(minFreq = 0.50))
    ## every call is below 50%, so nothing passes and nothing is drawn
    expect_message(p <- plotFrequencySpectrum(whe, usePassedOnly = TRUE),
                   "No data to plot after QC filtering")
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum handles many samples (>7)", {
    skip_if_not_installed("ggplot2")
    # Create a WHE with 8 samples to trigger the >7 code path
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(1:5, width = 1))
    freq <- matrix(seq(0.01, 0.5, length.out = 40), nrow = 5, ncol = 8)
    colnames(freq) <- paste0("S", 1:8)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(
            sample_id = paste0("S", 1:8)))
    p <- plotFrequencySpectrum(whe)
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencySpectrum returns empty ggplot on no data", {
    skip_if_not_installed("ggplot2")
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(1:3, width = 1))
    freq <- matrix(NA_real_, nrow = 3, ncol = 1)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    expect_message(
        p <- plotFrequencySpectrum(whe),
        "No data")
    expect_s3_class(p, "ggplot")
})

# ==========================================================================
# plotPairScatter
# ==========================================================================

test_that("plotPairScatter returns ggplot", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    p <- plotPairScatter(whe, pairId = "p1")
    expect_s3_class(p, "ggplot")
})

test_that("plotPairScatter with custom threshold", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    p <- plotPairScatter(whe, pairId = "p1", threshold = 0.05)
    expect_s3_class(p, "ggplot")
})

test_that("plotPairScatter correctly categorises shared/unique", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    p <- plotPairScatter(whe, pairId = "p1")
    # Build the plot data to verify categories exist
    built <- ggplot2::ggplot_build(p)
    expect_true(nrow(built$data[[2]]) > 0 ||
                nrow(built$data[[3]]) > 0)
})

# ==========================================================================
# plotQCDashboard
# ==========================================================================

test_that("plotQCDashboard returns patchwork composite", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    whe <- .make_single_whe()
    whe <- flagISNV(whe, ISNVFilter())
    p <- plotQCDashboard(whe)
    expect_s3_class(p, c("patchwork", "ggplot"))
})

test_that("plotQCDashboard with multiple samples", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    whe <- .make_pair_whe()
    whe <- flagISNV(whe, ISNVFilter(minDepth = 200L))
    p <- plotQCDashboard(whe)
    expect_s3_class(p, c("patchwork", "ggplot"))
})

test_that("plotQCDashboard without QC steps still works", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    whe <- .make_single_whe()
    # No flagISNV applied -- qcLog is empty
    p <- plotQCDashboard(whe)
    expect_s3_class(p, c("patchwork", "ggplot"))
})

test_that("plotQCDashboard without totalDepth assay still works", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(1:5, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.05, 0.12, 0.20, 0.31, 0.44),
                                       ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    whe <- flagISNV(whe, ISNVFilter())
    p <- plotQCDashboard(whe)
    expect_s3_class(p, c("patchwork", "ggplot"))
})

test_that("a recorded 0 is not drawn or counted as an iSNV (was: drawn in both plots and counted in panel d)", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    ## sample c holds only sites assessed without the allele
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.12, 0, 0.30, NA,
                                         0.15, 0.20, 0, 0,
                                         0, 0, NA, 0), ncol = 3),
                      totalDepth = matrix(1000L, 4, 3)),
        rowRanges = GenomicRanges::GRanges("seg",
            IRanges::IRanges(c(10, 20, 30, 40), width = 1)),
        colData = S4Vectors::DataFrame(sample_id = c("a", "b", "c")))
    calls <- c(0.12, 0.30, 0.15, 0.20)
    expect_equal(plotFrequencySpectrum(whe)$data$altFreq, calls)
    flagged <- flagISNV(whe, ISNVFilter(minAltReads = 0L))
    expect_equal(plotFrequencySpectrum(flagged,
                                       usePassedOnly = FALSE)$data$altFreq,
                 calls)
    dash <- plotQCDashboard(flagged)
    expect_equal(dash[[1]][[1]]$data$freq, calls)
    panel_d <- dash[[2]][[2]]
    expect_equal(panel_d$data$count, c(2L, 2L, 0L, 0L, 0L, 0L))
    expect_identical(panel_d$layers[[2]]$data$label, c("2/2", "2/2", "0/0"))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(print(dash))
})

# ==========================================================================
# plotFrequencyTrajectory
# ==========================================================================

test_that("plotFrequencyTrajectory returns ggplot", {
    skip_if_not_installed("ggplot2")
    whe <- .make_longitudinal_whe()
    p <- plotFrequencyTrajectory(whe, hostId = "H1")
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencyTrajectory with unknown host returns empty", {
    skip_if_not_installed("ggplot2")
    whe <- .make_longitudinal_whe()
    expect_message(
        p <- plotFrequencyTrajectory(whe, hostId = "NONEXISTENT"),
        "No trajectory")
    expect_s3_class(p, "ggplot")
})

test_that("plotFrequencyTrajectory errors on bad input", {
    skip_if_not_installed("ggplot2")
    whe <- .make_longitudinal_whe()
    expect_error(plotFrequencyTrajectory(whe),
                 "hostId")
})

# ==========================================================================
# plotSFS
# ==========================================================================

test_that("plotSFS returns ggplot for single SFS", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    p <- plotSFS(sfs)
    expect_s3_class(p, "ggplot")
})

test_that("plotSFS returns ggplot for normalized SFS", {
    skip_if_not_installed("ggplot2")
    whe <- .make_single_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    p <- plotSFS(sfs, normalize = TRUE)
    expect_s3_class(p, "ggplot")
})

test_that("plotSFS handles named list of SFS (overlay)", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    sfs1 <- buildSFS(whe, sampleIdx = 1L, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, sampleIdx = 2L, genomeLength = 1000L)
    p <- plotSFS(list(donor = sfs1, recipient = sfs2))
    expect_s3_class(p, "ggplot")
})

test_that("plotSFS errors on invalid input", {
    skip_if_not_installed("ggplot2")
    expect_error(plotSFS("not_sfs"), "WithinHostSFS")
})

# ==========================================================================
# Internal helpers
# ==========================================================================

test_that(".requireGgplot2 does not error when ggplot2 is installed", {
    skip_if_not_installed("ggplot2")
    expect_silent(WithinHostExperiment:::.requireGgplot2())
})

test_that(".theme_pub returns a ggplot theme", {
    skip_if_not_installed("ggplot2")
    th <- WithinHostExperiment:::.theme_pub()
    expect_s3_class(th, "theme")
})

test_that("colorblind palette has expected named elements", {
    pal <- WithinHostExperiment:::.whe_pal
    expect_true(is.list(pal))
    expect_true(all(c("blue", "vermillon", "green") %in% names(pal)))
})
