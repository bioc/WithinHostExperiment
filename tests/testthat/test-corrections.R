# tests/testthat/test-corrections.R
#
# One regression test per item in the "Corrections made before submission"
# section of NEWS.md, each named after the bug it pins, in the order of
# that section. The fixtures are built in code.

# -- fixtures ---------------------------------------------------------------

## A donor-recipient pair with read counts, on the sequences and positions
## given; donor frequencies 0.1 and 0.8, recipient 90 of 100 reads at the
## first site and 0 of 100 at the second.
.pair_with_counts <- function(seqnames, pos) {
    gr <- GenomicRanges::GRanges(seqnames, IRanges::IRanges(pos, width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("G", "T")
    freq  <- cbind(d = c(0.1, 0.8), r = c(0.9, NA))
    alt   <- cbind(d = c(10L, 80L), r = c(90L, 0L))
    depth <- cbind(d = c(100L, 100L), r = c(100L, 100L))
    WithinHostExperiment(
        assays = list(altFreq = freq, altCount = alt, totalDepth = depth),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
            role = c("donor", "recipient"), pair_id = c("p", "p"),
            row.names = c("d", "r")))
}

.folded_sfs <- function(counts, threshold = 0.03) {
    new("WithinHostSFS", counts = as.integer(counts),
        breaks = seq(threshold, 0.5, length.out = length(counts) + 1L),
        folded = TRUE, polarization = "none", sampleId = "s",
        genomeLength = 1000L, nSites = as.integer(sum(counts)),
        meanDepth = 1000, threshold = threshold, corrected = FALSE)
}

.single_whe <- function() {
    readWithinHost(system.file("extdata", "test_donor.vcf",
                               package = "WithinHostExperiment"),
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
}

# -- 1. compareSFS ----------------------------------------------------------

test_that("compareSFS reports a small expected count (was: silenced it)", {
    sfs1 <- .folded_sfs(c(3, 1, 0, 2))
    sfs2 <- .folded_sfs(c(1, 2, 1, 0))
    expect_no_warning(res <- compareSFS(sfs1, sfs2))
    expect_false(res$chisq_approx_ok)
    expect_lt(res$min_expected, 5)
    reference <- suppressWarnings(stats::chisq.test(
        rbind(c(3, 1, 0, 2), c(1, 2, 1, 0))))
    expect_equal(res$chisq_p, reference$p.value)
    big <- compareSFS(.folded_sfs(c(40, 30, 20, 10)),
                      .folded_sfs(c(35, 35, 15, 15)))
    expect_true(big$chisq_approx_ok)
})

# -- 2. flagTemporalInconsistency -------------------------------------------

test_that("transient sites are flagged per host (was: with the last host)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20),
                                                        width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("G", "T")
    ## site 1 is seen at every timepoint of H1 and once in H2
    freq <- cbind(H1_t1 = c(0.1, NA), H1_t2 = c(0.1, NA),
                  H2_t1 = c(0.1, 0.2), H2_t2 = c(NA, 0.2))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq), rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
            host_id = c("H1", "H1", "H2", "H2"), timepoint = c(1, 2, 1, 2),
            row.names = colnames(freq)))
    out <- flagTemporalInconsistency(whe, minTimepoints = 2L)
    qc <- SummarizedExperiment::assay(out, "qcPass")
    expect_true(all(qc[1, c("H1_t1", "H1_t2")]))
    expect_false(qc[1, "H2_t1"])
})

# -- 3. exactBottleneck site matching ----------------------------------------

test_that("exact_bb matches sites by chrom, pos and allele (was: pos only)", {
    ## two sites at the same position on different segments, ordered so
    ## that matching on position alone pairs the donor's 0.8 with the
    ## recipient's 90 reads instead of its 0
    two_segments <- .pair_with_counts(c("segB", "segA"), c(100, 100))
    one_segment  <- .pair_with_counts(c("seg", "seg"), c(200, 100))
    a <- exactBottleneck(two_segments, pairId = "p", maxNb = 20L,
                         method = "exact_bb")
    b <- exactBottleneck(one_segment, pairId = "p", maxNb = 20L,
                         method = "exact_bb")
    expect_identical(a$n_variants, 2L)
    expect_equal(a$loglik, b$loglik)
    expect_equal(a$Nb, b$Nb)
    ## the row order of the object does not matter either
    reversed <- one_segment[2:1, ]
    expect_equal(exactBottleneck(reversed, pairId = "p", maxNb = 20L,
                                 method = "exact_bb")$loglik, b$loglik)
})

test_that("exact_bb warns when the estimate reaches maxNb (was: silent)", {
    whe <- .pair_with_counts(c("seg", "seg"), c(100, 200))
    expect_warning(exactBottleneck(whe, pairId = "p", maxNb = 1L,
                                   method = "exact_bb"), "maxNb")
})

# -- 4. exactBottleneck documentation ---------------------------------------

test_that("wright_fisher is reported as a heuristic (was: an estimator)", {
    whe <- .pair_with_counts(c("seg", "seg"), c(100, 200))
    expect_message(exactBottleneck(whe, pairId = "p", maxNb = 10L,
                                   method = "wright_fisher"),
                   "heuristic score, not a likelihood")
    expect_no_message(exactBottleneck(whe, pairId = "p", maxNb = 10L,
                                      method = "exact_bb"))
})

# -- 5. unused arguments ----------------------------------------------------

test_that("unused arguments stop, named or not (was: absorbed by ...)", {
    whe <- .single_whe()
    expect_error(flagISNV(whe, ISNVFilter(), minDepht = 100L),
                 "Unused argument")
    expect_error(calcDiversity(whe, indices = "richness", TRUE, 5),
                 "<unnamed>")
    expect_error(calcNeutralityTests(whe, genomeLength = 1000L,
                                     usePassedonly = FALSE),
                 "Unused argument")
})

# -- 6. adaptiveFilter("ont") ----------------------------------------------

test_that("adaptiveFilter('ont') asks for 60 reads (was: 30)", {
    expect_identical(slot(adaptiveFilter("ont"), "minDepth"), 60L)
})

# -- 7. documented as implemented -------------------------------------------

test_that("truncation rescales the whole spectrum by max/(max - thr)", {
    sfs <- .folded_sfs(c(10, 20, 30), threshold = 0.03)
    corrected <- correctSFSBias(sfs, method = "truncation")
    expect_identical(corrected@counts,
                     as.integer(round(c(10, 20, 30) / ((0.5 - 0.03) / 0.5))))
    expect_true(corrected@corrected)
})

test_that("SFS bins are right-closed, the lowest break included", {
    breaks <- seq(0.03, 0.5, length.out = 3)
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:4 * 10, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.03, breaks[2], 0.5, 0.3),
                                       ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    sfs <- buildSFS(whe, nBins = 2L, threshold = 0.03, genomeLength = 1000L)
    ## 0.03 and the middle break fall in the first bin, 0.3 and 0.5 in
    ## the second
    expect_identical(sfs@counts, c(2L, 2L))
})

test_that("the SFS depth slot holds the median (was: documented as mean)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3 * 10, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1),
                      totalDepth = matrix(c(100L, 1000L, 5000L), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_equal(sfs@meanDepth, 1000)
    expect_output(show(sfs), "median depth")
})

test_that("the consensus object holds frequencies and qcPass only", {
    mk <- function(freq) {
        gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20),
                                                            width = 1))
        S4Vectors::mcols(gr)$ref <- c("A", "C")
        S4Vectors::mcols(gr)$alt <- c("G", "T")
        WithinHostExperiment(
            assays = list(altFreq = matrix(freq, ncol = 1,
                                           dimnames = list(NULL, "s")),
                          totalDepth = matrix(1000L, 2, 1,
                                              dimnames = list(NULL, "s"))),
            rowRanges = gr,
            colData = S4Vectors::DataFrame(sample_id = "s"))
    }
    cons <- consensusISNV(list(a = mk(c(0.10, 0.20)), b = mk(c(0.12, 0.21))))
    expect_identical(SummarizedExperiment::assayNames(cons),
                     c("altFreq", "qcPass"))
})

test_that("the plausibility ranges are the package's own heuristics", {
    res <- checkDiversityPlausibility(
        S4Vectors::DataFrame(pi = 1e-4, shannon = 1), pathogen = "tb")
    expect_equal(res$expected_min, c(1e-9, 0))
    expect_equal(res$expected_max, c(1e-4, 5))
})

# -- 8. ISNVFilter ------------------------------------------------------------

test_that("flagISNV warns on replicateConc and the constructor needs altFreq", {
    whe <- .single_whe()
    expect_warning(flagISNV(whe, ISNVFilter(replicateConc = TRUE)),
                   "flagReplicateDiscordance")
    expect_no_warning(flagISNV(whe, ISNVFilter()))
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(10, width = 1))
    expect_error(WithinHostExperiment(
        assays = list(totalDepth = matrix(100L, ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1")), "altFreq")
})

# -- 9. trackFrequency --------------------------------------------------------

test_that("trackFrequency keeps column-major order (was: a nested loop)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20),
                                                        width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("G", "T")
    freq <- cbind(t1 = c(0.1, NA), t2 = c(0.2, 0.3))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq), rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("t1", "t2"),
            host_id = "H", timepoint = c(1, 2), row.names = c("t1", "t2")))
    traj <- trackFrequency(whe)
    expect_identical(traj$frequency, c(0.1, 0.2, 0.3))
    expect_identical(traj$timepoint, c(1, 2, 2))
})

# -- 10. plots ----------------------------------------------------------------

test_that("plots clip the view instead of dropping points (was: a warning)", {
    skip_if_not_installed("ggplot2")
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    pair <- readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role = c("donor", "recipient"), pair_id = c("p1", "p1")),
        caller = "ivar")
    expect_no_warning(built <- ggplot2::ggplot_build(
        plotPairScatter(pair, pairId = "p1")))
    n_sites <- nrow(exportPairFrequencies(pair, pairId = "p1"))
    ## every site is still drawn: count the rows of the point layers
    is_point <- vapply(built$plot$layers,
                       function(l) inherits(l$geom, "GeomPoint"), logical(1))
    n_points <- sum(vapply(built$data[is_point], nrow, integer(1)))
    expect_identical(n_points, n_sites)
    expect_no_warning(ggplot2::ggplot_build(plotFrequencySpectrum(pair)))
})

test_that("trajectory and SFS plots use the Okabe-Ito palette, recycled", {
    skip_if_not_installed("ggplot2")
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:9 * 10, width = 1))
    S4Vectors::mcols(gr)$ref <- rep("A", 9)
    S4Vectors::mcols(gr)$alt <- rep("G", 9)
    freq <- cbind(t1 = seq(0.1, 0.5, length.out = 9),
                  t2 = seq(0.2, 0.6, length.out = 9))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq), rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("t1", "t2"),
            host_id = "H", timepoint = c(1, 2), row.names = c("t1", "t2")))
    okabe_ito <- unlist(WithinHostExperiment:::.whe_pal, use.names = FALSE)
    p <- plotFrequencyTrajectory(whe, hostId = "H")
    cols <- unique(ggplot2::ggplot_build(p)$data[[1]]$colour)
    expect_true(all(cols %in% okabe_ito))
    ## nine variants need the eight colours recycled
    expect_length(cols, 8L)
    sfs <- lapply(1:9, function(i) {
        s <- .folded_sfs(c(3, 2, 1))
        s@sampleId <- paste0("s", i)
        s
    })
    names(sfs) <- vapply(sfs, function(s) s@sampleId, "")
    fills <- unique(ggplot2::ggplot_build(plotSFS(sfs))$data[[1]]$fill)
    expect_true(all(fills %in% okabe_ito))
})

# -- README claim: QC only removes calls, so pi cannot rise ------------------

test_that("QC cannot raise pi, as the README states", {
    vcfs <- system.file("extdata",
        c("test_donor.vcf", "test_recipient.vcf", "test_independent.vcf"),
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcfs,
        colData = S4Vectors::DataFrame(sample_id = c("d", "r", "i")),
        caller = "ivar",
        genome = Biostrings::readDNAStringSet(system.file("extdata",
            "test_ref.fa", package = "WithinHostExperiment")))
    before <- calcDiversity(whe, indices = "pi", usePassedOnly = FALSE)$pi
    for (filt in list(ISNVFilter(), ISNVFilter(minDepth = 500L),
                      ISNVFilter(minFreq = 0.10),
                      ISNVFilter(maxFreq = 0.20),
                      ISNVFilter(minAltReads = 300L))) {
        after <- calcDiversity(flagISNV(whe, filt), indices = "pi")$pi
        expect_true(all(after <= before), info = format(slot(filt, "minFreq")))
    }
    ## and a filter that removes something lowers it for that sample
    strict <- calcDiversity(flagISNV(whe, ISNVFilter(minFreq = 0.30)),
                            indices = c("pi", "richness"))
    expect_true(any(strict$pi < before))
})
