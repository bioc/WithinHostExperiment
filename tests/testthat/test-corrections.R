# tests/testthat/test-corrections.R
#
# One regression test per item in the "Corrections made before submission"
# section of NEWS.md, each named after the bug it pins, in the order of
# that section. The fixtures are built in code.

# -- fixtures ---------------------------------------------------------------

## A donor-recipient pair with read counts, on the sequences and positions
## given; donor frequencies 0.1 and 0.8, recipient 90 of 100 reads at the
## first site and 0 of 100 at the second.
.pair_with_counts <- function(seqnames, pos, ref = c("A", "C"),
                              alt = c("G", "T")) {
    gr <- GenomicRanges::GRanges(seqnames, IRanges::IRanges(pos, width = 1))
    S4Vectors::mcols(gr)$ref <- ref
    S4Vectors::mcols(gr)$alt <- alt
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

## An object on one 1000-base segment from a frequency matrix with one
## column per sample, with read counts and depths to match and every call
## passing so far.
.whe_from_freq <- function(freq, depth, pos, colData) {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(pos, width = 1))
    GenomeInfoDb::seqlengths(gr) <- c(seg = 1000L)
    S4Vectors::mcols(gr)$ref <- rep("A", length(pos))
    S4Vectors::mcols(gr)$alt <- rep("G", length(pos))
    dimnames(depth) <- dimnames(freq)
    alt <- matrix(as.integer(round(freq * depth)), nrow(freq),
                  dimnames = dimnames(freq))
    WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth, altCount = alt,
                      qcPass = matrix(TRUE, nrow(freq), ncol(freq),
                                      dimnames = dimnames(freq))),
        rowRanges = gr, colData = colData)
}

## Six calls in one sample, all passing ISNVFilter().
.single_whe <- function() {
    freq <- matrix(c(0.35, 0.08, 0.12, 0.20, 0.05, 0.30), ncol = 1,
                   dimnames = list(NULL, "donor"))
    depth <- matrix(c(1500L, 800L, 1200L, 2000L, 900L, 1100L), ncol = 1)
    .whe_from_freq(freq, depth, pos = c(45, 112, 230, 401, 655, 889),
        colData = S4Vectors::DataFrame(sample_id = "donor",
                                       row.names = "donor"))
}

## A donor-recipient pair whose frequencies reach both ends of the range,
## where jitter or a histogram bin falls outside the axis.
.edge_pair <- function() {
    freq <- cbind(
        donor = c(0.01, 0.02, 0.99, 0.97, 0.50, 0.30, 0.04, 0.60, NA, NA,
                  0.25, 0.015),
        recipient = c(0.99, NA, 0.98, 0.03, 0.02, 0.31, NA, 0.01, 0.40,
                      0.97, 0.26, 0.99))
    depth <- matrix(1000L, nrow(freq), 2L)
    .whe_from_freq(freq, depth, pos = seq_len(nrow(freq)) * 50L,
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role = c("donor", "recipient"), pair_id = c("p1", "p1"),
            row.names = c("donor", "recipient")))
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

test_that("exact_bb matches sites by allele at a multi-allelic position (was: pos only)", {
    ## two alternative alleles at one position: G at 0.1 in the donor and
    ## 90 of 100 recipient reads, T at 0.8 and 0 reads. Matching on
    ## chromosome and position alone cannot tell the rows apart.
    one_position <- .pair_with_counts(c("seg", "seg"), c(100, 100),
                                      ref = c("A", "A"), alt = c("G", "T"))
    two_positions <- .pair_with_counts(c("seg", "seg"), c(100, 200),
                                       ref = c("A", "A"), alt = c("G", "T"))
    a <- exactBottleneck(one_position, pairId = "p", maxNb = 20L,
                         method = "exact_bb")
    b <- exactBottleneck(two_positions, pairId = "p", maxNb = 20L,
                         method = "exact_bb")
    expect_identical(a$n_variants, 2L)
    expect_equal(a$loglik, b$loglik)
    expect_equal(a$Nb, b$Nb)
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
    expect_error(calcNeutralityTests(whe, usePassedonly = FALSE),
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

test_that("replicate intersection flags only sites one replicate lacks, in every replicate", {
    ## the documented recipe: requireBothDetected = TRUE with an infinite
    ## tolerance. Site 1 differs by 14 points between the replicates,
    ## site 2 is seen in one replicate only, site 3 agrees; 'other' is
    ## in no replicate group.
    freq <- cbind(r1 = c(0.37, 0.04, 0.10), r2 = c(0.23, NA, 0.10),
                  other = c(0.37, 0.04, 0.10))
    whe <- .whe_from_freq(freq, matrix(1000L, 3, 3), pos = c(10, 20, 30),
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
            replicate_group = c("g", "g", NA), replicate_id = c("1", "2", NA),
            row.names = colnames(freq)))
    q <- SummarizedExperiment::assay(flagReplicateDiscordance(whe,
        freqTolerance = Inf, requireBothDetected = TRUE), "qcPass")
    expect_identical(unname(q[, "r1"]), c(TRUE, FALSE, TRUE))
    expect_identical(unname(q[, "r2"]), c(TRUE, FALSE, TRUE))
    expect_true(all(q[, "other"]))
    ## a finite tolerance flags the frequency difference too, in both
    q5 <- SummarizedExperiment::assay(flagReplicateDiscordance(whe,
        freqTolerance = 0.05), "qcPass")
    expect_identical(unname(q5[, "r1"]), c(FALSE, FALSE, TRUE))
    expect_identical(unname(q5[, "r2"]), c(FALSE, FALSE, TRUE))
})

test_that("the plausibility ranges are the package's own heuristics", {
    res <- checkDiversityPlausibility(
        S4Vectors::DataFrame(pi = 1e-4, shannon = 1), pathogen = "tb")
    expect_equal(res$expected_min, c(1e-9, 0))
    expect_equal(res$expected_max, c(1e-4, 5))
})

# -- 8. ISNVFilter ------------------------------------------------------------

test_that("flagISNV applies replicateConc (was: recorded, warned, ignored)", {
    ## r1 and r2 are replicates: site 1 differs by 14 points, site 2 is
    ## called in r1 only, site 3 agrees, and site 4 differs and fails
    ## minFreq in r2; 'other' has no replicate
    freq <- cbind(r1 = c(0.37, 0.04, 0.10, 0.20),
                  r2 = c(0.23, NA, 0.10, 0.02),
                  other = c(0.37, 0.04, 0.10, 0.20))
    cd <- S4Vectors::DataFrame(sample_id = colnames(freq),
        replicate_group = c("g", "g", NA), replicate_id = c("1", "2", NA),
        row.names = colnames(freq))
    whe <- .whe_from_freq(freq, matrix(1000L, 4, 3),
                          pos = c(10, 20, 30, 40), colData = cd)
    rule <- ConcordanceRule(freqTolerance = 0.05)
    expect_warning(got <- flagISNV(whe, ISNVFilter(replicateConc = rule)),
                   "1 sample\\(s\\) without a replicate \\(other\\)")
    q <- SummarizedExperiment::assay(got, "qcPass")
    expect_identical(unname(q[, "r1"]), c(FALSE, FALSE, TRUE, FALSE))
    expect_identical(unname(q[, "r2"]), c(FALSE, FALSE, TRUE, FALSE))
    expect_true(all(q[, "other"]))
    ## the same flags and log rows as the separate function after the
    ## per-sample thresholds; the rule travels with the filter
    want <- flagReplicateDiscordance(flagISNV(whe, ISNVFilter()),
                                     freqTolerance = 0.05)
    expect_identical(q, SummarizedExperiment::assay(want, "qcPass"))
    cols <- c("step", "parameter", "value", "n_flagged", "n_passed")
    expect_identical(as.data.frame(qcLog(got))[, cols],
                     as.data.frame(qcLog(want))[, cols])
    expect_identical(tail(qcLog(got)$n_flagged, 2), c(5L, 0L))
    expect_identical(slot(qcFilters(got)[[1]], "replicateConc"), rule)
    ## a rule with nothing to compare is an error, not a silent no-op
    expect_error(flagISNV(.single_whe(), ISNVFilter(replicateConc = rule)),
                 "replicate_group")
    lone <- whe
    SummarizedExperiment::colData(lone)$replicate_group <- c("a", "b", NA)
    expect_error(flagISNV(lone, ISNVFilter(replicateConc = rule)),
                 "two or more samples")
})

test_that("a difference equal to freqTolerance passes (was: flagged by rounding)", {
    ## sites 1 and 2 are 50 reads in 1000 apart, which the subtraction
    ## puts just above 0.05; site 3 is 51 reads apart; site 4 is 0.375
    ## against 0.125, exactly 0.25 apart in binary
    freq <- cbind(r1 = c(150, 350, 100, 375) / 1000,
                  r2 = c(200, 400, 151, 125) / 1000)
    whe <- .whe_from_freq(freq, matrix(1000L, 4, 2),
                          pos = c(10, 20, 30, 40),
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
            replicate_group = c("g", "g"), row.names = colnames(freq)))
    qc <- function(x) unname(SummarizedExperiment::assay(x, "qcPass"))
    both <- function(v) matrix(v, length(v), 2L)
    by_fun <- function(tol, rbd = TRUE)
        qc(flagReplicateDiscordance(whe, freqTolerance = tol,
                                    requireBothDetected = rbd))
    by_rule <- function(tol, rbd = TRUE)
        qc(flagISNV(whe, ISNVFilter(replicateConc =
                                        ConcordanceRule(tol, rbd))))
    expect_identical(by_fun(0.05), both(c(TRUE, TRUE, FALSE, FALSE)))
    expect_identical(by_rule(0.05), both(c(TRUE, TRUE, FALSE, FALSE)))
    ## the binary-exact boundary from both sides
    expect_identical(by_fun(0.25), both(rep(TRUE, 4)))
    expect_identical(by_rule(0.25, FALSE), both(rep(TRUE, 4)))
    expect_identical(by_fun(0.125), both(c(TRUE, TRUE, TRUE, FALSE)))
    expect_identical(by_rule(0.125), both(c(TRUE, TRUE, TRUE, FALSE)))
    expect_identical(
        qcLog(flagReplicateDiscordance(whe, freqTolerance = 0.05))$n_flagged,
        c(4L, 0L))
})

test_that("flagISNV() refuses a filter changed into an invalid one (was: applied it)", {
    ## identical replicates: a negative tolerance flagged every call, and a
    ## rule that flags nothing wrote rows that look like a check
    freq <- cbind(r1 = c(0.10, 0.20), r2 = c(0.10, 0.20))
    whe <- .whe_from_freq(freq, matrix(1000L, 2, 2), pos = c(10, 20),
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
            replicate_group = c("g", "g"), row.names = colnames(freq)))
    f <- ISNVFilter(replicateConc = ConcordanceRule())
    f@replicateConc@freqTolerance <- -1
    expect_error(flagISNV(whe, f), "'freqTolerance' must be positive")
    f <- ISNVFilter(replicateConc = ConcordanceRule(Inf, TRUE))
    f@replicateConc@requireBothDetected <- FALSE
    expect_error(flagISNV(whe, f), "flags nothing")
    f <- ISNVFilter()
    f@minDepth <- -5L
    expect_error(flagISNV(whe, f), "'minDepth' must be non-negative")
})

test_that("validObject() reports an invalid filter in qcFilters (was: TRUE)", {
    ## readRDS() restores a filter saved under 0.99.0 with the logical
    ## replicateConc it then held; show() and flagISNV() failed on it with
    ## a slot-access error while the object validated
    old <- ISNVFilter()
    attr(old, "replicateConc") <- FALSE
    x <- .single_whe()
    qcFilters(x) <- list(ISNVFilter(), old)
    expect_error(validObject(x),
                 "'qcFilters' element 2 is not a valid ISNVFilter")
    expect_error(validObject(x), "rebuild it with ISNVFilter\\(\\)")
    expect_error(flagISNV(.single_whe(), old), "replicateConc")
    ## a stored rule changed after construction is reported too
    f <- ISNVFilter(replicateConc = ConcordanceRule())
    f@replicateConc@freqTolerance <- -1
    qcFilters(x) <- list(f)
    expect_error(validObject(x), "'freqTolerance' must be positive")
    ## valid filters, and the message for anything else, are unchanged
    qcFilters(x) <- list(ISNVFilter(), f@replicateConc)
    expect_error(validObject(x),
                 "'qcFilters' must contain only ISNVFilter objects")
    qcFilters(x) <- list(ISNVFilter())
    expect_true(validObject(x))
})

test_that("the constructor needs altFreq", {
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
    pair <- .edge_pair()
    p <- plotPairScatter(pair, pairId = "p1")
    s <- plotFrequencySpectrum(pair)
    ## the limits sit on the coordinate system, which clips the view, and
    ## not on the scales, which would remove what falls outside
    for (g in list(p, s)) {
        expect_false(is.null(g$coordinates$limits$x))
        expect_null(g$scales$get_scales("x")$limits)
    }
    expect_false(is.null(p$coordinates$limits$y))
    expect_null(p$scales$get_scales("y")$limits)
    ## ggplot2 removes rows outside a scale when the plot is drawn, so the
    ## whole table is built, on a device that writes nothing
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_warning(ggplot2::ggplotGrob(p))
    expect_no_warning(ggplot2::ggplotGrob(s))
    ## every site is still drawn: count the point rows with a position
    built <- ggplot2::ggplot_build(p)
    is_point <- vapply(built$plot$layers,
                       function(l) inherits(l$geom, "GeomPoint"), logical(1))
    n_points <- sum(vapply(built$data[is_point], function(d)
        sum(is.finite(d$x) & is.finite(d$y)), integer(1)))
    expect_identical(n_points,
                     nrow(exportPairFrequencies(pair, pairId = "p1")))
})

test_that("plotQCDashboard clips the frequency panel instead of dropping calls (was: a warning)", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    pair <- .edge_pair()
    ## the layout is (a + b) / (c + d): panel a is the first of the first row
    panel_a <- plotQCDashboard(pair)[[1]][[1]]
    expect_true("freq" %in% names(panel_a$data))
    ## the limits sit on the coordinate system, which clips the view, and
    ## not on the scale, which would remove the edge bars when the panel
    ## is drawn
    expect_false(is.null(panel_a$coordinates$limits$x))
    expect_null(panel_a$scales$get_scales("x")$limits)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_warning(ggplot2::ggplotGrob(panel_a))
    ## every call is counted in some bar that has a position
    built <- ggplot2::ggplot_build(panel_a)
    bars <- built$data[[1]]
    n_calls <- sum(!is.na(SummarizedExperiment::assay(pair, "altFreq")))
    expect_equal(sum(bars$count[is.finite(bars$xmin) & is.finite(bars$xmax)]),
                 n_calls)
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
    ## twelve sites in three samples on a 1000-base segment; s1 and s2 are
    ## replicates of one sample, so the replicate step applies to them
    freq <- cbind(
        s1 = c(0.02, 0.05, 0.08, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
               0.45, 0.04),
        s2 = c(0.03, 0.05, 0.15, 0.12, 0.30, 0.20, NA, 0.30, 0.36, 0.40,
               0.44, NA),
        s3 = c(NA, 0.06, 0.09, NA, 0.18, 0.22, 0.27, NA, 0.33, 0.41, NA, 0.05))
    depth <- matrix(rep(c(60L, 150L, 400L, 800L, 1500L, 3000L),
                        length.out = 36L), 12L, 3L)
    whe <- .whe_from_freq(freq, depth, pos = seq(40, 920, by = 80),
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
            replicate_group = c("g", "g", NA), replicate_id = c("1", "2", NA),
            row.names = colnames(freq)))
    before <- calcDiversity(whe, indices = "pi", usePassedOnly = FALSE)$pi
    for (filt in list(ISNVFilter(), ISNVFilter(minDepth = 500L),
                      ISNVFilter(minFreq = 0.10),
                      ISNVFilter(maxFreq = 0.20),
                      ISNVFilter(minAltReads = 300L))) {
        after <- calcDiversity(flagISNV(whe, filt), indices = "pi")$pi
        expect_true(all(after <= before), info = format(slot(filt, "minFreq")))
    }
    ## the replicate step, which is what the README describes, at every
    ## tolerance
    lowered <- FALSE
    for (tol in c(0.01, 0.03, 0.05, 0.2)) {
        after <- calcDiversity(flagReplicateDiscordance(whe,
            freqTolerance = tol), indices = "pi")$pi
        expect_true(all(after <= before), info = paste("tolerance", tol))
        lowered <- lowered || any(after < before)
    }
    expect_true(lowered)
    ## and a filter that removes something lowers it for that sample
    strict <- calcDiversity(flagISNV(whe, ISNVFilter(minFreq = 0.30)),
                            indices = c("pi", "richness"))
    expect_true(any(strict$pi < before))
})
