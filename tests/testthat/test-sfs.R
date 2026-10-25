# tests/testthat/test-sfs.R

library(GenomicRanges)
library(S4Vectors)

## Helper: fixed frequencies and depths, so every run builds the same
## spectrum
.make_whe <- function(n = 10) {
    gr <- GRanges("seg", IRanges::IRanges(seq_len(n) * 100L, width = 1))
    freq <- matrix(seq(0.04, 0.44, length.out = n), ncol = 1)
    depth <- matrix(as.integer(seq(500, 3000, length.out = n)), ncol = 1)
    WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth),
        rowRanges = gr,
        colData = DataFrame(sample_id = "S1"))
}

# ---- buildSFS ----

test_that("buildSFS creates a valid WithinHostSFS", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_s4_class(sfs, "WithinHostSFS")
    expect_true(validObject(sfs))
    expect_equal(sfs@sampleId, "S1")
    expect_equal(sfs@genomeLength, 1000L)
    expect_true(sfs@folded)
    expect_false(sfs@corrected)
    expect_true(sum(sfs@counts) > 0)
})

test_that("buildSFS fold=FALSE gives unfolded SFS", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, fold = FALSE, genomeLength = 1000L)
    expect_false(sfs@folded)
    expect_true(max(sfs@breaks) > 0.5)
})

test_that("buildSFS respects threshold", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, threshold = 0.10, genomeLength = 1000L)
    expect_equal(sfs@threshold, 0.10)
    expect_true(min(sfs@breaks) >= 0.10)
})

.two_sample_whe <- function(a, b) {
    gr <- GRanges("seg", IRanges::IRanges(seq_along(a) * 100L, width = 1))
    freq <- cbind(a = a, b = b)
    WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(sample_id = c("a", "b"),
                            row.names = c("a", "b")))
}

test_that("buildSFS() folds a sample with no call to an empty spectrum (was: 'x' must be numeric)", {
    whe <- .two_sample_whe(a = c(0.1, 0.7), b = c(NA, NA))
    sfs <- buildSFS(whe, sampleIdx = 2L, nBins = 5L)
    expect_identical(sfs@counts, integer(5))
    expect_identical(sfs@nSites, 0L)
    ## buildSFSList() no longer fails for the whole object
    sfs_list <- buildSFSList(whe, nBins = 5L)
    expect_length(sfs_list, 2L)
    expect_identical(sum(sfs_list$b@counts), 0L)
    ## the fold itself is unchanged: 0.7 is binned as 0.3
    expect_identical(sfs_list$a@counts, c(1L, 0L, 1L, 0L, 0L))
})

test_that("buildSFS() names the argument it cannot use (was: errors from hist() or a pooled spectrum)", {
    data(example_whe)
    expect_error(buildSFS(example_whe, sampleIdx = "nope"),
                 "'sampleIdx'.*nope.*donor, recipient, independent")
    ## two samples used to be pooled into one spectrum without a word
    expect_error(buildSFS(example_whe, sampleIdx = 1:2), "'sampleIdx'")
    expect_error(buildSFS(example_whe, sampleIdx = 4L), "'sampleIdx'")
    expect_error(buildSFS(example_whe, sampleIdx = 1.5), "'sampleIdx'")
    expect_error(buildSFS(example_whe, sampleIdx = NA), "'sampleIdx'")
    expect_identical(buildSFS(example_whe, sampleIdx = "recipient"),
                     buildSFS(example_whe, sampleIdx = 2L))
    expect_identical(buildSFS(example_whe, sampleIdx = 2),
                     buildSFS(example_whe, sampleIdx = 2L))
    expect_error(buildSFS(example_whe, nBins = 0L), "'nBins'")
    expect_error(buildSFS(example_whe, nBins = 2.5), "'nBins'")
    expect_error(buildSFS(example_whe, nBins = c(5L, 10L)), "'nBins'")
    expect_error(buildSFS(example_whe, threshold = -1), "'threshold'")
    expect_error(buildSFS(example_whe, threshold = 0.7),
                 "'threshold'.*0.5")
    expect_error(buildSFS(example_whe, threshold = 1, fold = FALSE),
                 "'threshold'")
    expect_error(buildSFS(example_whe, threshold = NA_real_), "'threshold'")
    expect_s4_class(buildSFS(example_whe, threshold = 0.7, fold = FALSE),
                    "WithinHostSFS")
    expect_s4_class(buildSFS(example_whe, threshold = 0, nBins = 10),
                    "WithinHostSFS")
})

test_that("buildSFS() leaves out sites at frequency 0 or 1 (was: a fixed site in the top unfolded bin)", {
    whe <- .two_sample_whe(a = c(0.1, 0.6, 1, 1, 0),
                           b = c(0.2, NA, NA, NA, NA))
    unf <- buildSFS(whe, sampleIdx = "a", nBins = 5L, fold = FALSE)
    expect_identical(unf@nSites, 2L)
    expect_identical(unf@counts, c(1L, 0L, 1L, 0L, 0L))
    ## from a threshold of 0, neither a 0 nor a 1 folded to 0 is counted
    fol <- buildSFS(whe, sampleIdx = "a", nBins = 5L, threshold = 0)
    expect_identical(fol@nSites, 2L)
})

test_that("buildSFSList returns one SFS per sample", {
    data(example_whe)
    sfs_list <- buildSFSList(example_whe, genomeLength = 1000L)
    expect_equal(length(sfs_list), ncol(example_whe))
    expect_true(all(vapply(sfs_list, is, logical(1), "WithinHostSFS")))
})

# ---- correctSFSBias ----

test_that("correctSFSBias truncation method works", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    sfs_c <- correctSFSBias(sfs, method = "truncation")
    expect_true(sfs_c@corrected)
    expect_true(sum(sfs_c@counts) >= sum(sfs@counts))
})

test_that("correctSFSBias binomial method works", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    sfs_c <- correctSFSBias(sfs, method = "binomial")
    expect_true(sfs_c@corrected)
})

test_that("correctSFSBias() corrects a spectrum of small bins (was: every bin of 7 or fewer left as it was)", {
    ## the help page's example: 15 iSNVs, no bin above 2
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf, colData = DataFrame(sample_id = "donor"),
                          caller = "ivar")
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_identical(sfs@nSites, 15L)
    tr <- correctSFSBias(sfs, method = "truncation")
    ## 15 * 0.5 / (0.5 - 0.03) = 15.96
    expect_identical(tr@nSites, 16L)
    expect_identical(sum(tr@counts), tr@nSites)
    expect_true(all(tr@counts >= sfs@counts))
    ## one iSNV in each of 20 bins at a depth of 60
    flat <- new("WithinHostSFS", counts = rep(1L, 20L),
                breaks = seq(0.03, 0.5, length.out = 21L), folded = TRUE,
                polarization = "none", sampleId = "s", nSites = 20L,
                genomeLength = 1000L, meanDepth = 60, threshold = 0.03)
    mids <- (flat@breaks[-21L] + flat@breaks[-1L]) / 2
    p_detect <- pmax(1 - pbinom(ceiling(60 * 0.03) - 1, 60, mids), 0.01)
    bi <- correctSFSBias(flat, method = "binomial")
    expect_identical(bi@nSites, as.integer(round(sum(1 / p_detect))))
    expect_identical(sum(bi@counts), bi@nSites)
    expect_identical(correctSFSBias(flat)@nSites,
                     as.integer(round(20 * 0.5 / 0.47)))
})

# ---- neutralityFromSFS ----

test_that("neutralityFromSFS returns all statistics", {
    whe <- .make_whe(20)
    sfs <- buildSFS(whe, fold = FALSE, genomeLength = 1000L)
    res <- neutralityFromSFS(sfs)
    expect_true(is.list(res))
    expect_true(all(c("S", "n", "pi", "theta_W", "theta_H",
                       "tajimaD", "fusFs", "fayWuH") %in% names(res)))
    expect_true(is.numeric(res$tajimaD))
    expect_true(is.numeric(res$fayWuH))
})

test_that("neutralityFromSFS handles empty SFS", {
    sfs <- new("WithinHostSFS",
        counts = integer(20), breaks = seq(0.03, 1, length.out = 21),
        folded = FALSE, polarization = "reference", sampleId = "empty",
        genomeLength = 1000L, nSites = 0L, meanDepth = 1000,
        threshold = 0.03, corrected = FALSE)
    expect_warning(expect_warning(res <- neutralityFromSFS(sfs),
                                  "at least three segregating sites"),
                   "no segregating site")
    expect_equal(res$S, 0L)
    expect_true(is.na(res$tajimaD))
    ## Fu's Fs is undefined without a segregating site (was: 0)
    expect_identical(res$fusFs, NA_real_)
})

test_that("neutralityFromSFS() puts pi and theta_H on the same classes, so H is 0 under neutrality (was: biased by about -theta/n)", {
    ## one unfolded bin centred on 0.5 with four sites, n = 10
    one <- new("WithinHostSFS", counts = 4L, breaks = c(0.45, 0.55),
               folded = FALSE, polarization = "reference", sampleId = "s",
               genomeLength = 1000L, nSites = 4L, meanDepth = 10,
               threshold = 0.45, corrected = FALSE)
    res <- neutralityFromSFS(one, n = 10)
    expect_equal(res$pi, 2 * 5 * 5 * 4 / (10 * 9) / 1000)
    expect_equal(res$theta_H, 2 * 5^2 * 4 / (10 * 9) / 1000)
    expect_equal(res$fayWuH, 0)
    ## the standard neutral spectrum, xi_i proportional to 1 / i, on bins
    ## centred at i / n: pi, theta_W and theta_H all estimate theta
    n <- 20L
    i <- seq_len(n - 1L)
    xi <- as.integer(232792560 / i)
    sfs <- new("WithinHostSFS", counts = xi, breaks = (c(i, n) - 0.5) / n,
               folded = FALSE, polarization = "reference", sampleId = "s",
               genomeLength = 1000000000L, nSites = sum(xi),
               meanDepth = n, threshold = 0.5 / n, corrected = FALSE)
    ## so many sites exceed the haplotypes n sequences can carry
    expect_warning(res <- neutralityFromSFS(sfs, n = n),
                   "Fu's Fs is undefined")
    expect_equal(res$pi, res$theta_H)
    expect_equal(res$pi, res$theta_W)
    expect_lt(abs(res$fayWuH), 1e-12)
    expect_equal(res$tajimaD, 0)
})

.sweep_whe <- function() {
    ## two rare calls and five high-frequency derived ones
    freq <- c(0.05, 0.08, 0.85, 0.9, 0.92, 0.95, 0.97)
    gr <- GRanges("seg", IRanges::IRanges(seq_along(freq) * 100L, width = 1))
    WithinHostExperiment(
        assays = list(altFreq = matrix(freq, ncol = 1),
                      totalDepth = matrix(1000L, nrow = length(freq),
                                          ncol = 1)),
        rowRanges = gr,
        colData = DataFrame(sample_id = "S1"))
}

test_that("neutralityFromSFS() refuses an unfolded spectrum QC capped below 1 (was: H positive by construction)", {
    whe <- .sweep_whe()
    ## without QC nothing is capped, and the sweep gives a negative H
    sfs <- buildSFS(whe, fold = FALSE, nBins = 10L, genomeLength = 1000L)
    expect_identical(sfs@maxFreq, 1)
    expect_lt(neutralityFromSFS(sfs)$fayWuH, 0)
    ## the default filter removes every call above 0.5, the classes H
    ## weights most
    capped <- flagISNV(whe, ISNVFilter(minDepth = 100L))
    sfs <- buildSFS(capped, fold = FALSE, nBins = 10L, genomeLength = 1000L)
    expect_identical(sfs@maxFreq, 0.5)
    expect_error(neutralityFromSFS(sfs), "removed by QC")
    expect_error(neutralityFromSFS(sfs), "maxFreq = 1", fixed = TRUE)
    ## a filter that keeps them gives the negative H back
    kept <- flagISNV(whe, ISNVFilter(minDepth = 100L, maxFreq = 1))
    sfs <- buildSFS(kept, fold = FALSE, nBins = 10L, genomeLength = 1000L)
    expect_identical(sfs@maxFreq, 1)
    expect_lt(neutralityFromSFS(sfs)$fayWuH, 0)
    ## the strictest of several filters is recorded, and only when QC is
    ## applied
    twice <- flagISNV(capped, ISNVFilter(minDepth = 100L, maxFreq = 0.4))
    expect_identical(buildSFS(twice, genomeLength = 1000L)@maxFreq, 0.4)
    expect_identical(buildSFS(twice, usePassedOnly = FALSE,
                              genomeLength = 1000L)@maxFreq, 1)
    ## flagISNV() adds to earlier flags, so flagging the capped object
    ## again with maxFreq = 1 does not restore the removed calls, and the
    ## message says to start from the unfiltered object
    again <- flagISNV(capped, ISNVFilter(minDepth = 100L, maxFreq = 1))
    sfs <- buildSFS(again, fold = FALSE, nBins = 10L, genomeLength = 1000L)
    expect_identical(sfs@maxFreq, 0.5)
    expect_error(neutralityFromSFS(sfs), "Flag the unfiltered object")
})

test_that("buildSFS() names fold and a cap of 0 (was: a base-R or validity error)", {
    whe <- .make_whe(20)
    expect_error(buildSFS(whe, fold = NA, genomeLength = 1000L),
                 "'fold' must be TRUE or FALSE")
    none <- flagISNV(whe, ISNVFilter(minFreq = 0, maxFreq = 0))
    expect_error(buildSFS(none, genomeLength = 1000L),
                 "QC removed every call")
})

test_that("WithinHostSFS carries its maxFreq through validity, show, correction and comparison", {
    expect_error(new("WithinHostSFS", maxFreq = 0), "maxFreq")
    expect_error(new("WithinHostSFS", maxFreq = 1.5), "maxFreq")
    expect_error(new("WithinHostSFS", maxFreq = NA_real_), "maxFreq")
    expect_identical(new("WithinHostSFS")@maxFreq, 1)
    whe <- .make_whe(20)
    open <- buildSFS(whe, genomeLength = 1000L)
    capped <- buildSFS(flagISNV(whe, ISNVFilter(minDepth = 100L)),
                       genomeLength = 1000L)
    expect_output(show(capped), "maxFreq: 0.5")
    expect_false(any(grepl("maxFreq", capture.output(show(open)))))
    for (m in c("truncation", "binomial"))
        expect_identical(correctSFSBias(capped, method = m)@maxFreq, 0.5)
    expect_error(compareSFS(open, capped), "sfs1: 1; sfs2: 0.5")
    expect_no_error(compareSFS(capped, capped))
})

# ---- compareSFS ----

test_that("compareSFS compares two SFS objects", {
    whe <- .make_whe(20)
    sfs1 <- buildSFS(whe, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, genomeLength = 1000L)
    res <- compareSFS(sfs1, sfs2)
    expect_true(is.list(res))
    expect_true("chisq_p" %in% names(res))
})

.sfs_with_counts <- function(counts) {
    new("WithinHostSFS",
        counts = as.integer(counts),
        breaks = seq(0.03, 0.5, length.out = length(counts) + 1L),
        folded = TRUE, polarization = "none", sampleId = "s",
        genomeLength = 1000L, nSites = as.integer(sum(counts)),
        meanDepth = 1000, threshold = 0.03, corrected = FALSE)
}

test_that("compareSFS reports small expected counts instead of warning", {
    sfs1 <- .sfs_with_counts(c(3, 1, 0, 2))
    sfs2 <- .sfs_with_counts(c(1, 2, 1, 0))
    expect_no_warning(res <- compareSFS(sfs1, sfs2))
    mat <- rbind(c(3, 1, 0, 2), c(1, 2, 1, 0))
    reference <- suppressWarnings(stats::chisq.test(mat))
    expect_equal(res$chisq_stat, unname(reference$statistic))
    expect_equal(res$chisq_p, reference$p.value)
    expected <- outer(rowSums(mat), colSums(mat)) / sum(mat)
    expect_equal(res$min_expected, min(expected))
    expect_false(res$chisq_approx_ok)
})

test_that("compareSFS flags the approximation as reliable for large counts", {
    sfs1 <- .sfs_with_counts(c(40, 30, 20, 10))
    sfs2 <- .sfs_with_counts(c(35, 35, 15, 15))
    res <- compareSFS(sfs1, sfs2)
    expect_gte(res$min_expected, 5)
    expect_true(res$chisq_approx_ok)
})

test_that("compareSFS returns NA diagnostics when fewer than two bins remain", {
    sfs1 <- .sfs_with_counts(c(4, 0, 0))
    sfs2 <- .sfs_with_counts(c(2, 0, 0))
    res <- compareSFS(sfs1, sfs2)
    expect_true(is.na(res$chisq_p))
    expect_true(is.na(res$min_expected))
    expect_true(is.na(res$chisq_approx_ok))
})

test_that("compareSFS() runs no test when a spectrum is empty (was: NaN on positive df)", {
    sfs1 <- .sfs_with_counts(c(3, 2, 2, 1))
    sfs2 <- .sfs_with_counts(c(0, 0, 0, 0))
    expect_no_warning(res <- compareSFS(sfs1, sfs2))
    expect_identical(res$chisq_stat, NA_real_)
    expect_identical(res$chisq_p, NA_real_)
    expect_identical(res$df, NA_integer_)
    expect_identical(res$min_expected, NA_real_)
    expect_identical(res$chisq_approx_ok, NA)
    expect_equal(c(res$n1, res$n2), c(8, 0))
    expect_equal(res$proportions2, rep(0, 4))
    expect_identical(compareSFS(sfs2, sfs1)$df, NA_integer_)
    ## the same through buildSFS(): every call of b lies below threshold
    whe <- .two_sample_whe(a = c(0.1, 0.2, 0.3), b = c(0.01, 0.02, NA))
    res <- compareSFS(buildSFS(whe, sampleIdx = "a", nBins = 4L),
                      buildSFS(whe, sampleIdx = "b", nBins = 4L))
    expect_identical(res$chisq_p, NA_real_)
    expect_identical(res$chisq_approx_ok, NA)
})

test_that("compareSFS rejects mismatched bins", {
    whe <- .make_whe()
    sfs1 <- buildSFS(whe, nBins = 10, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, nBins = 20, genomeLength = 1000L)
    expect_error(compareSFS(sfs1, sfs2), "same number of bins")
})

# ---- show ----

test_that("show method prints without error", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_output(show(sfs), "WithinHostSFS")
})

# ---- polarisation ----

test_that("buildSFS records the polarisation the folding implies", {
    whe <- .make_whe()
    expect_identical(buildSFS(whe, genomeLength = 1000L)@polarization,
                     "none")
    expect_identical(buildSFS(whe, fold = FALSE,
                              genomeLength = 1000L)@polarization,
                     "reference")
})

test_that("correctSFSBias keeps the polarisation", {
    whe <- .make_whe()
    for (fold in c(TRUE, FALSE)) {
        sfs <- buildSFS(whe, fold = fold, genomeLength = 1000L)
        for (m in c("truncation", "binomial"))
            expect_identical(correctSFSBias(sfs, method = m)@polarization,
                             sfs@polarization, info = paste(fold, m))
    }
})

test_that("validity ties the polarisation to the folding", {
    expect_error(new("WithinHostSFS", counts = 0L, breaks = c(0.03, 0.5),
                     folded = TRUE, polarization = "reference"),
                 "folded spectrum carries no polarisation")
    expect_error(new("WithinHostSFS", counts = 0L, breaks = c(0.03, 1),
                     folded = FALSE, polarization = "none"),
                 "must say which allele is ancestral")
    expect_error(new("WithinHostSFS", counts = 0L, breaks = c(0.03, 1),
                     folded = FALSE, polarization = "derived"),
                 "polarization must be one of")
    ## the prototype is unfolded, so its polarisation must not be "none"
    proto <- new("WithinHostSFS")
    expect_false(proto@folded)
    expect_identical(proto@polarization, "reference")
    expect_true(validObject(new("WithinHostSFS", counts = 0L,
        breaks = c(0.03, 1), folded = FALSE, polarization = "outgroup")))
})

test_that("neutralityFromSFS refuses a folded spectrum (was: a message)", {
    whe <- .make_whe(20)
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_error(neutralityFromSFS(sfs),
                 "Fay and Wu's H needs derived-allele frequencies")
    expect_error(neutralityFromSFS(sfs), "calcNeutralityTests")
})

test_that("buildSFS takes the genome length from the object (was: NA unless typed)", {
    whe <- .make_whe(20)
    expect_true(is.na(buildSFS(whe)@genomeLength))
    GenomeInfoDb::seqlengths(whe) <- c(seg = 2500L)
    expect_identical(buildSFS(whe)@genomeLength, 2500L)
    ## an explicit length still wins
    expect_identical(buildSFS(whe, genomeLength = 1000L)@genomeLength, 1000L)
})

test_that("neutralityFromSFS stops without a genome length (was: per-genome values labelled per site)", {
    whe <- .make_whe(20)
    sfs <- buildSFS(whe, fold = FALSE)
    expect_error(neutralityFromSFS(sfs), "needs the genome length")
    GenomeInfoDb::seqlengths(whe) <- c(seg = 2500L)
    expect_equal(neutralityFromSFS(buildSFS(whe, fold = FALSE)),
                 neutralityFromSFS(buildSFS(whe, fold = FALSE,
                                            genomeLength = 2500L)))
})

test_that("neutralityFromSFS reports which allele it treats as ancestral", {
    whe <- .make_whe(20)
    res <- neutralityFromSFS(buildSFS(whe, fold = FALSE,
                                      genomeLength = 1000L))
    expect_identical(res$polarization, "reference")
})

test_that("compareSFS refuses spectra with different breaks", {
    whe <- .make_whe(20)
    ## same number of bins, different threshold, so different breaks
    sfs1 <- buildSFS(whe, threshold = 0.03, nBins = 5, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, threshold = 0.10, nBins = 5, genomeLength = 1000L)
    expect_identical(length(sfs1@counts), length(sfs2@counts))
    expect_error(compareSFS(sfs1, sfs2), "must share their bin breaks")
})

test_that("compareSFS refuses a folded and an unfolded spectrum", {
    folded <- new("WithinHostSFS", counts = c(3L, 2L),
                  breaks = c(0.03, 0.5, 1), folded = TRUE,
                  polarization = "none")
    unfolded <- new("WithinHostSFS", counts = c(1L, 4L),
                    breaks = c(0.03, 0.5, 1), folded = FALSE,
                    polarization = "reference")
    expect_error(compareSFS(folded, unfolded),
                 "sfs1: folded; sfs2: reference")
})

test_that("compareSFS refuses bias-corrected spectra", {
    whe <- .make_whe(20)
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_error(compareSFS(sfs, correctSFSBias(sfs)),
                 "bias-corrected spectrum holds rescaled counts")
})

test_that("show prints the polarisation and labels the depth a median", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, fold = FALSE, genomeLength = 1000L)
    expect_output(show(sfs), "polarization: reference")
    expect_output(show(sfs), "median depth")
})
