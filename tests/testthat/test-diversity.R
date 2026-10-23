test_that("shannonISNV computes correct value", {
    expect_equal(shannonISNV(0.5), log(2), tolerance = 1e-10)
    expect_equal(shannonISNV(numeric(0)), 0)
    expect_equal(shannonISNV(c(NA, NA)), 0)
    expect_equal(shannonISNV(c(0, 0.5, 1)), log(2), tolerance = 1e-10)
})

test_that("piISNV computes correct value", {
    expect_equal(piISNV(0.5, 1000), 0.0005, tolerance = 1e-10)
    expect_equal(piISNV(numeric(0), 1000), 0)
})

test_that("piISNV with finite-sample correction", {
    raw <- piISNV(0.5, 1000)
    corrected <- piISNV(0.5, 1000, meanDepth = 100)
    ## Correction factor = 100/99
    expect_equal(corrected, raw * 100 / 99, tolerance = 1e-10)
    ## Higher depth -> correction closer to 1
    deep <- piISNV(0.5, 1000, meanDepth = 10000)
    expect_true(abs(deep - raw) < abs(corrected - raw))
})

test_that("piISNV rejects a mean depth that rounds to one", {
    ## A depth above 1 that rounds to 1 used to give n/(n-1) = Inf
    expect_error(piISNV(0.5, 1000, meanDepth = 1.2), "rounds to 1")
    expect_error(piISNV(0.5, 1000, meanDepth = 1.4), "rounds to 1")
    expect_true(is.finite(piISNV(0.5, 1000, meanDepth = 1.5)))
    expect_error(piISNV(0.5, 1000, meanDepth = 1), "greater than 1|> 1")
})

test_that("piISNV(), shannonISNV() and simpsonISNV() refuse frequencies outside [0, 1] (was: a negative pi, values dropped)", {
    expect_error(piISNV(c(0.1, 2), 1000), "'freq'.*between 0 and 1")
    expect_error(shannonISNV(c(0.1, 1.5)), "'freq'.*between 0 and 1")
    expect_error(simpsonISNV(c(-0.1, 0.2)), "'freq'.*between 0 and 1")
    ## text was coerced to NA and gave 0 with only a coercion warning
    expect_error(shannonISNV("a"), "'freq' must be a numeric")
    expect_error(simpsonISNV("0.2"), "'freq' must be a numeric")
    expect_error(piISNV("0.2", 1000), "'freq' must be a numeric")
    ## 0, 1 and NA are still accepted, and 0 and 1 still add nothing
    expect_equal(shannonISNV(c(0, 0.5, 1, NA)), log(2))
    expect_equal(simpsonISNV(c(0, 0.5, 1, NA)), 0.5)
    expect_equal(piISNV(c(0, 0.5, 1, NA), 1000), 0.0005)
    expect_equal(shannonISNV(NA), 0)
})

test_that("piISNV rejects invalid genomeLength", {
    expect_error(piISNV(0.5, -1))
    expect_error(piISNV(0.5, 0))
})

test_that("wattersonISNV returns positive value", {
    w <- wattersonISNV(15, 13588, 2000)
    expect_true(is.numeric(w))
    expect_true(w > 0)
})

test_that("wattersonISNV returns 0 for 0 sites", {
    expect_equal(wattersonISNV(0, 1000, 500), 0)
})

## the reference the synthetic calls were made against, 1000 bp
test_ref_dna <- function() {
    Biostrings::readDNAStringSet(system.file("extdata", "test_ref.fa",
        package = "WithinHostExperiment"))
}

test_that("calcDiversity works on WHE", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar", genome = test_ref_dna())
    div <- calcDiversity(whe)
    expect_s4_class(div, "DataFrame")
    expect_true("shannon" %in% colnames(div))
    expect_true("pi" %in% colnames(div))
    expect_true("richness" %in% colnames(div))
    expect_true(div$shannon > 0)
    expect_true(div$richness > 0)
})

test_that("calcDiversity stops when the genome length is not recorded", {
    ## (was: a genomeLength argument to type)
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    expect_true(is.na(GenomeInfoDb::seqlengths(whe)[["test_segment"]]))
    expect_error(calcDiversity(whe, indices = "pi"),
                 "not recorded for test_segment")
    expect_error(calcDiversity(whe, indices = c("richness", "watterson")),
                 'seqlengths(x)[c("test_segment")] <- <length>', fixed = TRUE)
    ## indices that are not per site do not need it
    expect_no_error(calcDiversity(whe, indices = c("richness", "shannon")))
    expect_error(calcDiversity(whe, genomeLength = 1000L), "Unused argument")
    ## recording it on the object is enough
    GenomeInfoDb::seqlengths(whe) <- c(test_segment = 1000L)
    expect_no_error(div <- calcDiversity(whe, indices = "pi"))
    expect_true(div$pi > 0)
})

test_that("calcDiversity warns when totalDepth missing for watterson", {
    gr <- GenomicRanges::GRanges("s",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = S4Vectors::DataFrame(sample_id = "S1"))
    GenomeInfoDb::seqlengths(whe) <- c(s = 1000L)
    expect_warning(
        calcDiversity(whe, indices = "watterson"),
        "totalDepth"
    )
})

test_that("calcDiversity respects QC pass filter", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    div_before <- calcDiversity(whe, indices = "richness")

    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    div_after <- calcDiversity(whe, indices = "richness")
    expect_true(div_after$richness <= div_before$richness)
})

# --- simpsonISNV tests ---

test_that("simpsonISNV computes correct value", {
    ## Single site at p=0.5: 1 - (0.25 + 0.25) = 0.5
    expect_equal(simpsonISNV(0.5), 0.5, tolerance = 1e-10)
    ## Empty input
    expect_equal(simpsonISNV(numeric(0)), 0)
    expect_equal(simpsonISNV(c(NA, NA)), 0)
    ## 0 and 1 are excluded, only 0.5 counted
    expect_equal(simpsonISNV(c(0, 0.5, 1)), 0.5, tolerance = 1e-10)
})

test_that("simpsonISNV is higher for more even single site", {
    ## p=0.5 is maximally diverse, p=0.1 is less diverse
    expect_true(simpsonISNV(0.5) > simpsonISNV(0.1))
})

# --- chao1ISNV tests ---

test_that("chao1ISNV computes correct value", {
    ## 2 singletons, 1 doubleton: 5 + 4/2 = 7
    expect_equal(chao1ISNV(c(1, 1, 2, 5, 10), detected = 5), 7)
})

test_that("chao1ISNV handles zero doubletons", {
    ## f1=2, f2=0: bias-corrected = 5 + 2*1/2 = 6
    expect_equal(chao1ISNV(c(1, 1, 3, 5, 10), detected = 5), 6)
})

test_that("chao1ISNV returns detected when no singletons", {
    expect_equal(chao1ISNV(c(3, 4, 5), detected = 3), 3)
})

test_that("chao1ISNV handles empty counts", {
    expect_equal(chao1ISNV(integer(0), detected = 5), 5)
})

# --- calcDiversity with simpson and chao1 ---

test_that("calcDiversity supports simpson index", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    div <- calcDiversity(whe, indices = "simpson")
    expect_s4_class(div, "DataFrame")
    expect_true("simpson" %in% colnames(div))
    expect_true(is.numeric(div$simpson))
})

test_that("calcDiversity supports chao1 index", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    div <- calcDiversity(whe, indices = "chao1")
    expect_s4_class(div, "DataFrame")
    expect_true("chao1" %in% colnames(div))
    expect_true(div$chao1 >= 0)
})

test_that("calcDiversity warns when altCount missing for chao1", {
    gr <- GenomicRanges::GRanges("s",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = S4Vectors::DataFrame(sample_id = "S1"))
    expect_warning(
        calcDiversity(whe, indices = "chao1"),
        "altCount"
    )
})

test_that("simpsonISNV sums per-site diversity and is never negative", {
    freq <- c(0.05, 0.1, 0.2, 0.3)
    expect_equal(simpsonISNV(freq), sum(2 * freq * (1 - freq)),
                 tolerance = 1e-12)
    expect_true(simpsonISNV(c(0.1, 0.2)) > 0)
    expect_equal(simpsonISNV(c(0.5, 0.5)), 1, tolerance = 1e-12)
})

# --- provenance columns and the genome length ---

test_that("n_sites_assessed and detection_limit tell a zero apart", {
    gr <- GenomicRanges::GRanges("s", IRanges::IRanges(c(10, 20, 30),
                                                      width = 1))
    GenomeInfoDb::seqlengths(gr) <- c(s = 1000L)
    freq <- cbind(a = c(0.2, 0.3, NA), b = c(NA, NA, NA), c = c(0.01, 0.02, NA))
    depth <- matrix(1000L, nrow = 3, ncol = 3, dimnames = dimnames(freq))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
                                       row.names = colnames(freq)))
    ## before any QC: sites assessed are the non-missing calls
    div <- calcDiversity(whe, indices = c("richness", "pi"))
    expect_identical(div$n_sites_assessed, c(2L, 0L, 2L))
    expect_true(all(is.na(div$detection_limit)))
    expect_identical(div$richness, c(2L, 0L, 2L))

    whe <- flagISNV(whe, ISNVFilter(minDepth = 100L, minFreq = 0.03))
    div <- calcDiversity(whe, indices = c("richness", "pi"))
    ## sample c had two calls, both below the limit: assessed 2, kept 0,
    ## pi 0 for that reason; sample b had no call at all
    expect_identical(div$n_sites_assessed, c(2L, 0L, 2L))
    expect_identical(div$richness, c(2L, 0L, 0L))
    expect_equal(div$detection_limit, rep(0.03, 3))
    expect_equal(div$pi[2:3], c(0, 0))
    ## the strictest filter applied sets the limit
    whe2 <- flagISNV(whe, ISNVFilter(minDepth = 100L, minFreq = 0.10))
    expect_equal(calcDiversity(whe2, indices = "richness")$detection_limit,
                 rep(0.10, 3))
    ## without QC the limit is unknown
    no_qc <- calcDiversity(whe2, indices = "richness", usePassedOnly = FALSE)
    expect_true(all(is.na(no_qc$detection_limit)))
})

test_that("the genome length sums every recorded sequence (was: only those with sites)", {
    gr <- GenomicRanges::GRanges(c("seg1", "seg2"),
                                 IRanges::IRanges(c(10, 20), width = 1))
    GenomeInfoDb::seqlengths(gr) <- c(seg1 = 600L, seg2 = 400L)
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    div <- calcDiversity(whe, indices = "pi")
    expect_equal(div$pi, piISNV(c(0.2, 0.3), genomeLength = 1000))
    ## a segment without calls is sequence in which none was found, so
    ## it stays in the denominator
    only1 <- whe[1, ]
    expect_identical(GenomeInfoDb::seqlevelsInUse(only1), "seg1")
    expect_equal(calcDiversity(only1, indices = "pi")$pi,
                 piISNV(0.2, genomeLength = 1000))
    ## one unknown length is enough to stop, and the advice names it and
    ## works as written (was: it named the first sequence)
    GenomeInfoDb::seqlengths(whe) <- c(seg1 = 600L, seg2 = NA_integer_)
    expect_error(calcDiversity(whe, indices = "pi"), "not recorded for seg2")
    expect_error(calcDiversity(whe, indices = "pi"),
                 'seqlengths(x)[c("seg2")] <- <length>', fixed = TRUE)
    ## the second recipe applies to a sequence without calls: keep only
    ## the site on seg1, and drop seg2, whose length is unknown
    seg1_only <- whe[1, ]
    expect_error(calcDiversity(seg1_only, indices = "pi"),
                 "if it carries no call")
    pruned <- GenomeInfoDb::keepSeqlevels(seg1_only,
        GenomeInfoDb::seqlevelsInUse(seg1_only), pruning.mode = "coarse")
    expect_equal(calcDiversity(pruned, indices = "pi")$pi,
                 piISNV(0.2, genomeLength = 600))
    GenomeInfoDb::seqlengths(whe)[c("seg2")] <- 400L
    expect_equal(calcDiversity(whe, indices = "pi")$pi,
                 piISNV(c(0.2, 0.3), genomeLength = 1000))
})

test_that("compareDiversity takes the genome length from the object", {
    data(example_whe)
    ## donor and recipient: two groups of one sample each
    whe <- example_whe[, 1:2]
    GenomeInfoDb::seqlengths(whe) <- c(test_segment = 1000L)
    res <- compareDiversity(whe, groupBy = "role", index = "pi")
    expect_s4_class(res, "DataFrame")
    expect_true("note" %in% colnames(res))
    GenomeInfoDb::seqlengths(whe) <- c(test_segment = NA_integer_)
    expect_error(compareDiversity(whe, groupBy = "role", index = "pi"),
                 "not recorded")
    expect_error(compareDiversity(whe, groupBy = "role", index = "pi",
                                  genomeLength = 1000L), "unused argument")
})

test_that("a sample with no calls gives zeros, not an error (was: missing value where TRUE/FALSE needed)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20), width = 1))
    GenomeInfoDb::seqlengths(gr) <- c(seg = 1000L)
    freq <- cbind(a = c(0.2, 0.3), b = c(NA, NA))
    depth <- cbind(a = c(1000L, 1000L), b = c(NA_integer_, NA_integer_))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("a", "b"),
                                       row.names = c("a", "b")))
    div <- calcDiversity(whe, indices = c("richness", "pi", "watterson"))
    expect_identical(div$n_sites_assessed, c(2L, 0L))
    expect_identical(div$richness, c(2L, 0L))
    expect_equal(div$pi[2], 0)
    expect_equal(div$watterson[2], 0)
    expect_warning(nt <- calcNeutralityTests(whe), "fewer than three")
    expect_identical(nt$n_sites, c(2L, 0L))
    expect_true(all(is.na(nt$tajimaD)))
    expect_equal(nt$fusFs[2], 0)
})
