# tests/testthat/test-bottleneck-exact.R

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

# ------------------------------------------------------------------
# 1. presence_absence returns valid structure
# ------------------------------------------------------------------
test_that("exactBottleneck presence_absence returns valid result", {
    whe <- .make_pair_whe()
    res <- exactBottleneck(whe, pairId = "p1", maxNb = 50L,
                           method = "presence_absence")
    expect_true(is.list(res))
    expect_named(res, c("Nb", "ci", "loglik", "method",
                        "n_variants", "n_excluded", "pair_id"))
    expect_equal(res$method, "presence_absence")
    expect_true(is.finite(res$Nb))
    expect_true(res$Nb >= 1)
    expect_equal(length(res$ci), 2L)
    expect_true(res$ci[1] <= res$Nb)
    expect_true(res$ci[2] >= res$Nb)
    expect_true(res$n_variants > 0L)
    expect_equal(res$pair_id, "p1")
})

# ------------------------------------------------------------------
# 2. wright_fisher returns valid structure
# ------------------------------------------------------------------
test_that("exactBottleneck wright_fisher returns valid result", {
    whe <- .make_pair_whe()
    expect_message(res <- exactBottleneck(whe, pairId = "p1", maxNb = 30L,
                                          method = "wright_fisher"),
                   "heuristic")
    expect_true(is.list(res))
    expect_equal(res$method, "wright_fisher")
    expect_true(is.finite(res$Nb))
    expect_true(res$Nb >= 1)
    expect_equal(length(res$ci), 2L)
    expect_true(res$n_variants > 0L)
})

# ------------------------------------------------------------------
# 3. Input validation errors
# ------------------------------------------------------------------
test_that("exactBottleneck rejects invalid inputs", {
    whe <- .make_pair_whe()

    expect_error(exactBottleneck("not_whe", pairId = "p1"),
                 "must be a WithinHostExperiment")
    expect_error(exactBottleneck(whe, pairId = "p1", maxNb = 0L),
                 "must be a positive integer")
    expect_error(exactBottleneck(whe, pairId = "p1", threshold = -0.1),
                 "must be between 0 and 1")
    expect_error(exactBottleneck(whe, pairId = "p1", method = "bogus"),
                 "arg")
})

# ------------------------------------------------------------------
# 4. Internal helpers: .paSiteLogLik
# ------------------------------------------------------------------
test_that(".paSiteLogLik returns correct values", {
    ## Present variant with high donor freq -> high log-lik
    ll_present <- WithinHostExperiment:::.paSiteLogLik(0.3, TRUE, 5L)
    expect_true(is.finite(ll_present))
    expect_true(ll_present < 0)

    ## Absent variant -> Nb * log(1 - nu_d)
    ll_absent <- WithinHostExperiment:::.paSiteLogLik(0.3, FALSE, 5L)
    expect_equal(ll_absent, 5 * log(1 - 0.3), tolerance = 1e-10)

    ## Higher Nb increases detection probability
    ll_low  <- WithinHostExperiment:::.paSiteLogLik(0.1, TRUE, 2L)
    ll_high <- WithinHostExperiment:::.paSiteLogLik(0.1, TRUE, 20L)
    expect_true(ll_high > ll_low)
})

# ------------------------------------------------------------------
# 5. Internal helper: .logBetaBinomPMF
# ------------------------------------------------------------------
test_that(".logBetaBinomPMF returns finite values", {
    lp <- WithinHostExperiment:::.logBetaBinomPMF(5L, 20L, 2, 3)
    expect_true(is.finite(lp))
    expect_true(lp <= 0)

    ## k=0 should be valid
    lp0 <- WithinHostExperiment:::.logBetaBinomPMF(0L, 20L, 1, 10)
    expect_true(is.finite(lp0))

    ## k=n should be valid
    lpn <- WithinHostExperiment:::.logBetaBinomPMF(20L, 20L, 10, 1)
    expect_true(is.finite(lpn))
})

# ------------------------------------------------------------------
# 6. Internal helper: .exactBBSiteLogLik
# ------------------------------------------------------------------
test_that(".exactBBSiteLogLik returns finite value", {
    ll <- WithinHostExperiment:::.exactBBSiteLogLik(
        nu_d = 0.3, k_r = 10L, n_r = 100L, Nb = 5L, threshold = 0.03)
    expect_true(is.finite(ll))
    expect_true(is.numeric(ll))
})

# ------------------------------------------------------------------
# 7. .wfSiteLogLik returns finite for present/absent
# ------------------------------------------------------------------
test_that(".wfSiteLogLik returns finite values", {
    ll_p <- WithinHostExperiment:::.wfSiteLogLik(
        0.3, 0.25, TRUE, 10L, 0.03)
    expect_true(is.finite(ll_p))

    ll_a <- WithinHostExperiment:::.wfSiteLogLik(
        0.3, 0.0, FALSE, 10L, 0.03)
    expect_true(is.finite(ll_a))
})

# ------------------------------------------------------------------
# 8. Empty result when no variants above threshold
# ------------------------------------------------------------------
test_that("exactBottleneck handles no shared variants gracefully", {
    vcf  <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_independent.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("d", "ind"),
            role      = c("donor", "recipient"),
            pair_id   = c("p2", "p2")),
        caller = "ivar")
    res <- exactBottleneck(whe, pairId = "p2", maxNb = 20L,
                           method = "presence_absence")
    expect_true(is.list(res))
    expect_named(res, c("Nb", "ci", "loglik", "method",
                        "n_variants", "n_excluded", "pair_id"))
    expect_equal(res$method, "presence_absence")
})

test_that("exact_bb leaves out donor variants of unknown recipient depth (was: counted in n_variants)", {
    ## three donor variants: two with recipient read counts, and one at a
    ## site where the recipient's depth is not recorded
    pair <- function(n) {
        gr <- GenomicRanges::GRanges("seg",
                                     IRanges::IRanges(1:3 * 10, width = 1))
        S4Vectors::mcols(gr)$ref <- "A"
        S4Vectors::mcols(gr)$alt <- "G"
        whe <- WithinHostExperiment(
            assays = list(
                altFreq = cbind(d = c(0.3, 0.2, 0.25), r = c(0.8, 0.05, NA)),
                altCount = cbind(d = c(300L, 200L, 250L),
                                 r = c(800L, 50L, NA)),
                totalDepth = cbind(d = rep(1000L, 3),
                                   r = c(1000L, 1000L, NA))),
            rowRanges = gr,
            colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
                role = c("donor", "recipient"), pair_id = c("p", "p"),
                row.names = c("d", "r")))
        whe[seq_len(n), ]
    }
    expect_message(all3 <- exactBottleneck(pair(3L), pairId = "p",
                                           maxNb = 30L, method = "exact_bb"),
                   "1 donor variant site\\(s\\) whose recipient depth")
    expect_no_message(known <- exactBottleneck(pair(2L), pairId = "p",
                                               maxNb = 30L,
                                               method = "exact_bb"))
    expect_identical(all3$n_variants, 2L)
    expect_identical(all3$n_excluded[["unknown_depth"]], 1L)
    expect_identical(known$n_excluded[["unknown_depth"]], 0L)
    expect_equal(all3[c("Nb", "ci", "loglik")], known[c("Nb", "ci", "loglik")])
    ## the read-from-VCF pair: ten of its 15 donor variants are sites the
    ## recipient's file does not list
    vcf <- suppressMessages(exactBottleneck(.make_pair_whe(), pairId = "p1",
                                            maxNb = 200L, method = "exact_bb"))
    expect_identical(vcf$n_variants, 5L)
    expect_identical(vcf$n_excluded[["unknown_depth"]], 10L)
    ## the other estimators need no depth and leave nothing out for it
    pa <- exactBottleneck(.make_pair_whe(), pairId = "p1", maxNb = 50L,
                          method = "presence_absence")
    expect_identical(pa$n_variants, 15L)
    expect_identical(pa$n_excluded[["unknown_depth"]], 0L)
})

test_that("a recipient call removed by QC is left out of every estimator (was: scored as lost)", {
    ## four donor variants; the recipient carries the first at 70%, the
    ## second at 97%, which maxFreq = 0.95 removes, lacks the third at
    ## depth 1000 and carries the fourth at 5%
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:4 * 10, width = 1))
    S4Vectors::mcols(gr)$ref <- "A"
    S4Vectors::mcols(gr)$alt <- "G"
    r_freq <- c(0.7, 0.97, NA, 0.05)
    whe <- WithinHostExperiment(
        assays = list(
            altFreq = cbind(d = c(0.3, 0.4, 0.25, 0.2), r = r_freq),
            altCount = cbind(d = c(300L, 400L, 250L, 200L),
                             r = c(700L, 970L, 0L, 50L)),
            totalDepth = cbind(d = rep(1000L, 4), r = rep(1000L, 4))),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
            role = c("donor", "recipient"), pair_id = c("p", "p"),
            row.names = c("d", "r")))
    whe <- flagISNV(whe, ISNVFilter(maxFreq = 0.95))
    expect_identical(calcSharedVariants(whe, pairId = "p")$
                         n_donor_only_recipient_failed_qc, 1L)
    left_out <- "1 donor variant site\\(s\\) where QC removed the recipient"
    expect_message(vb <- asViralBottleneckInput(whe, pairId = "p"), left_out)
    expect_identical(vb$pos, c(10L, 30L, 40L))
    ## without the QC mask the call is used as it was made
    expect_no_message(raw <- asViralBottleneckInput(whe, pairId = "p",
                                                    usePassedOnly = FALSE))
    expect_identical(raw$recipient_freq, c(0.7, 0.97, 0, 0.05))
    ## the exported table keeps its columns
    expect_identical(colnames(exportPairFrequencies(whe, pairId = "p")),
                     c("chrom", "position", "ref", "alt", "donor_freq",
                       "recipient_freq", "shared", "pair_id"))
    ## each estimator gives what it gives on the pair without that site
    without <- whe[-2L, ]
    expect_message(q <- quickBottleneck(whe, pairId = "p", maxNb = 50L),
                   left_out)
    expect_equal(q[c("Nb", "ci", "loglik", "n_variants")],
                 quickBottleneck(without, pairId = "p",
                                 maxNb = 50L)[c("Nb", "ci", "loglik",
                                                "n_variants")])
    expect_identical(q$n_excluded[["recipient_failed_qc"]], 1L)
    for (m in c("presence_absence", "exact_bb")) {
        expect_message(e <- exactBottleneck(whe, pairId = "p", maxNb = 50L,
                                            method = m), left_out)
        expect_equal(e[c("Nb", "ci", "loglik", "n_variants")],
                     exactBottleneck(without, pairId = "p", maxNb = 50L,
                                     method = m)[c("Nb", "ci", "loglik",
                                                   "n_variants")])
        expect_identical(e$n_variants, 3L)
        expect_identical(e$n_excluded[["recipient_failed_qc"]], 1L)
    }
})

test_that("presence_absence leaves out a donor variant above 1 - threshold (was: Nb = 1 at log-likelihood -Inf)", {
    ## three donor variants: one the recipient carries, one it lacks,
    ## and one fixed in the donor that the recipient lacks
    pair <- function(n) {
        gr <- GenomicRanges::GRanges("seg",
                                     IRanges::IRanges(1:3 * 10, width = 1))
        S4Vectors::mcols(gr)$ref <- "A"
        S4Vectors::mcols(gr)$alt <- "G"
        whe <- WithinHostExperiment(
            assays = list(altFreq = cbind(d = c(0.3, 0.2, 1),
                                          r = c(0.4, NA, NA))),
            rowRanges = gr,
            colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
                role = c("donor", "recipient"), pair_id = c("p", "p"),
                row.names = c("d", "r")))
        whe[seq_len(n), ]
    }
    expect_message(pa <- exactBottleneck(pair(3L), pairId = "p", maxNb = 50L,
                                         method = "presence_absence"),
                   "1 donor variant site\\(s\\) above 1 - threshold")
    two <- exactBottleneck(pair(2L), pairId = "p", maxNb = 50L,
                           method = "presence_absence")
    expect_true(is.finite(pa$loglik))
    expect_equal(pa[c("Nb", "ci", "loglik", "n_variants")],
                 two[c("Nb", "ci", "loglik", "n_variants")])
    expect_identical(pa$n_excluded[["donor_near_fixation"]], 1L)
    expect_identical(two$n_excluded[["donor_near_fixation"]], 0L)
    ## wright_fisher keeps the fixed site, whose score is -Inf at every
    ## Nb, and stops rather than return Nb = 1
    expect_error(suppressMessages(
        exactBottleneck(pair(3L), pairId = "p", maxNb = 50L,
                        method = "wright_fisher")),
        "not finite at any Nb from 1 to 50")
})

test_that("exactBottleneck says that wright_fisher is a heuristic", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role = c("donor", "recipient"),
            pair_id = c("p1", "p1")),
        caller = "ivar")
    expect_message(exactBottleneck(whe, pairId = "p1", maxNb = 20L,
                                   method = "wright_fisher"),
                   "heuristic score, not a likelihood")
    expect_no_message(exactBottleneck(whe, pairId = "p1", maxNb = 20L,
                                      method = "presence_absence"))
})
