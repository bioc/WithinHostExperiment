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
                        "n_variants", "pair_id"))
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
    res <- exactBottleneck(whe, pairId = "p1", maxNb = 30L,
                           method = "wright_fisher")
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
                        "n_variants", "pair_id"))
    expect_equal(res$method, "presence_absence")
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
