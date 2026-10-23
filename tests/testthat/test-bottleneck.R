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

test_that("quickBottleneck returns valid result", {
    whe <- .make_pair_whe()
    nb  <- quickBottleneck(whe, pairId = "p1", maxNb = 50L)
    expect_true(is.list(nb))
    expect_true(is.finite(nb$Nb))
    expect_true(nb$Nb >= 1)
    expect_equal(length(nb$ci), 2L)
    expect_true(nb$ci[1] <= nb$Nb)
    expect_true(nb$ci[2] >= nb$Nb)
    expect_true(nb$n_variants > 0L)
})

test_that(".bbSiteLogLik returns finite value", {
    ll <- WithinHostExperiment:::.bbSiteLogLik(0.3, 0.25, 5L, 0.03)
    expect_true(is.finite(ll))
    # Log-density can be positive for beta distribution
    expect_true(is.numeric(ll))
})

test_that(".logSumExp is numerically stable", {
    x      <- c(-1000, -1001, -999)
    result <- WithinHostExperiment:::.logSumExp(x)
    expect_true(is.finite(result))
    expect_true(result > -1000)
})

test_that("quickBottleneck() and exactBottleneck() check threshold (was: base-R errors, or no window at all)", {
    whe <- .make_pair_whe()
    for (bad in list(0, -0.1, 0.5, 0.7, NA, c(0.03, 0.05), "0.03")) {
        expect_error(quickBottleneck(whe, pairId = "p1", threshold = bad),
                     "'threshold' must be a single number above 0 and below 0.5")
        expect_error(exactBottleneck(whe, pairId = "p1", threshold = bad,
                                     method = "presence_absence"),
                     "'threshold' must be a single number above 0 and below 0.5")
    }
})

test_that("quickBottleneck() refuses a maxNb below 1 or outside the integer range (was: Nb of length 0, truncation, text and logical accepted)", {
    whe <- .make_pair_whe()
    for (bad in list(0L, -3, 2.5, NA, c(10L, 20L), "50", TRUE, 3e9, Inf)) {
        expect_error(quickBottleneck(whe, pairId = "p1", maxNb = bad),
                     "'maxNb' must be a positive integer")
        expect_error(exactBottleneck(whe, pairId = "p1", maxNb = bad,
                                     method = "presence_absence"),
                     "'maxNb' must be a positive integer")
    }
    ## a whole number given as a double is accepted by both
    expect_identical(quickBottleneck(whe, pairId = "p1", maxNb = 50)$Nb,
                     quickBottleneck(whe, pairId = "p1", maxNb = 50L)$Nb)
})

test_that("quickBottleneck handles no shared variants gracefully", {
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
    nb <- quickBottleneck(whe, pairId = "p2", maxNb = 20L)
    expect_true(is.list(nb))
})
