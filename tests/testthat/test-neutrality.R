# tests/testthat/test-neutrality.R
# Neutrality test suite

test_that("tajimaD returns 0 when pi equals theta_W", {
    ## When pi_hat * L == S / a1, the numerator is zero => D = 0
    n <- 100L
    S <- 10L
    L <- 10000L
    a1 <- sum(1 / seq_len(n - 1L))
    theta_W <- S / a1
    ## Set pi_hat so that pi_hat * L == theta_W
    pi_hat <- theta_W / L
    D <- tajimaD(S = S, n = n, pi_hat = pi_hat, genomeLength = L)
    expect_equal(D, 0, tolerance = 1e-10)
})

test_that("tajimaD returns negative for excess rare variants", {
    ## Many segregating sites but low pi => excess rare variants => D < 0
    n <- 200L
    S <- 50L
    L <- 10000L
    ## pi_hat much smaller than theta_W / L
    pi_hat <- 1e-6
    D <- tajimaD(S = S, n = n, pi_hat = pi_hat, genomeLength = L)
    expect_true(is.numeric(D))
    expect_true(D < 0)
})

test_that("fusFs returns numeric", {
    result <- fusFs(S = 10, n = 50, theta_pi = 5)
    expect_true(is.numeric(result))
    expect_length(result, 1L)
})

test_that("calcNeutralityTests works on WHE from test VCF", {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    res <- calcNeutralityTests(whe)
    expect_s4_class(res, "DataFrame")
    expect_true("sample_id" %in% colnames(res))
    expect_true("tajimaD" %in% colnames(res))
    expect_true("fusFs" %in% colnames(res))
    expect_equal(nrow(res), 1L)
    expect_true(is.numeric(res$tajimaD))
    expect_true(is.numeric(res$fusFs))
})

test_that("calcNeutralityTests needs no genome length (was: required, without effect)", {
    ## five calls at a median depth of 1600, so the sample size is capped
    ## at 100
    freq <- c(0.35, 0.08, 0.12, 0.20, 0.30)
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(c(45, 112, 230, 401, 889), width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(freq, ncol = 1),
                      totalDepth = matrix(c(1500L, 800L, 1600L, 2000L,
                                            1700L), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "donor"))
    res <- calcNeutralityTests(whe)
    ## the statistics use pi summed over the genome, so a length that the
    ## old interface demanded only cancelled out
    S <- res$n_sites
    expect_identical(S, 5L)
    n <- 100L
    theta_pi <- 2 * sum(freq * (1 - freq)) * n / (n - 1L)
    expect_equal(res$tajimaD, tajimaD(S, n, theta_pi / 29903, 29903))
    expect_equal(res$fusFs, fusFs(S, n, theta_pi))
    expect_error(calcNeutralityTests(whe, genomeLength = 1000L),
                 "Unused argument")
})

test_that("tajimaD gives NA below three segregating sites (was: 0 at S = 0)", {
    ## With one site D is a monotone function of that site's frequency and
    ## with none it is 0/0; neither describes a spectrum.
    for (S in 0:2) {
        expect_warning(D <- tajimaD(S = S, n = 1000, pi_hat = 0.001,
                                    genomeLength = 10000),
                       "at least three segregating sites")
        expect_true(is.na(D), info = paste("S =", S))
    }
    expect_warning(tajimaD(S = 1, n = 1000, pi_hat = 0.001,
                           genomeLength = 10000), "with S = 1")
    ## three sites is enough
    expect_no_warning(D3 <- tajimaD(S = 3, n = 100, pi_hat = 1e-4,
                                    genomeLength = 10000))
    expect_true(is.finite(D3))
})

test_that("fusFs() returns NA with a warning at S = 0 and at theta_pi = 0 (was: 0)", {
    ## at S = 0, S' = P(K >= 1) = 1 and Fs = log(S' / (1 - S')) is infinite
    expect_warning(res <- fusFs(S = 0, n = 100, theta_pi = 5),
                   "no segregating site")
    expect_identical(res, NA_real_)
    ## with sites but no diversity there is no theta to condition on
    expect_warning(res <- fusFs(S = 3, n = 100, theta_pi = 0),
                   "theta_pi = 0")
    expect_identical(res, NA_real_)
    expect_warning(res <- fusFs(S = 0, n = 100, theta_pi = 0),
                   "no segregating site")
    expect_identical(res, NA_real_)
})

test_that("fusFs() is exact above n = 500, so it is continuous there (was: a normal approximation that jumped)", {
    expect_lt(abs(fusFs(40, 500, 3) - fusFs(40, 501, 3)), 0.1)
    expect_lt(abs(fusFs(20, 500, 3) - fusFs(20, 501, 3)), 0.1)
    expect_lt(abs(fusFs(5, 500, 3) - fusFs(5, 501, 3)), 0.1)
    ## an independent exact tail: under the Ewens sampling formula K is a
    ## sum of independent Bernoulli(theta / (theta + i - 1)), i = 1..n
    ewens_fs <- function(S, n, theta) {
        p <- theta / (theta + seq_len(n) - 1)
        dist <- 1
        for (pr in p) dist <- c(dist * (1 - pr), 0) + c(0, dist * pr)
        s_prime <- sum(dist[(S + 2):(n + 1)])
        log(s_prime / (1 - s_prime))
    }
    expect_equal(fusFs(20, 800, 3), ewens_fs(20, 800, 3), tolerance = 1e-6)
    expect_equal(fusFs(10, 2000, 5), ewens_fs(10, 2000, 5),
                 tolerance = 1e-6)
    ## the branch at or below 500 agrees with it too
    expect_equal(fusFs(10, 100, 5), ewens_fs(10, 100, 5), tolerance = 1e-6)
})

test_that("fusFs() lies well below 0 at the neutral expectation, the values its help page gives", {
    ## the haplotype count is S + 1, the infinite-sites maximum, so the
    ## value sits far from 0 even when S and theta_pi agree with neutrality
    a1 <- sum(1 / seq_len(99))
    expect_identical(round(c(2, 5, 10) * a1), c(10, 26, 52))
    expect_equal(fusFs(10, 100, 2), -1.463007, tolerance = 1e-6)
    expect_equal(fusFs(26, 100, 5), -6.998362, tolerance = 1e-6)
    expect_equal(fusFs(52, 100, 10), -25.01018, tolerance = 1e-6)
})

test_that("tajimaD validates inputs", {
    expect_error(tajimaD(S = -1, n = 10, pi_hat = 0.001,
        genomeLength = 1000), "non-negative")
    expect_error(tajimaD(S = 5, n = 1, pi_hat = 0.001,
        genomeLength = 1000), "must be a single number > 1")
    expect_error(tajimaD(S = 5, n = 100, pi_hat = 0.001,
        genomeLength = 0), "positive")
})

test_that("tajimaD returns NA with a warning when n rounds below four (was: Inf or NaN)", {
    ## with two or three sequences c1 = c2 = 0, so the variance of D is zero
    for (n in c(2, 3, 3.4)) {
        expect_warning(d <- tajimaD(S = 5, n = n, pi_hat = 0.001,
                                    genomeLength = 1000), "at least four")
        expect_true(is.na(d))
    }
    expect_true(is.finite(tajimaD(S = 5, n = 4, pi_hat = 0.001,
                                  genomeLength = 1000)))
})

test_that("fusFs returns NA with a warning when S + 1 exceeds n", {
    expect_warning(res <- fusFs(S = 100, n = 100, theta_pi = 5),
                   "undefined")
    expect_true(is.na(res))
    expect_true(is.finite(fusFs(S = 10, n = 100, theta_pi = 5)))
})

# ---- n_sites and the three-site rule in calcNeutralityTests ----

.neutrality_whe <- function(counts) {
    ## one column per entry of counts, with that many variant sites present
    n <- max(counts)
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(seq_len(n) * 10L, width = 1))
    freq <- vapply(counts, function(k)
        c(seq(0.05, 0.45, length.out = k), rep(NA_real_, n - k))[seq_len(n)],
        numeric(n))
    colnames(freq) <- paste0("s", seq_along(counts))
    WithinHostExperiment(
        assays = list(altFreq = freq,
                      totalDepth = matrix(1000L, nrow = n,
                                          ncol = length(counts),
                                          dimnames = dimnames(freq))),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = colnames(freq),
                                       row.names = colnames(freq)))
}

test_that("calcNeutralityTests reports n_sites and NA D below three sites", {
    whe <- .neutrality_whe(c(5, 2, 1))
    expect_warning(res <- calcNeutralityTests(whe),
                   "NA for 2 of 3 sample\\(s\\)")
    expect_identical(res$n_sites, c(5L, 2L, 1L))
    expect_true(is.finite(res$tajimaD[1]))
    expect_true(all(is.na(res$tajimaD[2:3])))
    ## Fu's Fs is still computed, since it is defined on one site
    expect_true(all(is.finite(res$fusFs)))
})

test_that("calcNeutralityTests is silent when every sample has three sites", {
    whe <- .neutrality_whe(c(4, 3))
    expect_no_warning(res <- calcNeutralityTests(whe))
    expect_identical(res$n_sites, c(4L, 3L))
    expect_true(all(is.finite(res$tajimaD)))
})

.fixed_site_whe <- function(freq) {
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(seq_along(freq) * 10L, width = 1))
    WithinHostExperiment(
        assays = list(altFreq = matrix(freq, ncol = 1),
                      totalDepth = matrix(1000L, nrow = length(freq),
                                          ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "s"))
}

test_that("calcNeutralityTests() counts only 0 < f < 1 as segregating (was: absent and fixed sites in S)", {
    ## no QC has been run, so qcPass is TRUE at every site
    whe <- .fixed_site_whe(c(0.05, 0.2, 0.3, 0, 1, 1, 1))
    res <- calcNeutralityTests(whe)
    three <- calcNeutralityTests(whe[1:3, ])
    expect_identical(res$n_sites, 3L)
    expect_equal(res$tajimaD, three$tajimaD)
    expect_equal(res$tajimaD, 0.8175039, tolerance = 1e-6)
    expect_equal(res$fusFs, three$fusFs)
    ## a sample with no polymorphism at all has no segregating site
    ## (was: n_sites 4 and D = -1.98)
    expect_warning(expect_warning(
        none <- calcNeutralityTests(.fixed_site_whe(c(0, 1, 1, 1))),
        "fewer than three"), "no segregating site")
    expect_identical(none$n_sites, 0L)
    expect_true(is.na(none$tajimaD))
    expect_true(is.na(none$fusFs))
})

test_that("calcNeutralityTests() gives one warning per cause of NA (was: per-sample warnings, and Fs = 0 at S = 0)", {
    ## s1: five sites; s2: none; s3: four sites at a depth of 3, so n = 3;
    ## s4: 100 sites at depth 1000, so S + 1 exceeds n = 100
    whe <- .neutrality_whe(c(5, 0, 4, 100))
    depth <- SummarizedExperiment::assay(whe, "totalDepth")
    depth[, 3] <- 3L
    SummarizedExperiment::assay(whe, "totalDepth") <- depth
    msgs <- character()
    res <- withCallingHandlers(calcNeutralityTests(whe),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_length(msgs, 4L)
    expect_match(msgs[1], "Tajima's D is NA for 1 of 4 .*fewer than three")
    expect_match(msgs[2], "Tajima's D is NA for 1 of 4 .*below four")
    expect_match(msgs[3], "Fu's Fs is NA for 1 of 4 .*no segregating site")
    expect_match(msgs[4],
                 "Fu's Fs is NA for 2 of 4 .*not below the effective sample")
    expect_identical(is.na(res$tajimaD), c(FALSE, TRUE, TRUE, FALSE))
    expect_identical(is.na(res$fusFs), c(FALSE, TRUE, TRUE, TRUE))
    expect_identical(res$n_sites, c(5L, 0L, 4L, 100L))
})

test_that("n_sites counts only the QC-passed sites it used", {
    whe <- .neutrality_whe(c(5))
    qc <- matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE), ncol = 1,
                 dimnames = dimnames(SummarizedExperiment::assay(whe,
                                                                 "altFreq")))
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <- qc
    expect_identical(calcNeutralityTests(whe)$n_sites,
                     3L)
    expect_identical(calcNeutralityTests(whe,
                                         usePassedOnly = FALSE)$n_sites, 5L)
})
