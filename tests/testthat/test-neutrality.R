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
    res <- calcNeutralityTests(whe, genomeLength = 1000L)
    expect_s4_class(res, "DataFrame")
    expect_true("sample_id" %in% colnames(res))
    expect_true("tajimaD" %in% colnames(res))
    expect_true("fusFs" %in% colnames(res))
    expect_equal(nrow(res), 1L)
    expect_true(is.numeric(res$tajimaD))
    expect_true(is.numeric(res$fusFs))
})

test_that("calcNeutralityTests errors on missing genomeLength", {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    expect_error(calcNeutralityTests(whe),
        "genomeLength.*required")
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

test_that("fusFs returns 0 when S is 0", {
    result <- fusFs(S = 0, n = 100, theta_pi = 5)
    expect_equal(result, 0)
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
    expect_warning(res <- calcNeutralityTests(whe, genomeLength = 1000L),
                   "NA for 2 of 3 sample\\(s\\)")
    expect_identical(res$n_sites, c(5L, 2L, 1L))
    expect_true(is.finite(res$tajimaD[1]))
    expect_true(all(is.na(res$tajimaD[2:3])))
    ## Fu's Fs is still computed, since it is defined on one site
    expect_true(all(is.finite(res$fusFs)))
})

test_that("calcNeutralityTests is silent when every sample has three sites", {
    whe <- .neutrality_whe(c(4, 3))
    expect_no_warning(res <- calcNeutralityTests(whe, genomeLength = 1000L))
    expect_identical(res$n_sites, c(4L, 3L))
    expect_true(all(is.finite(res$tajimaD)))
})

test_that("n_sites counts only the QC-passed sites it used", {
    whe <- .neutrality_whe(c(5))
    qc <- matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE), ncol = 1,
                 dimnames = dimnames(SummarizedExperiment::assay(whe,
                                                                 "altFreq")))
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <- qc
    expect_identical(calcNeutralityTests(whe, genomeLength = 1000L)$n_sites,
                     3L)
    expect_identical(calcNeutralityTests(whe, genomeLength = 1000L,
                                         usePassedOnly = FALSE)$n_sites, 5L)
})
