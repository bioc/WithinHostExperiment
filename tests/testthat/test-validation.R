test_that("checkDiversityPlausibility flags out-of-range pi", {
    div <- S4Vectors::DataFrame(pi = 0.5, shannon = 5.0)
    res <- checkDiversityPlausibility(div, pathogen = "sars_cov_2")
    expect_s4_class(res, "DataFrame")
    expect_true("in_range" %in% colnames(res))
    ## pi = 0.5 is way above SARS-CoV-2 range [1e-6, 1e-3]
    pi_row <- res[res$metric == "pi", ]
    expect_false(pi_row$in_range)
    ## shannon = 5.0 is within [0, 15]
    sh_row <- res[res$metric == "shannon", ]
    expect_true(sh_row$in_range)
})

test_that("checkDiversityPlausibility accepts long-format input", {
    div <- data.frame(
        metric = c("pi", "shannon"),
        value  = c(1e-4, 10.0),
        stringsAsFactors = FALSE
    )
    res <- checkDiversityPlausibility(div, pathogen = "influenza")
    expect_equal(nrow(res), 2L)
    expect_true(all(res$in_range))
    expect_true(all(is.na(res$sample_id)))
})

test_that("checkDiversityPlausibility checks every sample (was: first row)", {
    div <- S4Vectors::DataFrame(sample_id = c("a", "b", "c"),
        pi = c(1e-4, 0.3, 2e-5), shannon = c(5, 7.5, 40))
    res <- checkDiversityPlausibility(div, pathogen = "sars_cov_2")
    expect_equal(nrow(res), 6L)
    expect_identical(res$sample_id, rep(c("a", "b", "c"), 2L))
    expect_identical(res$in_range[res$metric == "pi"], c(TRUE, FALSE, TRUE))
    expect_identical(res$in_range[res$metric == "shannon"],
                     c(TRUE, TRUE, FALSE))
    ## no sample_id column: rows are numbered
    res2 <- checkDiversityPlausibility(div[, c("pi", "shannon")],
                                       pathogen = "sars_cov_2")
    expect_identical(res2$sample_id, rep(c("1", "2", "3"), 2L))
    ## NA is never in range
    div$pi[1] <- NA
    expect_false(checkDiversityPlausibility(div, "sars_cov_2")$in_range[1])
})

test_that("the ranges are the heuristics the man page lists", {
    res <- checkDiversityPlausibility(
        S4Vectors::DataFrame(pi = 1e-4, shannon = 1), pathogen = "tb")
    expect_equal(res$expected_min, c(1e-9, 0))
    expect_equal(res$expected_max, c(1e-4, 5))
    expect_error(validateDiversity(S4Vectors::DataFrame(pi = 1e-4)),
                 "could not find function")
})

test_that("benchmarkCallers computes TP/FP/FN with truth", {
    gr_truth <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300), width = 1))
    S4Vectors::mcols(gr_truth)$ref <- c("A", "C", "G")
    S4Vectors::mcols(gr_truth)$alt <- c("T", "G", "A")
    truth <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr_truth,
        colData = S4Vectors::DataFrame(sample_id = "S1"))

    ## Caller detects 2 of 3 truth sites plus 1 false positive
    gr_call <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 400), width = 1))
    S4Vectors::mcols(gr_call)$ref <- c("A", "C", "T")
    S4Vectors::mcols(gr_call)$alt <- c("T", "G", "C")
    caller <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.4), ncol = 1)),
        rowRanges = gr_call,
        colData = S4Vectors::DataFrame(sample_id = "S1"))

    res <- benchmarkCallers(list(mycaller = caller), truth = truth)
    expect_s4_class(res, "DataFrame")
    expect_equal(res$TP, 2L)
    expect_equal(res$FP, 1L)
    expect_equal(res$FN, 1L)
    expect_equal(res$precision, 2 / 3, tolerance = 1e-10)
    expect_equal(res$recall, 2 / 3, tolerance = 1e-10)
})

# ---- Multi-pathogen validation ----

test_that("the influenza fixture is in range", {
    flu <- system.file("extdata", "test_influenza_donor.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(flu, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "flu_donor"))
    GenomeInfoDb::seqlengths(whe) <- c(HA = 1701L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- checkDiversityPlausibility(div, pathogen = "influenza")
    expect_true(all(res$in_range))
})

test_that("the HIV fixture is in range", {
    hiv <- system.file("extdata", "test_hiv_sample.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(hiv, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "hiv_sample"))
    GenomeInfoDb::seqlengths(whe) <- c(pol = 1000L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- checkDiversityPlausibility(div, pathogen = "hiv")
    expect_true(all(res$in_range))
})

test_that("the TB fixture is in range", {
    tb <- system.file("extdata", "test_tb_sample.tsv",
                      package = "WithinHostExperiment")
    whe <- readWithinHostTable(tb, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "tb_sample"))
    GenomeInfoDb::seqlengths(whe) <- c(Rv0001 = 4411532L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- checkDiversityPlausibility(div, pathogen = "tb")
    expect_true(all(res$in_range))
})

test_that("calcNeutralityTests works on influenza data", {
    flu <- system.file("extdata", "test_influenza_donor.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(flu, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "flu"))
    nt <- calcNeutralityTests(whe)
    expect_s4_class(nt, "DataFrame")
    expect_true("tajimaD" %in% colnames(nt))
    expect_true("fusFs" %in% colnames(nt))
    expect_true(is.finite(nt$tajimaD[1]))
})

test_that("exactBottleneck presence_absence works on test pair", {
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
    nb <- exactBottleneck(whe, pairId = "p1", maxNb = 50L,
                          method = "presence_absence")
    expect_true(is.finite(nb$Nb))
    expect_equal(nb$method, "presence_absence")
    expect_true(nb$Nb >= 1)
})

test_that("benchmarkCallers falls back to pairwise without truth", {
    gr1 <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300), width = 1))
    S4Vectors::mcols(gr1)$ref <- c("A", "C", "G")
    S4Vectors::mcols(gr1)$alt <- c("T", "G", "A")
    whe1 <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr1,
        colData = S4Vectors::DataFrame(sample_id = "S1"))

    gr2 <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 400), width = 1))
    S4Vectors::mcols(gr2)$ref <- c("A", "C", "T")
    S4Vectors::mcols(gr2)$alt <- c("T", "G", "C")
    whe2 <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.4), ncol = 1)),
        rowRanges = gr2,
        colData = S4Vectors::DataFrame(sample_id = "S1"))

    res <- benchmarkCallers(list(ivar = whe1, lofreq = whe2))
    expect_s4_class(res, "DataFrame")
    expect_true("jaccard" %in% colnames(res))
})

test_that("a recorded frequency of 0 is not a call (was: counted, inflating FP and FN)", {
    site_whe <- function(pos, freq) {
        gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(pos, width = 1))
        S4Vectors::mcols(gr)$ref <- rep("A", length(pos))
        S4Vectors::mcols(gr)$alt <- rep("T", length(pos))
        WithinHostExperiment(
            assays = list(altFreq = matrix(freq, ncol = 1)),
            rowRanges = gr,
            colData = S4Vectors::DataFrame(sample_id = "S1"))
    }
    truth <- site_whe(c(100, 200, 300), c(0.11, 0.19, 0.31))
    ## two sites assessed without the allele, as a generic table records
    caller <- site_whe(c(100, 200, 300, 400, 500),
                       c(0.11, 0.19, 0.31, 0, 0))
    res <- benchmarkCallers(list(c1 = caller), truth = truth)
    expect_equal(c(res$n_total, res$TP, res$FP, res$FN), c(3L, 3L, 0L, 0L))
    expect_equal(c(res$precision, res$recall, res$F1), c(1, 1, 1))
    ## a truth set that lists a known negative at 0
    truth0 <- site_whe(c(100, 200, 300, 400), c(0.11, 0.19, 0.31, 0))
    res0 <- benchmarkCallers(list(c1 = truth), truth = truth0)
    expect_equal(c(res0$TP, res0$FP, res0$FN), c(3L, 0L, 0L))
    ## without a truth set
    conc <- benchmarkCallers(list(a = caller, b = truth))
    expect_equal(c(conc$n_shared, conc$n_only1, conc$n_only2), c(3L, 0L, 0L))
    expect_equal(conc$jaccard, 1)
})
