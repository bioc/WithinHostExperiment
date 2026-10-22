test_that("validateDiversity flags out-of-range pi", {
    div <- S4Vectors::DataFrame(pi = 0.5, shannon = 5.0)
    res <- validateDiversity(div, pathogen = "sars_cov_2")
    expect_s4_class(res, "DataFrame")
    expect_true("in_range" %in% colnames(res))
    ## pi = 0.5 is way above SARS-CoV-2 range [1e-6, 1e-3]
    pi_row <- res[res$metric == "pi", ]
    expect_false(pi_row$in_range)
    ## shannon = 5.0 is within [0, 15]
    sh_row <- res[res$metric == "shannon", ]
    expect_true(sh_row$in_range)
})

test_that("validateDiversity accepts long-format input", {
    div <- data.frame(
        metric = c("pi", "shannon"),
        value  = c(1e-4, 10.0),
        stringsAsFactors = FALSE
    )
    res <- validateDiversity(div, pathogen = "influenza")
    expect_equal(nrow(res), 2L)
    expect_true(all(res$in_range))
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

test_that("validateDiversity passes for influenza-range diversity", {
    flu <- system.file("extdata", "test_influenza_donor.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(flu, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "flu_donor"))
    GenomeInfoDb::seqlengths(whe) <- c(HA = 1701L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- validateDiversity(div, pathogen = "influenza")
    expect_true(all(res$in_range))
})

test_that("validateDiversity passes for HIV-range diversity", {
    hiv <- system.file("extdata", "test_hiv_sample.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(hiv, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "hiv_sample"))
    GenomeInfoDb::seqlengths(whe) <- c(pol = 1000L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- validateDiversity(div, pathogen = "hiv")
    expect_true(all(res$in_range))
})

test_that("validateDiversity passes for TB-range diversity", {
    tb <- system.file("extdata", "test_tb_sample.tsv",
                      package = "WithinHostExperiment")
    whe <- readWithinHostTable(tb, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "tb_sample"))
    GenomeInfoDb::seqlengths(whe) <- c(Rv0001 = 4411532L)
    div <- calcDiversity(whe, indices = c("pi", "shannon"))
    res <- validateDiversity(div, pathogen = "tb")
    expect_true(all(res$in_range))
})

test_that("calcNeutralityTests works on influenza data", {
    flu <- system.file("extdata", "test_influenza_donor.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(flu, format = "ivar",
        colData = S4Vectors::DataFrame(sample_id = "flu"))
    nt <- calcNeutralityTests(whe, genomeLength = 1701L)
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
