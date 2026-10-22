library(GenomicRanges)

test_that("flagPrimerSites flags overlapping sites", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    # Primers at pos 40-65 and 880-950
    # Donor has pos=45 (in primer F1) and pos=889,945 (in primer R1)
    primers <- GRanges("test_segment", IRanges::IRanges(
        start = c(40, 880), end = c(65, 950)))
    result <- flagPrimerSites(whe, primers, action = "flag")
    qc <- assay(result, "qcPass")[, 1]
    expect_true("primer_filter" %in% qcLog(result)$step)
    expect_true(any(!qc))
})

test_that("flagPrimerSites remove action drops rows", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    primers <- GRanges("test_segment", IRanges::IRanges(
        start = c(40, 880), end = c(65, 950)))
    result <- flagPrimerSites(whe, primers, action = "remove")
    expect_true(nrow(result) < nrow(whe))
})

test_that("flagPrimerSites stops when primers and sites share no sequence name (was: nothing flagged, silently)", {
    ## a site at 45 on test_segment, inside a primer if the names matched
    gr <- GenomicRanges::GRanges("test_segment",
                                 IRanges::IRanges(45, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(0.35, 1, 1),
                      qcPass = matrix(TRUE, 1, 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "donor"))
    primers <- GRanges("NC_045512.2", IRanges::IRanges(start = 40, end = 65))
    expect_error(flagPrimerSites(whe, primers),
                 "share no sequence names, so no site could be flagged")
    expect_error(flagPrimerSites(whe, primers), "test_segment")
})

test_that("flagPrimerSites flags nothing for an empty primer set (was: a naming error)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20), width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1),
                      qcPass = matrix(TRUE, 2, 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    expect_message(out <- flagPrimerSites(whe, GRanges()), "No variant sites")
    expect_true(all(SummarizedExperiment::assay(out, "qcPass")))
})

test_that("flagPrimerSites handles no overlap", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    primers <- GRanges("test_segment", IRanges::IRanges(start = 1, end = 5))
    expect_message(
        result <- flagPrimerSites(whe, primers),
        "No variant sites"
    )
    expect_true(all(assay(result, "qcPass")))
})
