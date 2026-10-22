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
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    primers <- GRanges("NC_045512.2", IRanges::IRanges(start = 40, end = 65))
    expect_error(flagPrimerSites(whe, primers),
                 "share no sequence names, so no site could be flagged")
    expect_error(flagPrimerSites(whe, primers), "test_segment")
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
