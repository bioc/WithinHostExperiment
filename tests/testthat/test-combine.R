# tests/testthat/test-combine.R

library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

# ---- combineRows ----

test_that("combineRows combines WHE objects by rows", {
    data(example_whe, package = "WithinHostExperiment")
    n <- nrow(example_whe)
    mid <- floor(n / 2)
    whe1 <- example_whe[seq_len(mid), ]
    whe2 <- example_whe[seq(mid + 1L, n), ]
    merged <- combineRows(whe1, whe2)
    expect_s4_class(merged, "WithinHostExperiment")
    expect_equal(nrow(merged), n)
    expect_equal(ncol(merged), ncol(example_whe))
})

test_that("combineRows merges qcLogs", {
    data(example_whe, package = "WithinHostExperiment")
    whe1 <- flagISNV(example_whe[1:5, ], ISNVFilter())
    whe2 <- flagISNV(example_whe[6:10, ], ISNVFilter(minDepth = 500L))
    merged <- combineRows(whe1, whe2)
    expect_equal(nrow(qcLog(merged)),
                 nrow(qcLog(whe1)) + nrow(qcLog(whe2)))
})

test_that("combineRows merges qcFilters", {
    data(example_whe, package = "WithinHostExperiment")
    whe1 <- flagISNV(example_whe[1:5, ], ISNVFilter())
    whe2 <- flagISNV(example_whe[6:10, ], ISNVFilter(minDepth = 500L))
    merged <- combineRows(whe1, whe2)
    expect_equal(length(qcFilters(merged)),
                 length(qcFilters(whe1)) + length(qcFilters(whe2)))
})

test_that("combineRows() matches samples by name (was: bound by position under the first object's names)", {
    data(example_whe, package = "WithinHostExperiment")
    a <- example_whe[1:5, c(1, 2)]
    b <- example_whe[6:10, c(2, 1)]
    merged <- combineRows(a, b)
    expect_identical(colnames(merged), colnames(a))
    expect_identical(unname(assay(merged, "altFreq")),
                     unname(assay(example_whe[1:10, 1:2], "altFreq")))
    expect_identical(unname(assay(merged, "qcPass")),
                     unname(assay(example_whe[1:10, 1:2], "qcPass")))
    ## other samples are an error that names them
    other <- example_whe[6:10, 1:2]
    colnames(other) <- c("X1", "X2")
    expect_error(combineRows(a, other), paste0(
        "argument 2 lacks 2 of the first object's samples \\(donor, ",
        "recipient\\) and has 2 samples the first object lacks \\(X1, X2\\)"))
    expect_error(combineRows(a, example_whe[6:10, 1]),
                 "argument 2 lacks 1 of the first object's samples \\(recipient\\)")
})

# ---- combineCols ----

test_that("combineCols combines WHE objects by columns", {
    data(example_whe, package = "WithinHostExperiment")
    whe1 <- example_whe[, 1]
    whe2 <- example_whe[, 2]
    merged <- combineCols(whe1, whe2)
    expect_s4_class(merged, "WithinHostExperiment")
    expect_equal(ncol(merged), 2L)
    expect_equal(nrow(merged), nrow(example_whe))
})

test_that("combineCols merges qcLogs from both objects", {
    data(example_whe, package = "WithinHostExperiment")
    whe1 <- flagISNV(example_whe[, 1], ISNVFilter())
    whe2 <- flagISNV(example_whe[, 2], ISNVFilter(minDepth = 500L))
    merged <- combineCols(whe1, whe2)
    expect_true(nrow(qcLog(merged)) >= 2L)
})

# ---- VRanges -> WHE ----

test_that("as(VRanges, 'WithinHostExperiment') works", {
    skip_if_not_installed("VariantAnnotation")
    vr <- VariantAnnotation::VRanges(
        seqnames = c("seg", "seg", "seg"),
        ranges = IRanges::IRanges(c(100, 200, 100), width = 1),
        ref = c("A", "T", "A"),
        alt = c("G", "C", "G"),
        totalDepth = c(1000L, 2000L, 1500L),
        altDepth = c(50L, 300L, 75L),
        refDepth = c(950L, 1700L, 1425L),
        sampleNames = c("S1", "S1", "S2"))
    whe <- as(vr, "WithinHostExperiment")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 2L)
    expect_equal(ncol(whe), 2L)
    expect_true("altFreq" %in% assayNames(whe))
    freq <- assay(whe, "altFreq")
    expect_equal(unname(freq[1, "S1"]), 50 / 1000, tolerance = 1e-6)
})

test_that("VRanges roundtrip preserves data", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe_orig <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")

    vr <- as(whe_orig, "VRanges")
    whe_back <- as(vr, "WithinHostExperiment")
    expect_equal(nrow(whe_back), nrow(whe_orig))
    expect_equal(ncol(whe_back), ncol(whe_orig))
})

test_that("VRanges -> WHE preserves depths", {
    skip_if_not_installed("VariantAnnotation")
    vr <- VariantAnnotation::VRanges(
        seqnames = "seg",
        ranges = IRanges::IRanges(100, width = 1),
        ref = "A", alt = "G",
        totalDepth = 5000L,
        altDepth = 250L,
        refDepth = 4750L,
        sampleNames = "S1")
    whe <- as(vr, "WithinHostExperiment")
    expect_equal(unname(assay(whe, "totalDepth")[1, 1]), 5000L)
    expect_equal(unname(assay(whe, "altCount")[1, 1]), 250L)
    expect_equal(unname(assay(whe, "refCount")[1, 1]), 4750L)
})

test_that("combine methods dispatch through the S4Vectors generics", {
    data(example_whe, package = "WithinHostExperiment")
    whe1 <- flagISNV(example_whe[1:5, ], ISNVFilter())
    whe2 <- flagISNV(example_whe[6:10, ], ISNVFilter(minDepth = 500L))
    merged <- S4Vectors::combineRows(whe1, whe2)
    expect_s4_class(merged, "WithinHostExperiment")
    expect_equal(nrow(qcLog(merged)),
                 nrow(qcLog(whe1)) + nrow(qcLog(whe2)))
    merged_cols <- S4Vectors::combineCols(example_whe[, 1],
                                          example_whe[, 2])
    expect_s4_class(merged_cols, "WithinHostExperiment")
    expect_equal(ncol(merged_cols), 2L)
})
