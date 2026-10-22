library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

# ---- as(whe, "data.frame") ----

test_that("coerce WHE to data.frame produces long format", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    df <- as(whe, "data.frame")
    expect_true(is.data.frame(df))
    expect_true(all(c("chrom", "position", "ref", "alt",
                      "sample_id", "altFreq") %in% colnames(df)))
    expect_equal(nrow(df), 15L)
    expect_true(all(df$altFreq >= 0 & df$altFreq <= 1))
})

test_that("coerce multi-sample WHE to data.frame", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    df <- as(whe, "data.frame")
    expect_true(is.data.frame(df))
    expect_equal(nrow(df), 25L)
    expect_true(all(c("D", "R") %in% df$sample_id))
})

test_that("coerce WHE to data.frame includes qcPass", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    df <- as(whe, "data.frame")
    expect_true("qcPass" %in% colnames(df))
    expect_true(any(!df$qcPass))
})

# ---- as(whe, "VRanges") ----

test_that("coerce WHE to VRanges", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    vr <- as(whe, "VRanges")
    expect_s4_class(vr, "VRanges")
    expect_equal(length(vr), 15L)
    expect_true(all(!is.na(VariantAnnotation::ref(vr))))
    expect_true(all(!is.na(VariantAnnotation::alt(vr))))
})

test_that("coerce multi-sample WHE to VRanges preserves sample names", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    vr <- as(whe, "VRanges")
    expect_s4_class(vr, "VRanges")
    expect_true(all(c("D", "R") %in% sampleNames(vr)))
    expect_equal(length(vr), 25L)
})

test_that("VRanges coercion keeps every recorded sequence (was: lengths dropped)", {
    gr <- GRanges(c("seg1", "seg1"), IRanges::IRanges(c(10, 20), width = 1),
                  ref = c("A", "C"), alt = c("G", "T"))
    GenomeInfoDb::seqlevels(gr) <- c("seg1", "seg2")
    GenomeInfoDb::seqlengths(gr) <- c(seg1 = 600L, seg2 = 400L)
    freq <- matrix(c(0.2, 0.3, 0.1, NA), ncol = 2,
                   dimnames = list(NULL, c("a", "b")))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq,
                      totalDepth = matrix(1000L, 2, 2, dimnames = dimnames(freq)),
                      altCount = matrix(c(200L, 300L, 100L, NA), 2, 2,
                                        dimnames = dimnames(freq))),
        rowRanges = gr,
        colData = DataFrame(sample_id = c("a", "b"), row.names = c("a", "b")))
    vr <- as(whe, "VRanges")
    expect_identical(GenomeInfoDb::seqlengths(vr), c(seg1 = 600L, seg2 = 400L))
    back <- as(vr, "WithinHostExperiment")
    expect_identical(GenomeInfoDb::seqlengths(back),
                     c(seg1 = 600L, seg2 = 400L))
    ## so the genome length, and pi, survive the round trip
    expect_equal(calcDiversity(back, indices = "pi")$pi,
                 calcDiversity(whe, indices = "pi")$pi)
})

test_that("VRanges coercion accepts a sequence name holding ':' (was: split on it)", {
    vr <- VariantAnnotation::VRanges(
        seqnames = c("chr:A", "chr:A"),
        ranges = IRanges::IRanges(c(10, 20), width = 1),
        ref = c("A", "C"), alt = c("G", "T"),
        totalDepth = c(1000L, 1000L), altDepth = c(100L, 200L),
        refDepth = c(900L, 800L), sampleNames = c("s1", "s1"))
    whe <- as(vr, "WithinHostExperiment")
    rr <- SummarizedExperiment::rowRanges(whe)
    expect_identical(as.character(GenomeInfoDb::seqnames(rr)),
                     c("chr:A", "chr:A"))
    expect_identical(BiocGenerics::start(rr), c(10L, 20L))
    expect_identical(mcols(rr)$alt, c("G", "T"))
})
