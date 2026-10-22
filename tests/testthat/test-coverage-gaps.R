# tests/testthat/test-coverage-gaps.R
# Targeted tests for coverage gaps

library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

# ---- .detectCaller: multiple headers ----

test_that(".detectCaller identifies lofreq from header", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##source=lofreq call",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ), tmp)
    result <- WithinHostExperiment:::.detectCaller(tmp)
    expect_equal(result, "lofreq")
    unlink(tmp)
})

test_that(".detectCaller identifies freebayes from header", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##source=freeBayes v1.3.6",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ), tmp)
    result <- WithinHostExperiment:::.detectCaller(tmp)
    expect_equal(result, "freebayes")
    unlink(tmp)
})

test_that(".detectCaller returns generic for unknown caller", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##source=MyCaller v9.9",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ), tmp)
    result <- WithinHostExperiment:::.detectCaller(tmp)
    expect_equal(result, "generic")
    unlink(tmp)
})

# ---- compareDiversity: edge cases ----

test_that("compareDiversity handles single group gracefully", {
    gr <- GRanges("s", IRanges::IRanges(1:5, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(runif(10), nrow = 5, ncol = 2)),
        rowRanges = gr,
        colData   = DataFrame(
            sample_id = c("A", "B"),
            group     = c("g1", "g1")))
    expect_message(
        compareDiversity(whe, groupBy = "group"),
        "2 groups"
    )
})

test_that("compareDiversity works with 2 groups (insufficient n)", {
    gr <- GRanges("s", IRanges::IRanges(seq_len(10) * 50, width = 1))
    freq <- matrix(runif(20, 0.03, 0.45), nrow = 10, ncol = 2)
    whe <- WithinHostExperiment(
        assays    = list(altFreq = freq),
        rowRanges = gr,
        colData   = DataFrame(
            sample_id = c("A", "B"),
            group     = c("g1", "g2")))
    result <- compareDiversity(whe, groupBy = "group")
    expect_true(is(result, "DataFrame"))
    expect_true("note" %in% colnames(result))
})

# ---- flagISNV: strand bias path ----

test_that("flagISNV applies strand bias filter when assays present", {
    gr <- GRanges("s", IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays = list(
            altFreq     = matrix(c(0.1, 0.2, 0.3), ncol = 1),
            totalDepth  = matrix(c(1000L, 1000L, 1000L), ncol = 1),
            altCount    = matrix(c(100L, 200L, 300L), ncol = 1),
            fwdAltCount = matrix(c(50L, 190L, 150L), ncol = 1),
            revAltCount = matrix(c(50L, 10L, 150L), ncol = 1)
        ),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "S1"))
    result <- flagISNV(whe, ISNVFilter(maxStrandBias = 10))
    qc <- assay(result, "qcPass")[, 1]
    expect_false(qc[2])
    expect_true(qc[1])
    expect_true(qc[3])
    expect_true("strand_bias_filter" %in% qcLog(result)$step)
})

# ---- transmissionPairs: multiple pairs ----

test_that("transmissionPairs handles multiple pairs", {
    gr <- GRanges("s", IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(runif(12), nrow = 3, ncol = 4)),
        rowRanges = gr,
        colData   = DataFrame(
            sample_id = c("D1", "R1", "D2", "R2"),
            role      = c("donor", "recipient", "donor", "recipient"),
            pair_id   = c("p1", "p1", "p2", "p2")))
    pairs <- transmissionPairs(whe)
    expect_equal(nrow(pairs), 2L)
    expect_true("p1" %in% pairs$pair_id)
    expect_true("p2" %in% pairs$pair_id)
})

# ---- WithinHostExperiment: metadata parameter ----

test_that("WithinHostExperiment stores metadata", {
    gr <- GRanges("s", IRanges::IRanges(1, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(0.1, ncol = 1)),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "S1"),
        metadata  = list(reference_genome = "test_ref.fa",
                         genome_length = 1000L))
    expect_equal(metadata(whe)$reference_genome, "test_ref.fa")
    expect_equal(metadata(whe)$genome_length, 1000L)
})

# ---- asViralBottleneckInput: no shared ----

test_that("asViralBottleneckInput handles no shared variants", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_independent.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(
            sample_id = c("d", "ind"),
            role      = c("donor", "recipient"),
            pair_id   = c("px", "px")),
        caller = "ivar")
    vb <- asViralBottleneckInput(whe, pairId = "px")
    expect_true(is.data.frame(vb))
    expect_true(all(vb$recipient_freq == 0))
})

# ---- readWithinHost: empty VCF ----

test_that("readWithinHost handles VCF with no data lines", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##source=iVar",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ), tmp)
    whe <- readWithinHost(tmp,
        colData = DataFrame(sample_id = "empty"),
        caller = "ivar")
    expect_equal(nrow(whe), 0L)
    expect_equal(ncol(whe), 1L)
    unlink(tmp)
})

# ---- passedISNV: all flagged ----

test_that("passedISNV returns 0 rows when all flagged", {
    gr <- GRanges("s", IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays = list(
            altFreq    = matrix(c(0.01, 0.01, 0.01), ncol = 1),
            totalDepth = matrix(c(50L, 50L, 50L), ncol = 1),
            altCount   = matrix(c(1L, 1L, 1L), ncol = 1)),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "S1"))
    whe    <- flagISNV(whe, ISNVFilter(minDepth = 1000L))
    passed <- passedISNV(whe)
    expect_equal(nrow(passed), 0L)
})

# ---- qcSummary: multiple samples ----

test_that("qcSummary works with multiple samples", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter(minDepth = 200L))
    s   <- qcSummary(whe)
    expect_equal(nrow(s), 2L)
    expect_true(all(c("sample_id", "total_sites", "passed_sites") %in%
                        colnames(s)))
})

# ---- calcDiversity: usePassedOnly flag ----

test_that("calcDiversity usePassedOnly flag works", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))

    div_with <- calcDiversity(whe, indices = "richness",
                              usePassedOnly = TRUE)
    div_without <- calcDiversity(whe, indices = "richness",
                                 usePassedOnly = FALSE)
    expect_true(div_with$richness <= div_without$richness)
})

# ---- readWithinHost: genome parameter ----

test_that("readWithinHost records the genome's sequence lengths", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    ref <- Biostrings::readDNAStringSet(system.file("extdata",
        "test_ref.fa", package = "WithinHostExperiment"))
    whe <- readWithinHost(vcf, colData = DataFrame(sample_id = "d"),
                          caller = "ivar", genome = ref)
    expect_identical(GenomeInfoDb::seqlengths(rowRanges(whe))[["test_segment"]],
                     1000L)
    expect_error(readWithinHost(vcf, colData = DataFrame(sample_id = "d"),
                                caller = "ivar", genome = "my_ref.fa"),
                 "readDNAStringSet")
})

# ---- plotFrequencySpectrum: density type ----

test_that("plotFrequencySpectrum density type works", {
    skip_if_not_installed("ggplot2")
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    p <- plotFrequencySpectrum(whe, type = "density")
    expect_s3_class(p, "ggplot")
})

# ---- readWithinHostTable: colData validation ----

test_that("readWithinHostTable rejects mismatched colData", {
    tsv <- system.file("extdata", "test_donor_ivar.tsv",
                       package = "WithinHostExperiment")
    expect_error(
        readWithinHostTable(tsv, format = "ivar",
                            colData = DataFrame(sample_id = c("a", "b"))),
        "must match"
    )
})

test_that("readWithinHostTable rejects colData without sample_id", {
    tsv <- system.file("extdata", "test_donor_ivar.tsv",
                       package = "WithinHostExperiment")
    expect_error(
        readWithinHostTable(tsv, format = "ivar",
                            colData = DataFrame(name = "x")),
        "sample_id"
    )
})

# ---- ISNVFilter: edge values ----

test_that("ISNVFilter accepts boundary values", {
    f <- ISNVFilter(minFreq = 0, maxFreq = 1, minDepth = 0L,
                    minAltReads = 0L)
    expect_s4_class(f, "ISNVFilter")
    expect_equal(f@minFreq, 0)
    expect_equal(f@maxFreq, 1)
})
