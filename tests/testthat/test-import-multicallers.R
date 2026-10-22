# tests/testthat/test-import-multicallers.R
# Tests for LoFreq and Freebayes VCF format support

library(S4Vectors)
library(SummarizedExperiment)

# ---- LoFreq ----

test_that("readWithinHost parses LoFreq VCF correctly", {
    vcf <- system.file("extdata", "test_donor_lofreq.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "lofreq")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    freq <- assay(whe, "altFreq")[, 1]
    expect_true(all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("LoFreq frequencies match iVar for same data", {
    vcf_lf <- system.file("extdata", "test_donor_lofreq.vcf",
                          package = "WithinHostExperiment")
    vcf_iv <- system.file("extdata", "test_donor.vcf",
                          package = "WithinHostExperiment")
    whe_lf <- readWithinHost(vcf_lf,
        colData = DataFrame(sample_id = "d"), caller = "lofreq")
    whe_iv <- readWithinHost(vcf_iv,
        colData = DataFrame(sample_id = "d"), caller = "ivar")
    freq_lf <- sort(assay(whe_lf, "altFreq")[, 1])
    freq_iv <- sort(assay(whe_iv, "altFreq")[, 1])
    expect_equal(freq_lf, freq_iv, tolerance = 0.01)
})

test_that("readWithinHost auto-detects LoFreq", {
    vcf <- system.file("extdata", "test_donor_lofreq.vcf",
                       package = "WithinHostExperiment")
    expect_message(
        whe <- readWithinHost(vcf,
            colData = DataFrame(sample_id = "d"),
            caller = "auto"),
        "lofreq"
    )
    expect_s4_class(whe, "WithinHostExperiment")
})

# ---- Freebayes ----

test_that("readWithinHost parses Freebayes VCF correctly", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "freebayes")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    freq <- assay(whe, "altFreq")[, 1]
    expect_true(all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("readWithinHost auto-detects Freebayes", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    expect_message(
        whe <- readWithinHost(vcf,
            colData = DataFrame(sample_id = "d"),
            caller = "auto"),
        "freebayes"
    )
})

test_that("Freebayes depths are correct", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "freebayes")
    dp <- assay(whe, "totalDepth")[, 1]
    expect_true(all(dp > 0, na.rm = TRUE))
    expect_true(is.integer(dp))
})

# ---- Multi-allelic records ----

test_that("multi-allelic LoFreq records are split per allele", {
    vcf <- tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.0",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "seg\t100\t.\tA\tG,T\t100\tPASS\tDP=1000;AF=0.10,0.05",
        "seg\t200\t.\tC\tT\t100\tPASS\tDP=900;AF=0.2"), vcf)
    expect_no_warning(whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"), caller = "lofreq"))
    expect_equal(nrow(whe), 3L)
    mc <- mcols(rowRanges(whe))
    freq <- assay(whe, "altFreq")[, 1]
    expect_equal(freq[mc$alt == "G"], 0.10)
    expect_equal(freq[mc$alt == "T" & start(rowRanges(whe)) == 100], 0.05)
    expect_false(anyNA(freq))
})

test_that("multi-allelic Freebayes records are split per allele", {
    vcf <- tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
        "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"RO\">",
        "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"AO\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\td",
        "seg\t100\t.\tA\tG,T\t100\t.\tDP=1000\tGT:DP:RO:AO\t0/1/2:1000:850:100,50"),
        vcf)
    whe <- readWithinHost(vcf, colData = DataFrame(sample_id = "d"),
                          caller = "freebayes")
    expect_equal(nrow(whe), 2L)
    expect_equal(sort(unname(assay(whe, "altFreq")[, 1])), c(0.05, 0.10))
    expect_equal(sort(unname(assay(whe, "altCount")[, 1])), c(50L, 100L))
})
