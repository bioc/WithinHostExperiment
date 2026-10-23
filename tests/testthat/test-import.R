# tests/testthat/test-import.R

library(S4Vectors)

# ---- readWithinHost: single VCF ----

test_that("readWithinHost reads single iVar VCF", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")

    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    expect_equal(ncol(whe), 1L)
    expect_true("altFreq" %in% assayNames(whe))
    expect_true("totalDepth" %in% assayNames(whe))
    expect_true("refCount" %in% assayNames(whe))
    expect_true("altCount" %in% assayNames(whe))
    expect_true("qcPass" %in% assayNames(whe))
})

test_that("altFreq values are in [0, 1]", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    freq <- SummarizedExperiment::assay(whe, "altFreq")[, 1]
    expect_true(all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("totalDepth values are positive integers", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    dp <- SummarizedExperiment::assay(whe, "totalDepth")[, 1]
    expect_true(all(dp > 0, na.rm = TRUE))
    expect_true(is.integer(dp))
})

# ---- readWithinHost: multiple VCFs ----

test_that("readWithinHost reads multiple VCFs and merges positions", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(
            sample_id = c("donor", "recipient"),
            role      = c("donor", "recipient"),
            pair_id   = c("p1", "p1")),
        caller = "ivar")

    expect_equal(ncol(whe), 2L)
    # Union: donor has 15 unique, recipient has 10, 5 shared -> 20 total
    expect_equal(nrow(whe), 20L)
    # Shared sites have values in both columns
    freq <- SummarizedExperiment::assay(whe, "altFreq")
    both_present <- sum(!is.na(freq[, 1]) & !is.na(freq[, 2]))
    expect_equal(both_present, 5L)
})

test_that("non-shared sites have NA in the other sample", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    freq <- SummarizedExperiment::assay(whe, "altFreq")
    # Donor-only sites should be NA in recipient column
    donor_only <- is.na(freq[, 2]) & !is.na(freq[, 1])
    expect_equal(sum(donor_only), 10L)  # 15 - 5 shared
})

# ---- readWithinHost: auto-detect ----

test_that("caller auto-detection works for iVar VCF", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    expect_message(
        whe <- readWithinHost(vcf,
            colData = DataFrame(sample_id = "d"),
            caller = "auto"),
        "ivar"
    )
    expect_s4_class(whe, "WithinHostExperiment")
})

# ---- readWithinHost: colData metadata preserved ----

test_that("colData columns are preserved", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(
            sample_id = "donor", host_id = "P1",
            role = "donor", pair_id = "pair_1",
            ct_value = 18.5),
        caller = "ivar")
    cd <- SummarizedExperiment::colData(whe)
    expect_true("host_id" %in% colnames(cd))
    expect_true("ct_value" %in% colnames(cd))
    expect_equal(cd$ct_value, 18.5)
})

# ---- readWithinHost: rowRanges has ref/alt ----

test_that("rowRanges contains ref and alt mcols", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    rr <- SummarizedExperiment::rowRanges(whe)
    expect_true("ref" %in% colnames(S4Vectors::mcols(rr)))
    expect_true("alt" %in% colnames(S4Vectors::mcols(rr)))
})

# ---- readWithinHost: error handling ----

test_that("readWithinHost rejects non-existent files", {
    expect_error(
        readWithinHost("no_such_file.vcf",
                       colData = DataFrame(sample_id = "x")),
        "not found"
    )
})

test_that("readWithinHost rejects mismatched colData length", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    expect_error(
        readWithinHost(vcf,
                       colData = DataFrame(sample_id = c("a", "b"))),
        "must match"
    )
})

test_that("readWithinHost rejects colData without sample_id", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    expect_error(
        readWithinHost(vcf, colData = DataFrame(name = "x")),
        "sample_id"
    )
})

# ---- readWithinHostTable: iVar TSV ----

test_that("readWithinHostTable reads iVar TSV", {
    tsv <- system.file("extdata", "test_donor_ivar.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(tsv, format = "ivar")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    expect_true("altFreq" %in% assayNames(whe))
})

test_that("iVar's PASS and strand columns are not read, as the man page says", {
    tsv <- tempfile(fileext = ".tsv")
    writeLines(c(paste("REGION", "POS", "REF", "ALT", "REF_DP", "ALT_DP",
                       "ALT_FREQ", "TOTAL_DP", "PVAL", "PASS", "REF_RV",
                       "ALT_RV", sep = "\t"),
        paste("seg", 100, "A", "G", 900, 100, 0.1, 1000, 0.001, "TRUE",
              450, 0, sep = "\t"),
        paste("seg", 200, "C", "T", 960, 40, 0.04, 1000, 0.2, "FALSE",
              480, 20, sep = "\t"),
        paste("seg", 300, "A", "+TTG", 900, 100, 0.1, 1000, 0.001, "TRUE",
              450, 50, sep = "\t")), tsv)
    whe <- flagISNV(readWithinHostTable(tsv, format = "ivar"))
    ## the PASS FALSE row and the insertion are calls that pass QC
    expect_true(all(assay(whe, "qcPass")))
    expect_identical(BiocGenerics::width(rowRanges(whe)), rep(1L, 3L))
    expect_identical(rowRanges(whe)$alt[3], "+TTG")
    expect_false("strand_bias_filter" %in% qcLog(whe)$step)
    ## the documented recipe drops the PASS FALSE row
    calls <- read.delim(tsv, colClasses = "character")
    passed <- tempfile(fileext = ".tsv")
    write.table(calls[calls$PASS == "TRUE", ], passed, sep = "\t",
        quote = FALSE, row.names = FALSE)
    kept <- readWithinHostTable(passed, format = "ivar")
    expect_identical(BiocGenerics::start(rowRanges(kept)), c(100L, 300L))
    expect_identical(rowRanges(kept)$alt, c("G", "+TTG"))
})

test_that("readWithinHostTable() keeps an allele of T when the column holds only T (was: TRUE)", {
    ## read.delim() and read.csv() read a column of "T" as logical
    tsv <- tempfile(fileext = ".tsv")
    writeLines(c(paste("REGION", "POS", "REF", "ALT", "REF_DP", "ALT_DP",
                       "ALT_FREQ", "TOTAL_DP", sep = "\t"),
        paste("seg", 100, "C", "T", 900, 100, 0.1, 1000, sep = "\t"),
        paste("seg", 200, "C", "T", 960, 40, 0.04, 1000, sep = "\t")), tsv)
    rr <- rowRanges(readWithinHostTable(tsv, format = "ivar"))
    expect_identical(rr$alt, c("T", "T"))
    expect_identical(rr$ref, c("C", "C"))
    csv <- tempfile(fileext = ".csv")
    writeLines(c("chrom,pos,ref,alt,aaf,depth", "T,100,T,G,0.1,1000",
                 "T,200,T,A,0.2,1000"), csv)
    rr <- rowRanges(readWithinHostTable(csv, format = "generic"))
    expect_identical(rr$ref, c("T", "T"))
    expect_identical(as.character(GenomicRanges::seqnames(rr)), c("T", "T"))
})

# ---- readWithinHostTable: generic CSV ----

test_that("readWithinHostTable reads generic CSV", {
    csv <- system.file("extdata", "test_generic.csv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(csv, format = "generic")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
})

test_that("readWithinHostTable auto-generates sample_id from filename", {
    tsv <- system.file("extdata", "test_donor_ivar.tsv",
                       package = "WithinHostExperiment")
    whe <- readWithinHostTable(tsv, format = "ivar")
    cd <- SummarizedExperiment::colData(whe)
    expect_equal(as.character(cd$sample_id), "test_donor_ivar")
})

# ---- readWithinHostTable: error handling ----

test_that("readWithinHostTable rejects non-existent files", {
    expect_error(
        readWithinHostTable("nope.tsv", format = "ivar"),
        "not found"
    )
})

# ---- callerInfo in metadata ----

test_that("caller is stored in metadata", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    ci <- S4Vectors::metadata(whe)$calling_info
    expect_s4_class(ci, "DataFrame")
    expect_equal(ci$caller, "ivar")
})
