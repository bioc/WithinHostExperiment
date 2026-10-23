.make_pair_whe <- function() {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role      = c("donor", "recipient"),
            pair_id   = c("p1", "p1")),
        caller = "ivar")
}

test_that("exportPairFrequencies returns correct structure", {
    whe <- .make_pair_whe()
    pf  <- exportPairFrequencies(whe, pairId = "p1")
    expect_s4_class(pf, "DataFrame")
    expect_true(all(c("donor_freq", "recipient_freq", "shared") %in%
                        colnames(pf)))
    expect_true(nrow(pf) > 0L)
})

test_that("exportPairFrequencies identifies shared variants", {
    whe <- .make_pair_whe()
    pf  <- exportPairFrequencies(whe)
    expect_equal(sum(pf$shared), 5L)
})

test_that("asViralBottleneckInput returns data.frame", {
    whe <- .make_pair_whe()
    vb  <- asViralBottleneckInput(whe, pairId = "p1")
    expect_true(is.data.frame(vb))
    expect_true(all(c("pos", "donor_freq", "recipient_freq") %in%
                        colnames(vb)))
    expect_true(all(vb$donor_freq >= 0.03))
    expect_true(!any(is.na(vb$recipient_freq)))
})

test_that("calcSharedVariants returns correct counts", {
    whe <- .make_pair_whe()
    sv  <- calcSharedVariants(whe, pairId = "p1")
    expect_equal(sv$n_shared, 5L)
    expect_equal(sv$n_donor_only, 10L)
    expect_equal(sv$n_recipient_only, 5L)
    expect_equal(sv$n_total, 20L)
    ## read from VCFs, the recipient has no depth at the sites it did
    ## not call, so every donor-only site is of unknown depth
    expect_equal(sv$n_donor_only_unknown_depth, 10L)
    expect_equal(sv$n_donor_only_adequate_depth, 0L)
    expect_equal(sv$n_donor_only_low_depth, 0L)
    expect_equal(sv$n_donor_only_recipient_failed_qc, 0L)
    expect_equal(sv$n_donor_only_adequate_depth + sv$n_donor_only_low_depth +
                     sv$n_donor_only_unknown_depth +
                     sv$n_donor_only_recipient_failed_qc, sv$n_donor_only)
    ## no filter applied: the ISNVFilter default depth is the threshold
    expect_identical(sv$min_depth, 100L)
})

test_that("calcSharedVariants splits donor-only sites by recipient depth", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:5 * 10, width = 1))
    S4Vectors::mcols(gr)$ref <- rep("A", 5)
    S4Vectors::mcols(gr)$alt <- rep("G", 5)
    ## sites 1-4 donor-only, site 5 shared; the recipient's depth at the
    ## donor-only sites is 2000, 50, NA and 100
    freq <- cbind(d = c(0.2, 0.2, 0.2, 0.2, 0.3),
                  r = c(NA, NA, NA, NA, 0.4))
    depth <- cbind(d = rep(1000L, 5),
                   r = c(2000L, 50L, NA, 100L, 1500L))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
            role = c("donor", "recipient"), pair_id = c("p", "p"),
            row.names = c("d", "r")))
    sv <- calcSharedVariants(whe, pairId = "p")
    expect_equal(sv$n_donor_only, 4L)
    expect_equal(sv$n_donor_only_adequate_depth, 2L)
    expect_equal(sv$n_donor_only_low_depth, 1L)
    expect_equal(sv$n_donor_only_unknown_depth, 1L)
    expect_equal(sv$n_shared, 1L)
    ## the strictest applied minDepth sets the threshold
    whe_qc <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    sv_qc <- calcSharedVariants(whe_qc, pairId = "p")
    expect_identical(sv_qc$min_depth, 500L)
    expect_equal(sv_qc$n_donor_only_adequate_depth, 1L)
    expect_equal(sv_qc$n_donor_only_low_depth, 2L)
    ## an explicit threshold overrides it
    sv_x <- calcSharedVariants(whe_qc, pairId = "p", minDepth = 10L)
    expect_equal(sv_x$n_donor_only_adequate_depth, 3L)
    expect_error(calcSharedVariants(whe, pairId = "p", minDepth = -1),
                 "non-negative")
})


test_that("a recipient call removed by QC is not an absence at adequate depth (was: counted as one)", {
    gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(c(10, 20), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "A")
    S4Vectors::mcols(gr)$alt <- c("G", "G")
    ## the donor carries both at 20%; the recipient carries the first at
    ## 60%, which maxFreq = 0.5 removes, and has no call at the second
    freq <- cbind(d = c(0.2, 0.2), r = c(0.6, NA))
    depth <- cbind(d = c(1000L, 1000L), r = c(5000L, NA))
    alt <- cbind(d = c(200L, 200L), r = c(3000L, NA))
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth, altCount = alt),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
            role = c("donor", "recipient"), pair_id = c("p", "p"),
            row.names = c("d", "r")))
    whe <- flagISNV(whe, ISNVFilter())
    sv <- calcSharedVariants(whe, pairId = "p")
    expect_equal(sv$n_donor_only, 2L)
    expect_equal(sv$n_donor_only_recipient_failed_qc, 1L)
    expect_equal(sv$n_donor_only_adequate_depth, 0L)
    expect_equal(sv$n_donor_only_unknown_depth, 1L)
})

## A pair read from two generic CSV tables that list every assessed site:
## an assessed site without alternative reads is a row with aaf = 0
.pair_from_tables <- function() {
    write_table <- function(aaf, depth) {
        f <- tempfile(fileext = ".csv")
        utils::write.csv(data.frame(chrom = "seg", pos = 1:5 * 100L,
                                    ref = "A", alt = "G", aaf = aaf,
                                    depth = depth),
                         f, row.names = FALSE)
        f
    }
    ## site 100 is in both hosts, 200 and 300 are donor variants the
    ## recipient lacks at depths 1000 and 20, 400 is in the recipient
    ## alone, and 500 is in neither
    donor <- write_table(c(0.10, 0.20, 0.30, 0, 0), rep(1000L, 5))
    recipient <- write_table(c(0.15, 0, 0, 0.20, 0),
                             c(1000L, 1000L, 20L, 1000L, 1000L))
    readWithinHostTable(c(donor, recipient), format = "generic",
        colData = S4Vectors::DataFrame(sample_id = c("d", "r"),
            role = c("donor", "recipient"), pair_id = c("p", "p")))
}

test_that("a frequency of 0 is not a call in a pair (was: counted as shared)", {
    whe <- .pair_from_tables()
    pf <- exportPairFrequencies(whe, pairId = "p")
    ## the site carried by neither host is not listed
    expect_identical(pf$position, c(100L, 200L, 300L, 400L))
    expect_identical(pf$shared, c(TRUE, FALSE, FALSE, FALSE))
    ## the recorded zeros are kept as they are
    expect_identical(pf$recipient_freq, c(0.15, 0, 0, 0.20))
    for (x in list(whe, flagISNV(whe, ISNVFilter()))) {
        sv <- calcSharedVariants(x, pairId = "p")
        expect_identical(sv$n_shared, 1L)
        expect_identical(sv$n_donor_only, 2L)
        expect_identical(sv$n_donor_only_adequate_depth, 1L)
        expect_identical(sv$n_donor_only_low_depth, 1L)
        expect_identical(sv$n_donor_only_unknown_depth, 0L)
        expect_identical(sv$n_donor_only_recipient_failed_qc, 0L)
        expect_identical(sv$n_recipient_only, 1L)
        expect_identical(sv$n_total, 4L)
        expect_equal(sv$sharing_rate, 1 / 4)
    }
})

test_that("exportPairFrequencies errors without pairs", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    expect_error(exportPairFrequencies(whe), "No transmission pairs")
})
