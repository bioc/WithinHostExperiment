## Helper: load two caller WHEs from test VCFs
.load_caller_whes <- function() {
    vcf_ivar <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    vcf_lofreq <- system.file("extdata", "test_donor_lofreq.vcf",
        package = "WithinHostExperiment")
    vcf_fb <- system.file("extdata", "test_donor_freebayes.vcf",
        package = "WithinHostExperiment")

    whe_ivar <- readWithinHost(vcf_ivar,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    whe_lofreq <- readWithinHost(vcf_lofreq,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "lofreq")
    whe_fb <- readWithinHost(vcf_fb,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "freebayes")

    list(ivar = whe_ivar, lofreq = whe_lofreq, freebayes = whe_fb)
}

test_that("consensusISNV returns WHE with all shared sites from test VCFs", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes[c("ivar", "lofreq")], min_callers = 2L)
    expect_s4_class(result, "WithinHostExperiment")
    ## Both VCFs have the same 15 sites with identical alleles
    expect_equal(nrow(result), 15L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    expect_true(all(mc$n_callers == 2L))
    expect_true(all(grepl("ivar", mc$caller_names)))
    expect_true(all(grepl("lofreq", mc$caller_names)))
})

test_that("consensusISNV freq_spread is zero for identical-frequency callers", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes[c("ivar", "lofreq")], min_callers = 2L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    ## ivar and lofreq report the same AF values, so spread should be ~0
    expect_true(all(mc$freq_spread < 1e-6))
})

test_that("consensusISNV works with three callers and min_callers = 3", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes, min_callers = 3L)
    expect_s4_class(result, "WithinHostExperiment")
    ## All three callers share the same 15 sites
    expect_equal(nrow(result), 15L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    expect_true(all(mc$n_callers == 3L))
})

## Helper: one caller's object from sites and a frequency matrix
.caller_whe <- function(seqnames, pos, alt, freq, ref = rep("A", length(pos)),
                        seqinfo = NULL, colData = NULL) {
    freq <- as.matrix(freq)
    if (is.null(colnames(freq)))
        colnames(freq) <- paste0("s", seq_len(ncol(freq)))
    gr <- GenomicRanges::GRanges(seqnames, IRanges::IRanges(pos, width = 1L))
    if (!is.null(seqinfo)) {
        GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevels(seqinfo)
        GenomeInfoDb::seqinfo(gr) <- seqinfo
    }
    S4Vectors::mcols(gr)$ref <- ref
    S4Vectors::mcols(gr)$alt <- alt
    if (is.null(colData))
        colData <- S4Vectors::DataFrame(sample_id = colnames(freq))
    WithinHostExperiment(assays = list(altFreq = freq), rowRanges = gr,
                         colData = colData)
}

test_that("consensusISNV() counts a caller once when it lists a site twice (was: twice)", {
    ## caller a reports seg:10 T twice, caller b never reports it
    a <- .caller_whe("seg", c(10, 10, 20), c("T", "T", "G"),
                     c(0.10, 0.11, 0.30))
    b <- .caller_whe("seg", c(20, 10), c("G", "T"), c(0.31, NA))
    cons <- consensusISNV(list(a = a, b = b))
    expect_identical(GenomicRanges::start(
        SummarizedExperiment::rowRanges(cons)), 20L)
    ## with b reporting it too, the site has two callers, each named once
    b2 <- .caller_whe("seg", c(20, 10), c("G", "T"), c(0.31, 0.12))
    cons2 <- consensusISNV(list(a = a, b = b2))
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(cons2))
    expect_identical(mc$n_callers, c(2L, 2L))
    expect_identical(mc$caller_names, c("a,b", "a,b"))
    ## all three frequencies enter the mean and the spread
    expect_equal(SummarizedExperiment::assay(cons2, "altFreq")[, 1],
                 c(mean(c(0.10, 0.11, 0.12)), mean(c(0.30, 0.31))))
    expect_equal(mc$freq_spread, c(0.02, 0.01))
})

test_that("consensusISNV() keeps a sequence name containing ':' (was: cut at the colon)", {
    a <- .caller_whe("HLA-A*01:01", c(10, 20), c("T", "G"), c(0.10, 0.20))
    b <- .caller_whe("HLA-A*01:01", c(20, 10), c("G", "T"), c(0.21, 0.11))
    rr <- SummarizedExperiment::rowRanges(consensusISNV(list(a = a, b = b)))
    expect_identical(as.character(GenomeInfoDb::seqnames(rr)),
                     c("HLA-A*01:01", "HLA-A*01:01"))
    expect_identical(GenomicRanges::start(rr), c(10L, 20L))
    expect_identical(S4Vectors::mcols(rr)$alt, c("T", "G"))
})

test_that("consensusISNV() matches a worked three-caller example with gaps", {
    a <- .caller_whe("seg", c(10, 20, 30), c("T", "G", "C"),
                     cbind(s1 = c(0.10, 0.20, NA), s2 = c(NA, 0.40, 0.50)))
    b <- .caller_whe("seg", c(20, 10, 30), c("G", "T", "C"),
                     cbind(s1 = c(0.22, 0.11, 0.30), s2 = c(0.41, NA, 0.52)))
    cc <- .caller_whe("seg", c(30, 10), c("C", "T"),
                      cbind(s1 = c(NA, 0.12), s2 = c(0.51, 0.50)))
    cons <- consensusISNV(list(a = a, b = b, c = cc))
    rr <- SummarizedExperiment::rowRanges(cons)
    expect_identical(GenomicRanges::start(rr), c(10L, 20L, 30L))
    expect_equal(SummarizedExperiment::assay(cons, "altFreq"),
                 cbind(s1 = c(0.11, 0.21, NA), s2 = c(NA, 0.405, 0.51)))
    mc <- S4Vectors::mcols(rr)
    expect_identical(mc$ref, c("A", "A", "A"))
    expect_identical(mc$alt, c("T", "G", "C"))
    ## the site columns describe the site's last kept record: at 20 that
    ## is sample s2, where a and b differ by 0.01
    expect_identical(mc$n_callers, c(3L, 2L, 3L))
    expect_identical(mc$caller_names, c("a,b,c", "a,b", "a,b,c"))
    expect_equal(mc$freq_spread, c(0.02, 0.01, 0.02))
    ## 30 in s1 and 10 in s2 were seen by one caller only
    expect_error(consensusISNV(list(a = a, b = b, c = cc), min_callers = 4L),
                 "No sites passed consensus filters")
})

test_that("consensusISNV() handles 2,000 sites in 20 samples from two callers", {
    set.seed(1)
    n <- 2000L
    freq <- matrix(stats::runif(n * 20L, 0.05, 0.5), n, 20L,
                   dimnames = list(NULL, paste0("s", seq_len(20L))))
    a <- .caller_whe("seg", seq_len(n), rep("T", n), freq)
    b <- .caller_whe("seg", seq_len(n), rep("T", n), freq + 0.01)
    cons <- consensusISNV(list(a = a, b = b))
    expect_identical(dim(cons), c(n, 20L))
    expect_equal(unname(SummarizedExperiment::assay(cons, "altFreq")),
                 unname(freq) + 0.005)
})

test_that("consensusISNV() keeps the inputs' seqinfo (was: seqlengths dropped)", {
    si <- GenomeInfoDb::Seqinfo(c("seg1", "seg2"), c(1000L, 800L),
                                genome = "toy")
    a <- .caller_whe("seg1", c(10, 20), c("T", "G"), c(0.10, 0.20),
                     seqinfo = si)
    b <- .caller_whe("seg1", c(20, 10), c("G", "T"), c(0.21, 0.11),
                     seqinfo = si)
    cons <- consensusISNV(list(a = a, b = b))
    ## seg2 holds no consensus site but keeps its length
    expect_identical(GenomeInfoDb::seqinfo(cons), si)
    ## an input without lengths leaves the other input's lengths in place
    b0 <- .caller_whe("seg1", c(20, 10), c("G", "T"), c(0.21, 0.11))
    expect_identical(GenomeInfoDb::seqlengths(
        consensusISNV(list(a = a, b = b0)))[["seg1"]], 1000L)
    ## lengths that disagree stop the consensus
    b2 <- .caller_whe("seg1", c(20, 10), c("G", "T"), c(0.21, 0.11),
                      seqinfo = GenomeInfoDb::Seqinfo("seg1", 999L))
    expect_error(consensusISNV(list(a = a, b = b2)),
                 "disagree on their sequences.*seg1")
})

test_that("consensusISNV() keeps colData and every sample (was: sample_id only)", {
    data(example_whe, package = "WithinHostExperiment")
    cons <- consensusISNV(list(a = example_whe, b = example_whe))
    expect_identical(as.data.frame(SummarizedExperiment::colData(cons)),
                     as.data.frame(SummarizedExperiment::colData(example_whe)))
    expect_identical(transmissionPairs(cons), transmissionPairs(example_whe))
    expect_identical(GenomeInfoDb::seqlengths(cons),
                     GenomeInfoDb::seqlengths(example_whe))

    ## s2 has no consensus site and stays as an empty column; s3 is
    ## known to the second caller only and takes its row from there
    cd_a <- S4Vectors::DataFrame(sample_id = c("s1", "s2"),
                                 host_id = c("h1", "h2"))
    cd_b <- S4Vectors::DataFrame(sample_id = c("s3", "s1"),
                                 host_id = c("h3", "h1-b"),
                                 role = c("recipient", "donor"))
    a <- .caller_whe("seg", c(10, 20), c("T", "G"),
                     cbind(s1 = c(0.10, 0.20), s2 = c(NA, 0.90)),
                     colData = cd_a)
    b <- .caller_whe("seg", c(10, 20), c("T", "G"),
                     cbind(s3 = c(0.50, NA), s1 = c(0.11, 0.21)),
                     colData = cd_b)
    cons <- consensusISNV(list(a = a, b = b))
    expect_identical(colnames(cons), c("s1", "s2", "s3"))
    cd <- SummarizedExperiment::colData(cons)
    expect_identical(cd$sample_id, c("s1", "s2", "s3"))
    expect_identical(cd$host_id, c("h1", "h2", "h3"))
    expect_identical(cd$role, c(NA, NA, "recipient"))
    expect_identical(unname(SummarizedExperiment::assay(cons, "altFreq")),
                     cbind(c(mean(c(0.10, 0.11)), mean(c(0.20, 0.21))),
                           NA_real_, NA_real_))
})

test_that("calcCallerConcordance returns correct pairwise stats from VCFs", {
    whes <- .load_caller_whes()
    result <- calcCallerConcordance(whes)
    expect_s4_class(result, "DataFrame")
    ## 3 callers -> 3 pairs
    expect_equal(nrow(result), 3L)
    expect_true(all(c("caller1", "caller2", "n_shared",
                      "n_only1", "n_only2", "jaccard") %in%
                    colnames(result)))
    ## All callers share all 15 sites -> jaccard = 1
    expect_true(all(result$jaccard == 1.0))
    expect_true(all(result$n_shared == 15L))
    expect_true(all(result$n_only1 == 0L))
    expect_true(all(result$n_only2 == 0L))
})

test_that("consensusISNV rejects unnamed or single-element lists", {
    whes <- .load_caller_whes()
    unnamed <- whes[c("ivar", "lofreq")]
    names(unnamed) <- NULL
    expect_error(consensusISNV(unnamed), "named list")
    expect_error(consensusISNV(whes[1]), "at least 2")
})
