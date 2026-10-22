# tests/testthat/test-consensus-sequences.R

cons_ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))

## Sites on seg1 with one frequency column per host.
cons_whe <- function(pos, ref, alt, freqs, seqnames = "seg1", qc = NULL) {
    freqs <- as.matrix(freqs)
    if (is.null(colnames(freqs)))
        colnames(freqs) <- paste0("host", seq_len(ncol(freqs)))
    gr <- GenomicRanges::GRanges(seqnames,
        IRanges::IRanges(pos, width = nchar(ref)))
    S4Vectors::mcols(gr)$ref <- ref
    S4Vectors::mcols(gr)$alt <- alt
    assays <- list(altFreq = freqs)
    if (!is.null(qc)) assays$qcPass <- qc
    WithinHostExperiment(
        assays = assays, rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = colnames(freqs),
                                       row.names = colnames(freqs)))
}

test_that("consensusSequences writes in each host's majority substitutions", {
    whe <- cons_whe(c(1, 6, 9), c("A", "G", "A"), c("G", "A", "G"),
                    cbind(host1 = c(0.95, 0.10, 0.02),
                          host2 = c(0.03, 0.88, 0.04),
                          host3 = c(0.01, 0.02, 0.97)))
    cons <- consensusSequences(whe, cons_ref)
    expect_s4_class(cons, "DNAStringSet")
    expect_identical(names(cons), c("host1", "host2", "host3"))
    expect_identical(as.character(cons),
                     c(host1 = "GTGCTGAAAGGGTAA",
                       host2 = "ATGCTAAAAGGGTAA",
                       host3 = "ATGCTGAAGGGGTAA"))
    expect_identical(S4Vectors::mcols(cons)$n_substitutions, c(1L, 1L, 1L))
    expect_true(all(Biostrings::width(cons) == 15L))
})

test_that("minFreq is a lower bound that includes its own value", {
    whe <- cons_whe(c(1, 6), c("A", "G"), c("G", "A"),
                    cbind(host1 = c(0.5, 0.49)))
    expect_identical(as.character(consensusSequences(whe, cons_ref)),
                     c(host1 = "GTGCTGAAAGGGTAA"))
    expect_identical(as.character(consensusSequences(whe, cons_ref,
                                                     minFreq = 0.4)),
                     c(host1 = "GTGCTAAAAGGGTAA"))
    expect_error(consensusSequences(whe, cons_ref, minFreq = 0),
                 "single number in \\(0, 1\\]")
    expect_error(consensusSequences(whe, cons_ref, minFreq = 1.5),
                 "single number in \\(0, 1\\]")
})

test_that("QC-failed and missing calls are not written in", {
    whe <- cons_whe(c(1, 6, 9), c("A", "G", "A"), c("G", "A", "G"),
                    cbind(host1 = c(0.9, 0.9, NA)),
                    qc = matrix(c(FALSE, TRUE, TRUE), ncol = 1))
    expect_identical(as.character(consensusSequences(whe, cons_ref)),
                     c(host1 = "ATGCTAAAAGGGTAA"))
    expect_identical(as.character(consensusSequences(whe, cons_ref,
                                                     usePassedOnly = FALSE)),
                     c(host1 = "GTGCTAAAAGGGTAA"))
})

test_that("multi-nucleotide substitutions are written in, indels are not", {
    whe <- cons_whe(c(4, 7, 10), c("CT", "A", "GGG"), c("GA", "AT", "G"),
                    cbind(host1 = c(0.9, 0.9, 0.9)))
    expect_warning(cons <- consensusSequences(whe, cons_ref),
                   "2 call\\(s\\) at or above minFreq are insertions")
    expect_identical(as.character(cons), c(host1 = "ATGGAGAAAGGGTAA"))
    expect_identical(S4Vectors::mcols(cons)$n_substitutions, 1L)
})

test_that("a site with the wrong ref allele is not written in", {
    whe <- cons_whe(c(1, 6), c("C", "G"), c("G", "A"),
                    cbind(host1 = c(0.9, 0.9)))
    expect_warning(cons <- consensusSequences(whe, cons_ref),
                   "seg1:1, C against A")
    expect_identical(as.character(cons), c(host1 = "ATGCTAAAAGGGTAA"))
})

test_that("of two overlapping qualifying alleles the more frequent wins", {
    ## two alternative alleles at position 6, both above minFreq = 0.3
    whe <- cons_whe(c(6, 6), c("G", "G"), c("A", "C"),
                    cbind(host1 = c(0.35, 0.60)))
    expect_warning(cons <- consensusSequences(whe, cons_ref, minFreq = 0.3),
                   "1 substitution\\(s\\) overlapped a more frequent one")
    expect_identical(as.character(cons), c(host1 = "ATGCTCAAAGGGTAA"))
})

test_that("region returns the joined coding sequence in its own direction", {
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGCCCAAAGGGTAA"))
    ## coding parts 1-6 and 10-18 skip CCC at 7-9
    whe <- cons_whe(c(1, 8, 12), c("A", "C", "A"), c("G", "T", "G"),
                    cbind(host1 = c(0.9, 0.9, 0.9)))
    plus <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(10, 1), c(18, 6)), strand = "+")
    cons <- consensusSequences(whe, ref, region = plus)
    expect_identical(as.character(cons), c(host1 = "GTGCTGAAGGGGTAA"))
    ## the site at 8 lies outside the region
    expect_identical(S4Vectors::mcols(cons)$n_substitutions, 2L)
    minus <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(1, 10), c(6, 18)), strand = "-")
    rc <- consensusSequences(whe, ref, region = minus)
    expect_identical(as.character(rc), c(host1 = as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(
            "GTGCTGAAGGGGTAA")))))
})

test_that("region must lie on one sequence and strand, inside it", {
    whe <- cons_whe(1, "A", "G", cbind(host1 = 0.9))
    expect_error(consensusSequences(whe, cons_ref,
        region = GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 20))),
        "beyond the 15 bases")
    expect_error(consensusSequences(whe, cons_ref,
        region = GenomicRanges::GRanges("seg1",
            IRanges::IRanges(c(1, 7), c(6, 12)), strand = c("+", "-"))),
        "one strand")
    expect_error(consensusSequences(whe, cons_ref, seqname = "seg1",
        region = GenomicRanges::GRanges("seg2", IRanges::IRanges(1, 6))),
        "lies on seg2")
})

test_that("a multi-sequence reference needs a seqname", {
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA",
                                      seg2 = "CCCCCCCCCC"))
    whe <- cons_whe(c(1, 5), c("A", "C"), c("G", "T"),
                    cbind(host1 = c(0.9, 0.9)), seqnames = c("seg1", "seg2"))
    expect_error(consensusSequences(whe, ref), "seqname.*seg1, seg2")
    expect_identical(as.character(consensusSequences(whe, ref,
                                                     seqname = "seg2")),
                     c(host1 = "CCCCTCCCCC"))
    expect_error(consensusSequences(whe, ref, seqname = "seg3"),
                 "must name one sequence")
})

test_that("consensusSequences stops when sites and ref share no name", {
    whe <- cons_whe(1, "A", "G", cbind(host1 = 0.9),
                    seqnames = "MN908947.3")
    ref <- Biostrings::DNAStringSet(c(NC_045512.2 = "ATGCTGAAAGGGTAA"))
    expect_error(consensusSequences(whe, ref), "MN908947\\.3")
})

test_that("a host with nothing above minFreq returns the reference", {
    whe <- cons_whe(1, "A", "G", cbind(host1 = 0.2, host2 = 0.9))
    cons <- consensusSequences(whe, cons_ref)
    expect_identical(as.character(cons[["host1"]]), "ATGCTGAAAGGGTAA")
    expect_identical(S4Vectors::mcols(cons)$n_substitutions, c(0L, 1L))
})

test_that("consensusSequences rejects paths, single sequences and extras", {
    whe <- cons_whe(1, "A", "G", cbind(host1 = 0.9))
    expect_error(consensusSequences(whe, "ref.fa"), "readDNAStringSet")
    expect_error(consensusSequences(whe, cons_ref[[1L]]),
                 "DNAStringSet\\(c\\(chr1 = ref\\)\\)")
    expect_error(consensusSequences(whe, cons_ref, indels = TRUE),
                 "Unused argument")
})

test_that("consensus sequences can be passed to MSA2dist", {
    skip_if_not_installed("MSA2dist")
    whe <- cons_whe(c(1, 6, 9), c("A", "G", "A"), c("G", "A", "G"),
                    cbind(host1 = c(0.95, 0.10, 0.02),
                          host2 = c(0.03, 0.88, 0.04),
                          host3 = c(0.01, 0.02, 0.97)))
    cons <- consensusSequences(whe, cons_ref)
    kaks <- suppressMessages(MSA2dist::dnastring2kaks(cons, model = "Li"))
    expect_identical(nrow(kaks), 3L)
    expect_setequal(c(kaks$seq1, kaks$seq2), c("host1", "host2", "host3"))
})
