## Helper: a 15-base CDS (ATG CTG AAA GGG TAA) with a GFF3 and FASTA.
## Nei-Gojobori synonymous sites: ATG 0, CTG 4/3, AAA 1/3, GGG 1; the stop
## codon is skipped, so S = 8/3 and N = 12 - 8/3 = 28/3.
.ng_fixture <- function(positions, refs, alts,
                        freqs = rep(0.2, length(positions))) {
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "ATGCTGAAAGGGTAA", ">seg2", "CCCCCCCCCC"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
    gr <- GenomicRanges::GRanges(names(positions),
        IRanges::IRanges(unname(positions), width = 1))
    S4Vectors::mcols(gr)$ref <- refs
    S4Vectors::mcols(gr)$alt <- alts
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(freqs, ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    skip_if_not_installed("txdbmaker")
    cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
        strand = "+", ID = "cds-A", gene = "geneA")
    ## seg2 carries no CDS; annotateCodonChange() says so with a warning
    whe <- withCallingHandlers(
        annotateCodonChange(whe, cds, Biostrings::readDNAStringSet(ref)),
        warning = function(w) {
            if (grepl("no ranges on variant seqlevel", conditionMessage(w)))
                invokeRestart("muffleWarning")
        })
    list(whe = whe, gff = gff, ref = ref)
}
.jc <- function(p) -0.75 * log(1 - 4 * p / 3)

test_that("NG86 synonymous sites match hand-computed values", {
    sites <- WithinHostExperiment:::.NG_SYN_SITES
    expect_equal(unname(sites["ATG"]), 0)
    expect_equal(unname(sites["CTG"]), 4 / 3)
    expect_equal(unname(sites["AAA"]), 1 / 3)
    expect_equal(unname(sites["GGG"]), 1)
    expect_false("TAA" %in% names(sites))
    expect_length(sites, 61L)
})

test_that("dndsWithinHost normalises by synonymous and nonsynonymous sites", {
    ## position 1 A>G: ATG -> GTG (M -> V, nonsynonymous)
    ## position 6 G>A: CTG -> CTA (L -> L, synonymous)
    fx <- .ng_fixture(c(seg1 = 1, seg1 = 6), c("A", "G"), c("G", "A"))
    res <- dndsWithinHost(fx$whe, gff = fx$gff, refFasta = fx$ref)
    expect_s4_class(res, "DataFrame")
    expect_equal(res$nS[1], 1L)
    expect_equal(res$nN[1], 1L)
    expect_equal(res$S_sites[1], 8 / 3)
    expect_equal(res$N_sites[1], 28 / 3)
    expect_equal(res$pS[1], 1 / (8 / 3))
    expect_equal(res$pN[1], 1 / (28 / 3))
    expect_equal(res$dNdS[1], .jc(3 / 28) / .jc(3 / 8), tolerance = 1e-12)
    expect_equal(res$gene[1], "geneA")
    expect_equal(res$gene_dNdS[1], res$dNdS[1])
})

test_that("variants outside CDS features are not counted as nonsynonymous", {
    fx <- .ng_fixture(c(seg1 = 1, seg1 = 6, seg2 = 5), c("A", "G", "C"),
                      c("G", "A", "T"))
    res <- dndsWithinHost(fx$whe, gff = fx$gff, refFasta = fx$ref)
    expect_equal(res$nN[1], 1L)
    expect_equal(res$nS[1], 1L)
})

test_that("dndsWithinHost returns NA when no synonymous iSNV is present", {
    fx <- .ng_fixture(c(seg1 = 1), "A", "G")
    res <- dndsWithinHost(fx$whe, gff = fx$gff, refFasta = fx$ref)
    expect_equal(res$nS[1], 0L)
    expect_true(is.na(res$dNdS[1]))
})

test_that("dndsWithinHost respects qcPass", {
    fx <- .ng_fixture(c(seg1 = 1, seg1 = 6), c("A", "G"), c("G", "A"))
    whe <- fx$whe
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <-
        matrix(FALSE, nrow = 2, ncol = 1)
    res <- dndsWithinHost(whe, gff = fx$gff, refFasta = fx$ref)
    expect_equal(res$nS[1] + res$nN[1], 0L)
})

test_that("dndsWithinHost requires the GFF3, FASTA and annotation", {
    fx <- .ng_fixture(c(seg1 = 1), "A", "G")
    expect_error(dndsWithinHost(fx$whe), "required")
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(100, width = 1))
    S4Vectors::mcols(gr)$ref <- "A"
    S4Vectors::mcols(gr)$alt <- "T"
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(0.1, ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    expect_error(dndsWithinHost(whe, gff = fx$gff, refFasta = fx$ref),
                 "REFAA")
})

test_that("codons shared by overlapping CDS features are counted once", {
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "ATGCTGAAAGGGTAA"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA",
        "seg1\t.\tCDS\t1\t9\t.\t+\t0\tgene=geneA"), gff)
    sites <- WithinHostExperiment:::.cds_site_counts(
        WithinHostExperiment:::.read_gff_cds(gff),
        WithinHostExperiment:::.read_fasta_simple(ref))
    expect_equal(sites$S_sites, 8 / 3)
    expect_equal(sites$N_sites, 28 / 3)
})

test_that("estimateSelectionCoefficient basic computation", {
    ## logit(0.5) = 0, logit(~0.731) = 1
    p1 <- 0.5
    p2 <- exp(1) / (1 + exp(1))  # logit(p2) = 1
    s <- estimateSelectionCoefficient(c(p1, p2))
    expect_equal(s, 1.0, tolerance = 1e-10)
})

test_that("estimateSelectionCoefficient with custom times", {
    freqs <- c(0.1, 0.3, 0.5)
    s <- estimateSelectionCoefficient(freqs, times = c(0, 5, 10))
    expect_length(s, 2L)
    expect_true(all(is.finite(s)))
})

test_that("estimateSelectionCoefficient errors on invalid input", {
    expect_error(estimateSelectionCoefficient(0.5), "at least 2")
    expect_error(estimateSelectionCoefficient(c(0, 0.5)),
                 "strictly between 0 and 1")
    expect_error(estimateSelectionCoefficient(c(0.5, NA)),
                 "NA")
})

test_that("estimateSelectionCoefficient respects generation_time", {
    freqs <- c(0.2, 0.4)
    s1 <- estimateSelectionCoefficient(freqs, generation_time = 1)
    s2 <- estimateSelectionCoefficient(freqs, generation_time = 2)
    expect_equal(s2, s1 * 2, tolerance = 1e-10)
})

test_that("estimateSelectionCoefficient rejects decreasing times", {
    expect_error(estimateSelectionCoefficient(c(0.1, 0.2, 0.3),
                                              times = c(0, 5, 3)),
                 "strictly increasing")
})

test_that("the FASTA reader ignores whitespace and empty final records", {
    ## .read_fasta_simple() still feeds the dN/dS site counts; stripping
    ## whitespace inside a record keeps every later coordinate in place.
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">seg1 description", "ATG CTG ", "AAA\t", "", ">empty"), fa)
    seqs <- WithinHostExperiment:::.read_fasta_simple(fa)
    expect_identical(seqs[["seg1"]], "ATGCTGAAA")
    expect_identical(seqs[["empty"]], "")
})

test_that("dndsWithinHost needs one transcript per site", {
    skip_if_not_installed("txdbmaker")
    ## collapse = "all" stores a CharacterList, which cannot be read as
    ## one amino acid per site
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "ATGCTGAAAGGGTAA"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(c(1, 6), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "G")
    S4Vectors::mcols(gr)$alt <- c("G", "A")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.2, 0.2), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
        strand = "+", ID = "cds-A", gene = "geneA")
    all_tx <- annotateCodonChange(whe, cds,
        Biostrings::readDNAStringSet(ref), collapse = "all")
    expect_error(dndsWithinHost(all_tx, gff = gff, refFasta = ref),
                 'collapse = "longest"')
})

test_that("dndsWithinHost warns when the gene labels do not meet", {
    skip_if_not_installed("txdbmaker")
    ## geneCol = NULL reports transcript names in GENEID, which the GFF3
    ## gene attribute does not contain; every count would otherwise be
    ## zero with no explanation
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "ATGCTGAAAGGGTAA"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(c(1, 6), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "G")
    S4Vectors::mcols(gr)$alt <- c("G", "A")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.2, 0.2), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
        strand = "+", ID = "cds-A")
    ann <- annotateCodonChange(whe, cds,
        Biostrings::readDNAStringSet(ref), geneCol = NULL)
    expect_warning(res <- dndsWithinHost(ann, gff = gff, refFasta = ref),
                   "cds-A")
    expect_equal(res$nS[1], 0L)
    expect_equal(res$nN[1], 0L)
})

test_that("a codon with an ambiguous base is counted in neither nS nor nN", {
    skip_if_not_installed("txdbmaker")
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "ATGCNGAAAGGGTAA"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(c(6, 7), width = 1))
    S4Vectors::mcols(gr)$ref <- c("G", "A")
    S4Vectors::mcols(gr)$alt <- c("A", "T")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.2, 0.2), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
        strand = "+", ID = "cds-A", gene = "geneA")
    ann <- suppressWarnings(annotateCodonChange(whe, cds,
        Biostrings::readDNAStringSet(ref)))
    res <- dndsWithinHost(ann, gff = gff, refFasta = ref)
    ## only the nonsense change at position 7 is counted
    expect_equal(res$nS[1], 0L)
    expect_equal(res$nN[1], 1L)
})
