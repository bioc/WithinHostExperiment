# tests/testthat/test-annotate-genes.R

genes_whe <- function(pos, seqnames = "chr1") {
    gr <- GenomicRanges::GRanges(seqnames,
        IRanges::IRanges(pos, width = 1))
    S4Vectors::mcols(gr)$ref <- rep("A", length(pos))
    S4Vectors::mcols(gr)$alt <- rep("G", length(pos))
    WithinHostExperiment(
        assays = list(altFreq = matrix(0.1, nrow = length(pos), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
}

gene_col <- function(whe) {
    S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$GENE
}

test_that("annotateGenes labels each site with the gene it falls in", {
    whe <- genes_whe(c(5, 50, 150, 300))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 101), c(100, 200)),
        strand = c("+", "-"), Name = c("gA", "gB"))
    out <- annotateGenes(whe, genes)
    expect_s4_class(out, "WithinHostExperiment")
    ## the minus-strand gene still labels its unstranded site
    expect_identical(gene_col(out), c("gA", "gA", "gB", NA))
    rr <- SummarizedExperiment::rowRanges(out)
    expect_identical(S4Vectors::mcols(rr)$ref, rep("A", 4L))
    expect_identical(qcLog(out), qcLog(whe))
})

test_that("annotateGenes stops when sites and genes share no sequence name", {
    whe <- genes_whe(c(266, 21563), "MN908947.3")
    genes <- GenomicRanges::GRanges("NC_045512.2",
        IRanges::IRanges(c(266, 21563), c(21555, 25384)),
        Name = c("ORF1ab", "S"))
    expect_error(annotateGenes(whe, genes), "MN908947\\.3")
    expect_error(annotateGenes(whe, genes), "NC_045512\\.2")
    expect_error(annotateGenes(whe, genes), "share no sequence names")
})

test_that("annotateGenes warns about variant sequences with no gene range", {
    whe <- genes_whe(c(10, 10), c("seg1", "seg2"))
    genes <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(1, 100), Name = "g1")
    expect_warning(out <- annotateGenes(whe, genes), "seg2")
    expect_identical(gene_col(out), c("g1", NA))
})

test_that("annotateGenes names the available columns when nameCol is wrong", {
    whe <- genes_whe(5)
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100),
        Name = "gA", ID = "gene-A")
    expect_error(annotateGenes(whe, genes, nameCol = "gene_name"),
                 "Available columns: Name, ID")
    out <- annotateGenes(whe, genes, nameCol = "ID")
    expect_identical(gene_col(out), "gene-A")
})

test_that("annotateGenes rejects a label column with several values per range", {
    whe <- genes_whe(5)
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
    S4Vectors::mcols(genes)$Parent <- IRanges::CharacterList("a")
    expect_error(annotateGenes(whe, genes, nameCol = "Parent"),
                 "one label per range")
})

test_that("annotateGenes keeps the longest of overlapping genes by default", {
    whe <- genes_whe(c(95, 110, 50))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 90), c(100, 120)),
        Name = c("gLong", "gShort"))
    expect_message(out <- annotateGenes(whe, genes),
                   "1 site\\(s\\) fall in more than one gene")
    expect_identical(gene_col(out), c("gLong", "gShort", "gLong"))
})

test_that("annotateGenes breaks a tie in length by the order of 'genes'", {
    whe <- genes_whe(15)
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(11, 1), width = 20), Name = c("second", "first"))
    out <- suppressMessages(annotateGenes(whe, genes))
    expect_identical(gene_col(out), "second")
})

test_that("annotateGenes keeps every overlapping label with collapse = 'all'", {
    whe <- genes_whe(c(95, 110, 50, 500))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 90), c(100, 120)),
        Name = c("gLong", "gShort"))
    out <- annotateGenes(whe, genes, collapse = "all")
    g <- gene_col(out)
    expect_s4_class(g, "CharacterList")
    expect_identical(as.list(g),
        list(c("gLong", "gShort"), "gShort", "gLong", character()))
})

test_that("annotateGenes ignores unlabelled ranges with a warning", {
    whe <- genes_whe(c(5, 150))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 101), c(100, 200)), Name = c("gA", NA))
    expect_warning(out <- annotateGenes(whe, genes),
                   "1 range\\(s\\) in 'genes' have no 'Name'")
    expect_identical(gene_col(out), c("gA", NA))
})

test_that("annotateGenes rejects empty genes and unused arguments", {
    whe <- genes_whe(5)
    expect_error(annotateGenes(whe, GenomicRanges::GRanges()),
                 "contains no ranges")
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100),
        Name = "gA")
    expect_error(annotateGenes(whe, genes, featureType = "gene"),
                 "Unused argument")
})

test_that("annotateGenes counts one gene split into several ranges once", {
    ## two ranges of one gene, as for a CDS split by a frameshift
    whe <- genes_whe(c(50, 150))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 40), c(100, 60)), Name = c("ORF1ab", "ORF1ab"))
    expect_no_message(out <- annotateGenes(whe, genes))
    expect_identical(gene_col(out), c("ORF1ab", NA))
    out_all <- annotateGenes(whe, genes, collapse = "all")
    expect_identical(as.list(gene_col(out_all)), list("ORF1ab", character()))
})

test_that("annotateGenes reports every site that carries several genes", {
    whe <- genes_whe(c(50, 55, 150))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 40, 45), c(100, 60, 70)),
        Name = c("gA", "gB", "gC"))
    expect_message(annotateGenes(whe, genes),
                   "^2 site\\(s\\) fall in more than one gene")
})

test_that("annotateGenes stops when no range has a label", {
    whe <- genes_whe(5)
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100),
        Name = NA_character_)
    expect_error(annotateGenes(whe, genes),
                 "No range in 'genes' has a 'Name' label")
})

test_that("annotateGenes reads labels stored as Rle or factor", {
    whe <- genes_whe(c(5, 150))
    genes <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 101), c(100, 200)))
    S4Vectors::mcols(genes)$Name <- S4Vectors::Rle(c("gA", "gB"))
    expect_identical(gene_col(annotateGenes(whe, genes)), c("gA", "gB"))
    S4Vectors::mcols(genes)$Name <- factor(c("zeta", "alpha"))
    expect_identical(gene_col(annotateGenes(whe, genes)),
                     c("zeta", "alpha"))
})

test_that("annotateGenes validates nameCol", {
    whe <- genes_whe(5)
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100),
        Name = "gA", ID = "gene-A")
    expect_error(annotateGenes(whe, genes, nameCol = NA_character_),
                 "single column name")
    expect_error(annotateGenes(whe, genes, nameCol = c("Name", "ID")),
                 "single column name")
})

test_that("annotateGenes handles an object with no sites", {
    whe <- genes_whe(5)[0, ]
    genes <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100),
        Name = "gA")
    expect_identical(gene_col(annotateGenes(whe, genes)), character())
    expect_length(gene_col(annotateGenes(whe, genes, collapse = "all")), 0L)
})

test_that("annotateGenes takes genes from a TxDb", {
    skip_if_not_installed("txdbmaker")
    transcripts <- data.frame(tx_id = 1:2, tx_name = c("t1", "t2"),
        tx_chrom = "chr1", tx_strand = "+",
        tx_start = c(1L, 101L), tx_end = c(100L, 200L))
    splicings <- data.frame(tx_id = 1:2, exon_rank = 1L,
        exon_start = c(1L, 101L), exon_end = c(100L, 200L))
    genes <- data.frame(tx_id = 1:2, gene_id = c("gA", "gB"))
    chrominfo <- data.frame(chrom = "chr1", length = 300L,
        is_circular = FALSE)
    ## makeTxDb warns that no genome version was given; that is expected
    txdb <- suppressWarnings(txdbmaker::makeTxDb(transcripts, splicings,
        genes, chrominfo = chrominfo))
    out <- annotateGenes(genes_whe(c(5, 150, 250)), txdb)
    expect_identical(gene_col(out), c("gA", "gB", NA))
})
