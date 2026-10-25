# tests/testthat/test-selection.R

## One sample of single-base iSNVs on the given sequences.
pnps_whe <- function(pos, ref, alt, seqnames = "seg1",
                     freqs = rep(0.2, length(pos))) {
    gr <- GenomicRanges::GRanges(seqnames, IRanges::IRanges(pos, width = 1))
    S4Vectors::mcols(gr)$ref <- ref
    S4Vectors::mcols(gr)$alt <- alt
    WithinHostExperiment(
        assays = list(altFreq = matrix(freqs, ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
}

pnps_cds <- function(start, end, strand = "+", ID = "cds-A",
                     gene = "geneA", seqnames = "seg1") {
    GenomicRanges::GRanges(seqnames, IRanges::IRanges(start, end),
        strand = strand, ID = ID, gene = gene)
}

## ATG CTG AAA GGG TAA. Nei-Gojobori synonymous sites: ATG 0, CTG 4/3,
## AAA 1/3, GGG 1; the stop codon carries none, so S = 8/3 and
## N = 12 - 8/3 = 28/3.
ref15 <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))

## 30 sense codons, no stop: (ATG CTG AAA GGG TTT CCC) x 5. The six
## codons carry 0 + 4/3 + 1/3 + 1 + 1/3 + 1 = 4 synonymous sites.
seq90 <- paste(rep(c("ATG", "CTG", "AAA", "GGG", "TTT", "CCC"), 5L),
               collapse = "")

test_that("NG86 synonymous sites match hand-computed values", {
    sites <- WithinHostExperiment:::.NG_SYN_SITES
    expect_equal(unname(sites["ATG"]), 0)
    expect_equal(unname(sites["CTG"]), 4 / 3)
    expect_equal(unname(sites["AAA"]), 1 / 3)
    expect_equal(unname(sites["GGG"]), 1)
    expect_equal(unname(sites["TTT"]), 1 / 3)
    expect_equal(unname(sites["CCC"]), 1)
    expect_false(any(c("TAA", "TAG", "TGA") %in% names(sites)))
    expect_length(sites, 61L)
})

test_that("NG86 sites of the 61 sense codons equal those of MSA2dist", {
    skip_if_not_installed("MSA2dist")
    sites <- WithinHostExperiment:::.NG_SYN_SITES
    ## a one-codon sequence does not make a codon matrix, so each codon
    ## is read twice and its sites halved
    msa <- t(vapply(names(sites), function(codon) {
        two <- strrep(codon, 2L)
        cm <- MSA2dist::dnastring2codonmat(
            Biostrings::DNAStringSet(c(a = two, b = two)))
        r <- unlist(MSA2dist::codonmat2pnps(cm))
        as.numeric(r[c("S", "N")]) / 2
    }, numeric(2)))
    expect_equal(unname(msa[, 1L]), unname(sites))
    expect_equal(unname(msa[, 2L]), unname(3 - sites))
})

test_that("pNpSWithinHost normalises by synonymous and nonsynonymous sites", {
    skip_if_not_installed("txdbmaker")
    ## 1 A>G: ATG -> GTG (M -> V); 6 G>A: CTG -> CTA (L -> L)
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("G", "A"))
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_s4_class(res, "DataFrame")
    expect_equal(res$nS[1], 1L)
    expect_equal(res$nN[1], 1L)
    expect_equal(res$S_sites[1], 8 / 3)
    expect_equal(res$N_sites[1], 28 / 3)
    expect_equal(res$pS[1], 3 / 8)
    expect_equal(res$pN[1], 3 / 28)
    ## no multiple-hit correction: the ratio is (3/28) / (3/8)
    expect_equal(res$pNpS[1], 2 / 7)
    expect_identical(res$gene[1], "geneA")
    expect_equal(res$gene_pNpS[1], 2 / 7)
    expect_false(any(grepl("dNdS|Jukes", colnames(res))))
})

test_that("pNpSWithinHost counts a nonsense change as nonsynonymous", {
    skip_if_not_installed("txdbmaker")
    ## 7 A>T: AAA -> TAA (K -> stop)
    whe <- pnps_whe(c(6, 7), c("G", "A"), c("A", "T"))
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(res$nS[1], 1L)
    expect_equal(res$nN[1], 1L)
})

test_that("iSNVs in a stop codon are not counted (was: counted in nS and nN)", {
    skip_if_not_installed("txdbmaker")
    ## 13 T>C is TAA -> CAA (stop lost) and 15 A>G is TAA -> TAG; the
    ## stop codon carries no sites, so neither counts
    whe <- pnps_whe(c(1, 6, 13, 15), c("A", "G", "T", "A"),
                    c("G", "A", "C", "G"))
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(c(res$nS[1], res$nN[1]), c(1L, 1L))
    expect_equal(res$S_sites[1], 8 / 3)
    expect_equal(res$pNpS[1], 2 / 7)
    expect_equal(c(res$gene_nS[1], res$gene_nN[1]), c(1L, 1L))
    expect_equal(res$gene_pNpS[1], 2 / 7)
    ## on the minus strand, plus 15 A>G is minus 1 T>C (TAA -> TAG)
    minus <- Biostrings::reverseComplement(ref15)
    res_minus <- pNpSWithinHost(pnps_whe(c(1, 10), c("T", "C"), c("C", "T")),
                                pnps_cds(1, 15, strand = "-"), minus)
    expect_equal(c(res_minus$nS[1], res_minus$nN[1]), c(1L, 0L))
})

test_that("ATG -> CTG at the first base counts as nonsynonymous (was: synonymous by the initiator rule)", {
    skip_if_not_installed("txdbmaker")
    ## ATG carries no synonymous site, so a change at its first base is
    ## nonsynonymous for pN/pS, although annotateCodonChange() follows
    ## predictCoding() and calls it synonymous
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("C", "A"))
    ann <- suppressMessages(annotateCodonChange(whe, pnps_cds(1, 15), ref15))
    expect_identical(
        SummarizedExperiment::rowRanges(ann)$CONSEQUENCE[1], "synonymous")
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(c(res$nS[1], res$nN[1]), c(1L, 1L))
    expect_equal(c(res$gene_nS[1], res$gene_nN[1]), c(1L, 1L))
    expect_equal(res$pNpS[1], 2 / 7)
})

test_that("an iSNV in one transcript's stop codon counts by a transcript reading it in a sense codon", {
    skip_if_not_installed("txdbmaker")
    ## tx1 (g1): ATG AAA TAA GGG CCC AAA; tx2 (g2), 5-13 in another
    ## frame: AAT AAG GGC. 8 A>G is TAA -> TGA in tx1, the longer, and
    ## AAG -> GAG (K -> E) in tx2
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGAAATAAGGGCCCAAA"))
    cds <- pnps_cds(c(1, 5), c(18, 13), ID = c("tx1", "tx2"),
                    gene = c("g1", "g2"))
    whe <- pnps_whe(c(8, 12), c("A", "G"), c("G", "A"))
    expect_warning(res <- pNpSWithinHost(whe, cds, ref), "0 codon")
    res <- as.data.frame(res)
    ## 12 G>A is GGG -> GGA in tx1 (synonymous) and GGC -> GAC in tx2
    expect_identical(res$nS[1], 1L)
    expect_identical(res$nN[1], 1L)
    expect_identical(res$gene, c("g1", "g2"))
    expect_identical(res$gene_nS, c(1L, 0L))
    expect_identical(res$gene_nN, c(0L, 2L))
})

test_that("iSNVs outside every coding sequence are not counted", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA",
                                      seg2 = "CCCCCCCCCC"))
    whe <- pnps_whe(c(1, 6, 5), c("A", "G", "C"), c("G", "A", "T"),
                    seqnames = c("seg1", "seg1", "seg2"))
    expect_warning(res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref),
                   "seg2")
    expect_equal(res$nN[1], 1L)
    expect_equal(res$nS[1], 1L)
})

test_that("pNpSWithinHost returns NA when no synonymous iSNV is present", {
    skip_if_not_installed("txdbmaker")
    res <- pNpSWithinHost(pnps_whe(1, "A", "G"), pnps_cds(1, 15), ref15)
    expect_equal(res$nS[1], 0L)
    expect_equal(res$nN[1], 1L)
    expect_true(is.na(res$pNpS[1]))
    expect_true(is.na(res$gene_pNpS[1]))
})

test_that("pNpSWithinHost respects qcPass and usePassedOnly", {
    skip_if_not_installed("txdbmaker")
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("G", "A"))
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <-
        matrix(c(FALSE, TRUE), nrow = 2, ncol = 1)
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(c(res$nS[1], res$nN[1]), c(1L, 0L))
    all_calls <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15,
                                usePassedOnly = FALSE)
    expect_equal(c(all_calls$nS[1], all_calls$nN[1]), c(1L, 1L))
    expect_error(pNpSWithinHost(whe, pnps_cds(1, 15), ref15,
                                usePassedOnly = NA), "TRUE or FALSE")
})

test_that("an iSNV whose qcPass is NA is not counted (was: counted as passing)", {
    skip_if_not_installed("txdbmaker")
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("G", "A"))
    SummarizedExperiment::assay(whe, "qcPass", withDimnames = FALSE) <-
        matrix(c(NA, TRUE), nrow = 2, ncol = 1)
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(c(res$nS[1], res$nN[1]), c(1L, 0L))
})

test_that("an iSNV at zero frequency is not present", {
    skip_if_not_installed("txdbmaker")
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("G", "A"), freqs = c(0, 0.2))
    res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    expect_equal(c(res$nS[1], res$nN[1]), c(1L, 0L))
})

test_that("the site count equals three times the distinct sense codons", {
    skip_if_not_installed("txdbmaker")
    ## tx1 reads codons 1-20 and tx2 codons 11-30 in the same frame; the
    ## ten codons they share are counted once in the total and once in
    ## each gene.
    ref <- Biostrings::DNAStringSet(c(seg1 = seq90))
    cds <- pnps_cds(c(1, 31), c(60, 90), ID = c("tx1", "tx2"),
                    gene = c("g1", "g2"))
    ## 10 G>A is the first base of codon 4 (GGG), in tx1 only; 70 C>A is
    ## the first base of codon 24 (CCC), in tx2 only
    whe <- pnps_whe(c(10, 70), c("G", "C"), c("A", "A"))
    expect_warning(res <- pNpSWithinHost(whe, cds, ref),
                   "tx1 \\(g1\\) / tx2 \\(g2\\): 30 nt, 10 codon\\(s\\) shared")
    expect_equal(res$S_sites[1] + res$N_sites[1], 3 * 30)
    expect_equal(res$S_sites[1], 20)
    expect_equal(res$N_sites[1], 70)
    by_gene <- as.data.frame(res)[, c("gene", "gene_S_sites",
                                      "gene_N_sites")]
    expect_identical(by_gene$gene, c("g1", "g2"))
    expect_equal(by_gene$gene_S_sites + by_gene$gene_N_sites, c(60, 60))
    expect_equal(by_gene$gene_S_sites, c(40 / 3, 40 / 3))
    ## the genes' sites add to more than the total: the shared codons
    ## belong to both
    expect_gt(sum(by_gene$gene_S_sites + by_gene$gene_N_sites),
              res$S_sites[1] + res$N_sites[1])
})

test_that("an iSNV in a codon two genes share counts in both genes (was: only in the longest transcript's)", {
    skip_if_not_installed("txdbmaker")
    ## as above, tx1 (g1) reads codons 1-20 and tx2 (g2) codons 11-30.
    ## 10 G>A: GGG -> AGG in codon 4, tx1 only, nonsynonymous;
    ## 42 G>A: CTG -> CTA in codon 14, shared, synonymous;
    ## 75 G>A: ATG -> ATA in codon 25, tx2 only, nonsynonymous
    ref <- Biostrings::DNAStringSet(c(seg1 = seq90))
    cds <- pnps_cds(c(1, 31), c(60, 90), ID = c("tx1", "tx2"),
                    gene = c("g1", "g2"))
    whe <- pnps_whe(c(10, 42, 75), c("G", "G", "G"), c("A", "A", "A"))
    expect_warning(res <- pNpSWithinHost(whe, cds, ref), "10 codon")
    res <- as.data.frame(res)
    expect_identical(res$gene, c("g1", "g2"))
    expect_identical(res$gene_nS, c(1L, 1L))
    expect_identical(res$gene_nN, c(1L, 1L))
    ## and once in the total
    expect_identical(res$nS[1], 1L)
    expect_identical(res$nN[1], 2L)
    ## each gene's ratio rests on its own sites and its own iSNVs
    expect_equal(res$gene_pNpS,
                 (1 / res$gene_N_sites) / (1 / res$gene_S_sites))
})

test_that("coding sequences in different frames add their codons", {
    skip_if_not_installed("txdbmaker")
    ## tx1: ATG GAT GCT GCT GCT; tx2, one frame along: ATG CTG CTG
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGGATGCTGCTGCT"))
    cds <- pnps_cds(c(1, 5), c(15, 13), ID = c("tx1", "tx2"),
                    gene = c("g1", "g2"))
    whe <- pnps_whe(14, "C", "A")
    expect_warning(res <- pNpSWithinHost(whe, cds, ref),
                   "9 nt, 0 codon\\(s\\) shared")
    ## 5 + 3 distinct sense codons
    expect_equal(res$S_sites[1] + res$N_sites[1], 3 * 8)
})

test_that("a frameshifted transcript is counted as one coding sequence", {
    skip_if_not_installed("txdbmaker")
    ## parts 1-30 and 30-59 share nucleotide 30, as ORF1ab does at its
    ## -1 frameshift; the joined transcript is 60 nt of 20 sense codons
    ## and overlaps no other transcript
    ref <- Biostrings::DNAStringSet(c(seg1 = seq90))
    cds <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(1, 30), c(30, 59)), strand = "+",
        ID = "tx1", gene = "g1")
    whe <- pnps_whe(10, "G", "A")
    expect_no_warning(res <- pNpSWithinHost(whe, cds, ref))
    expect_equal(res$S_sites[1] + res$N_sites[1], 3 * 20)
})

test_that("minus-strand sites match the plus strand of the complement", {
    skip_if_not_installed("txdbmaker")
    minus <- Biostrings::reverseComplement(ref15)
    ## plus 1 A>G (ATG -> GTG) is minus 15 T>C
    plus_res <- pNpSWithinHost(pnps_whe(1, "A", "G"), pnps_cds(1, 15),
                               ref15)
    minus_res <- pNpSWithinHost(pnps_whe(15, "T", "C"),
                                pnps_cds(1, 15, strand = "-"), minus)
    expect_equal(minus_res$S_sites[1], plus_res$S_sites[1])
    expect_equal(minus_res$N_sites[1], plus_res$N_sites[1])
    expect_equal(minus_res$nN[1], 1L)
})

test_that("a codon with an ambiguous base carries no sites and no count", {
    skip_if_not_installed("txdbmaker")
    ## ATG CNG AAA GGG TAA: S = 0 + 1/3 + 1 = 4/3 over three codons
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCNGAAAGGGTAA"))
    whe <- pnps_whe(c(6, 7), c("G", "A"), c("A", "T"))
    expect_warning(res <- pNpSWithinHost(whe, pnps_cds(1, 15), ref),
                   "ambiguous base")
    expect_equal(res$S_sites[1], 4 / 3)
    expect_equal(res$N_sites[1], 9 - 4 / 3)
    ## only the nonsense change at 7 is counted
    expect_equal(c(res$nS[1], res$nN[1]), c(0L, 1L))
})

test_that("pNpSWithinHost takes the same inputs as annotateCodonChange", {
    skip_if_not_installed("txdbmaker")
    whe <- pnps_whe(c(1, 6), c("A", "G"), c("G", "A"))
    from_gr <- pNpSWithinHost(whe, pnps_cds(1, 15), ref15)
    grl <- GenomicRanges::GRangesList(cds_A = GenomicRanges::GRanges("seg1",
        IRanges::IRanges(1, 15), strand = "+"))
    S4Vectors::mcols(grl)$gene <- "geneA"
    from_grl <- pNpSWithinHost(whe, grl, ref15, geneCol = "gene")
    expect_equal(as.data.frame(from_grl), as.data.frame(from_gr))
    txdb <- suppressWarnings(txdbmaker::makeTxDb(
        data.frame(tx_id = 1L, tx_name = "cds-A", tx_chrom = "seg1",
                   tx_strand = "+", tx_start = 1L, tx_end = 15L),
        data.frame(tx_id = 1L, exon_rank = 1L, exon_start = 1L,
                   exon_end = 15L, cds_start = 1L, cds_end = 15L),
        data.frame(tx_id = 1L, gene_id = "geneA"),
        chrominfo = data.frame(chrom = "seg1", length = 15L,
                               is_circular = FALSE)))
    from_txdb <- pNpSWithinHost(whe, txdb, ref15)
    expect_equal(as.data.frame(from_txdb), as.data.frame(from_gr))
})

test_that("pNpSWithinHost rejects file paths and unused arguments", {
    whe <- pnps_whe(1, "A", "G")
    expect_error(pNpSWithinHost(whe, pnps_cds(1, 15), "ref.fa"),
                 "readDNAStringSet")
    expect_error(pNpSWithinHost(whe, "genes.gff3", ref15),
                 "rtracklayer::import")
    expect_error(pNpSWithinHost(whe, pnps_cds(1, 15), ref15[[1L]]),
                 "DNAStringSet\\(c\\(chr1 = ref\\)\\)")
    skip_if_not_installed("txdbmaker")
    expect_error(pNpSWithinHost(whe, pnps_cds(1, 15), ref15,
                                jukesCantor = TRUE),
                 "Unused argument")
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
