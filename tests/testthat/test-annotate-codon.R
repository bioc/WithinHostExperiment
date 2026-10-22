# tests/testthat/test-annotate-codon.R

## A sequence of 20 sense codons ending in TAA, so every frame-1 codon
## but the last is a sense codon.
codon_seq20 <- paste0("ATGCTGAAAGGGTTTCCCAAACTGATGGTTCTGAAAGGG",
                      "TTTCCCAAACTGATGGTTTAA")

codon_whe <- function(pos, ref, alt, seqnames = "seg1") {
    gr <- GenomicRanges::GRanges(seqnames,
        IRanges::IRanges(pos, width = nchar(ref)))
    S4Vectors::mcols(gr)$ref <- ref
    S4Vectors::mcols(gr)$alt <- alt
    WithinHostExperiment(
        assays = list(altFreq = matrix(0.1, nrow = length(pos), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
}

codon_mc <- function(whe) {
    S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))
}

cds_gr <- function(start, end, strand = "+", ID = "cds-A", gene = "geneA",
                   seqnames = "seg1", ...) {
    GenomicRanges::GRanges(seqnames, IRanges::IRanges(start, end),
        strand = strand, ID = ID, gene = gene, ...)
}

test_that("annotateCodonChange reports codon, amino acid and consequence", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    ## 1 A>G: ATG -> GTG (M -> V); 6 G>A: CTG -> CTA (L -> L);
    ## 7 A>T: AAA -> TAA (K -> stop)
    whe <- codon_whe(c(1, 6, 7), c("A", "G", "A"), c("G", "A", "T"))
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))
    expect_identical(mc$GENEID, rep("geneA", 3L))
    expect_identical(mc$TXNAME, rep("cds-A", 3L))
    expect_identical(mc$REFCODON, c("ATG", "CTG", "AAA"))
    expect_identical(mc$VARCODON, c("GTG", "CTA", "TAA"))
    expect_identical(mc$REFAA, c("M", "L", "K"))
    expect_identical(mc$VARAA, c("V", "L", "*"))
    expect_identical(mc$CONSEQUENCE,
                     c("nonsynonymous", "synonymous", "nonsense"))
    expect_identical(mc$CDSLOC, c(1L, 6L, 7L))
    expect_identical(mc$PROTEINLOC, c(1L, 2L, 3L))
})

test_that("annotateCodonChange records sequence lengths from 'ref'", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(1, "A", "G")
    expect_true(is.na(GenomeInfoDb::seqlengths(
        SummarizedExperiment::rowRanges(whe))[["seg1"]]))
    out <- annotateCodonChange(whe, cds_gr(1, 15), ref)
    expect_identical(GenomeInfoDb::seqlengths(
        SummarizedExperiment::rowRanges(out))[["seg1"]], 15L)
})

test_that("annotateCodonChange stops when a recorded length disagrees", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(1, "A", "G")
    rr <- SummarizedExperiment::rowRanges(whe)
    GenomeInfoDb::seqlengths(rr) <- c(seg1 = 99L)
    SummarizedExperiment::rowRanges(whe) <- rr
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), ref),
                 "seg1 \\(99 vs 15\\)")
})

test_that("annotateCodonChange complements alleles on the minus strand", {
    skip_if_not_installed("txdbmaker")
    ## reverse complement of ATG CTG AAA GGG TAA
    ref <- Biostrings::DNAStringSet(c(seg1 = "TTACCCTTTCAGCAT"))
    ## plus-strand T>C at 15 is A>G at codon 1 on the minus strand
    whe <- codon_whe(15, "T", "C")
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15, strand = "-"),
                                       ref))
    expect_identical(mc$REFCODON, "ATG")
    expect_identical(mc$VARCODON, "GTG")
    expect_identical(mc$CONSEQUENCE, "nonsynonymous")
})

test_that("annotateCodonChange stops when sites and CDS share no seqname", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(NC_045512.2 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(1, "A", "G", seqnames = "MN908947.3")
    cds <- cds_gr(1, 15, seqnames = "NC_045512.2")
    expect_error(annotateCodonChange(whe, cds, ref), "MN908947\\.3")
    expect_error(annotateCodonChange(whe, cds, ref), "NC_045512\\.2")
})

test_that("annotateCodonChange matches FASTA names cut at the first space", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(
        "seg1 a description as in an NCBI header" = "ATGCTGAAAGGGTAA"))
    mc <- codon_mc(annotateCodonChange(codon_whe(1, "A", "G"),
                                       cds_gr(1, 15), ref))
    expect_identical(mc$VARAA, "V")
})

test_that("annotateCodonChange stops when 'ref' lacks a coding sequence", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(other = "ATGCTGAAAGGGTAA"))
    expect_error(annotateCodonChange(codon_whe(1, "A", "G"),
                                     cds_gr(1, 15), ref),
                 "no sequence named seg1")
})

test_that("annotateCodonChange skips a site whose ref allele is wrong", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    ## the reference base at 1 is A, not C
    whe <- codon_whe(c(1, 6), c("C", "G"), c("G", "A"))
    expect_warning(out <- annotateCodonChange(whe, cds_gr(1, 15), ref),
                   "seg1:1, C against A")
    mc <- codon_mc(out)
    expect_identical(mc$REFAA, c(NA, "L"))
})

test_that("annotateCodonChange leaves indels and multi-base changes NA", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(c(1, 4, 7), c("A", "CT", "A"), c("AT", "GA", "G"))
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))
    expect_identical(is.na(mc$REFAA), c(TRUE, TRUE, FALSE))
    expect_identical(mc$VARAA[3], "E")
})

test_that("annotateCodonChange leaves sites outside every CDS NA", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA",
                                      seg2 = "CCCCCCCCCC"))
    whe <- codon_whe(c(1, 5), c("A", "C"), c("G", "T"),
                     seqnames = c("seg1", "seg2"))
    expect_warning(out <- annotateCodonChange(whe, cds_gr(1, 15), ref),
                   "seg2")
    mc <- codon_mc(out)
    expect_identical(mc$GENEID, c("geneA", NA))
    expect_identical(GenomeInfoDb::seqlengths(
        SummarizedExperiment::rowRanges(out))[["seg2"]], 10L)
})

test_that("annotateCodonChange shortens the first CDS part by its phase", {
    skip_if_not_installed("txdbmaker")
    ## one leading base before the first complete codon (phase 1)
    ref <- Biostrings::DNAStringSet(c(seg1 = "CATGCTGAAA"))
    whe <- codon_whe(c(2, 7, 1), c("A", "G", "C"), c("G", "A", "T"))
    cds <- cds_gr(1, 10, phase = 1L)
    mc <- codon_mc(annotateCodonChange(whe, cds, ref))
    expect_identical(mc$REFCODON[1:2], c("ATG", "CTG"))
    expect_identical(mc$VARCODON[1:2], c("GTG", "CTA"))
    expect_identical(mc$CONSEQUENCE[2], "synonymous")
    ## the phase base is not part of any codon
    expect_true(is.na(mc$REFCODON[3]))
})

test_that("annotateCodonChange stops on a later CDS part out of frame", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    ## part 1 is 10 nt, so part 2 must have phase 2, not 0
    cds <- cds_gr(c(1, 21), c(10, 60), ID = "tx1", gene = "g1",
                  phase = c(0L, 0L))
    expect_error(annotateCodonChange(codon_whe(1, "A", "G"), cds, ref),
                 "shift the reading frame")
})

test_that("annotateCodonChange joins the parts of a frameshifted CDS", {
    skip_if_not_installed("txdbmaker")
    ## two parts sharing nucleotide 30, as at the ORF1ab -1 frameshift:
    ## 1-30 (10 codons) then 30-60 read on as codons 11-20 shifted by one
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    cds <- cds_gr(c(1, 30), c(30, 59), ID = "tx1", gene = "g1",
                  phase = c(0L, 0L))
    ## site 31 is the second base of part 2: CDS position 32, codon 11
    whe <- codon_whe(31, "C", "A")
    mc <- codon_mc(annotateCodonChange(whe, cds, ref))
    expect_identical(mc$CDSLOC, 32L)
    expect_identical(mc$PROTEINLOC, 11L)
    joined <- paste0(substr(codon_seq20, 1, 30), substr(codon_seq20, 30, 59))
    expect_identical(mc$REFCODON, substr(joined, 31, 33))
})

test_that("a nucleotide read twice by one transcript keeps its first reading (was: row order, reported as two transcripts)", {
    skip_if_not_installed("txdbmaker")
    ## parts 1-30 and 30-59 share nucleotide 30, which tx1 reads as the
    ## third base of codon 10 and the first base of codon 11
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    cds <- cds_gr(c(1, 30), c(30, 59), ID = "tx1", gene = "g1",
                  phase = c(0L, 0L))
    whe <- codon_whe(30, substr(codon_seq20, 30, 30), "C")
    msgs <- character()
    out <- withCallingHandlers(annotateCodonChange(whe, cds, ref),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_false(any(grepl("more than one transcript", msgs)))
    expect_true(any(grepl("1 site\\(s\\) lie on a nucleotide that the kept transcript reads twice",
                          msgs)))
    all <- codon_mc(annotateCodonChange(whe, cds, ref, collapse = "all"))
    expect_identical(as.list(all$CDSLOC), list(c(30L, 31L)))
    mc <- codon_mc(out)
    expect_identical(mc$CDSLOC, 30L)
    expect_identical(mc$CONSEQUENCE, as.list(all$CONSEQUENCE)[[1L]][1L])
    ## the reference annotator reported the second reading instead
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", codon_seq20), fa)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        paste("seg1", ".", "CDS", c(1, 30), c(30, 59), ".", "+", "0",
              "ID=tx1;gene=g1", sep = "\t")), gff)
    old <- refAnnotateCodonChange(whe, gff, fa)
    expect_identical(old$REF_CODON, as.list(all$REFCODON)[[1L]][2L])
    expect_false(identical(old$REF_CODON, mc$REFCODON))
})

test_that("annotateCodonChange keeps the longest transcript by default", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    ## tx2 lies inside tx1 in the same frame, starting at the ATG at 25
    cds <- cds_gr(c(1, 25), c(60, 60), ID = c("tx1", "tx2"),
                  gene = c("g1", "g2"))
    ## 4 C>A: codon 2 of tx1 only; 31 C>A: codon 11 of tx1, codon 3 of
    ## tx2; CTG -> ATG in both
    whe <- codon_whe(c(4, 31), c("C", "C"), c("A", "A"))
    expect_message(out <- annotateCodonChange(whe, cds, ref),
                   "1 site\\(s\\) fall in more than one transcript.*at 0 of")
    mc <- codon_mc(out)
    expect_identical(mc$TXNAME, c("tx1", "tx1"))
    expect_identical(mc$GENEID, c("g1", "g1"))
    expect_identical(mc$PROTEINLOC, c(2L, 11L))
})

test_that("collapse = 'all' keeps every transcript", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    cds <- cds_gr(c(1, 25), c(60, 60), ID = c("tx1", "tx2"),
                  gene = c("g1", "g2"))
    whe <- codon_whe(c(4, 31), c("C", "C"), c("A", "A"))
    mc <- codon_mc(annotateCodonChange(whe, cds, ref, collapse = "all"))
    expect_s4_class(mc$TXNAME, "CharacterList")
    expect_s4_class(mc$PROTEINLOC, "IntegerList")
    expect_identical(as.list(mc$TXNAME), list("tx1", c("tx1", "tx2")))
    expect_identical(as.list(mc$PROTEINLOC), list(2L, c(11L, 3L)))
    expect_identical(as.list(mc$CONSEQUENCE),
                     list("nonsynonymous", rep("nonsynonymous", 2L)))
})

test_that("a consequence that differs by transcript is reported", {
    skip_if_not_installed("txdbmaker")
    ## tx1 reads ATG GAT GCT GCT GCT; tx2, in the next frame, reads
    ## ATG CTG CTG. G>A at 10 is GCT -> ACT (A -> T) in tx1 and
    ## CTG -> CTA (L -> L) in tx2.
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGGATGCTGCTGCT"))
    cds <- cds_gr(c(1, 5), c(15, 13), ID = c("tx1", "tx2"),
                  gene = c("g1", "g2"))
    whe <- codon_whe(10, "G", "A")
    expect_message(out <- annotateCodonChange(whe, cds, ref),
                   "differs between transcripts at 1 of")
    expect_identical(codon_mc(out)$CONSEQUENCE, "nonsynonymous")
    all <- codon_mc(annotateCodonChange(whe, cds, ref, collapse = "all"))
    expect_identical(as.list(all$CONSEQUENCE),
                     list(c("nonsynonymous", "synonymous")))
    expect_identical(as.list(all$REFCODON), list(c("GCT", "CTG")))
})

test_that("annotateCodonChange reads an alternative start codon as M", {
    skip_if_not_installed("txdbmaker")
    ## Biostrings::translate() treats CTG at codon 1 as an initiator, so
    ## CTG -> ATG at the start of a coding sequence is synonymous (M -> M)
    ## although the same change elsewhere is L -> M.
    ref <- Biostrings::DNAStringSet(c(seg1 = "CTGCTGAAAGGGTAA"))
    whe <- codon_whe(c(1, 4), c("C", "C"), c("A", "A"))
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))
    expect_identical(mc$REFAA, c("M", "L"))
    expect_identical(mc$VARAA, c("M", "M"))
    expect_identical(mc$CONSEQUENCE, c("synonymous", "nonsynonymous"))
})

test_that("losing a stop codon is nonsynonymous, not nonsense", {
    skip_if_not_installed("txdbmaker")
    ## The removed in-house annotator called this change "Nonsense";
    ## predictCoding() calls it nonsynonymous. pN/pS counts it as
    ## nonsynonymous either way.
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(13, "T", "C")
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))
    expect_identical(mc$REFAA, "*")
    expect_identical(mc$VARAA, "Q")
    expect_identical(mc$CONSEQUENCE, "nonsynonymous")
})

test_that("annotateCodonChange accepts a TxDb and a GRangesList", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(c(1, 6), c("A", "G"), c("G", "A"))
    from_gr <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))

    grl <- GenomicRanges::GRangesList(cds_A = GenomicRanges::GRanges("seg1",
        IRanges::IRanges(1, 15), strand = "+"))
    S4Vectors::mcols(grl)$gene <- "geneA"
    from_grl <- codon_mc(annotateCodonChange(whe, grl, ref,
                                             geneCol = "gene"))
    expect_identical(from_grl$CONSEQUENCE, from_gr$CONSEQUENCE)
    expect_identical(from_grl$GENEID, from_gr$GENEID)

    txdb <- suppressWarnings(txdbmaker::makeTxDb(
        data.frame(tx_id = 1L, tx_name = "cds-A", tx_chrom = "seg1",
                   tx_strand = "+", tx_start = 1L, tx_end = 15L),
        data.frame(tx_id = 1L, exon_rank = 1L, exon_start = 1L,
                   exon_end = 15L, cds_start = 1L, cds_end = 15L),
        data.frame(tx_id = 1L, gene_id = "geneA"),
        chrominfo = data.frame(chrom = "seg1", length = 15L,
                               is_circular = FALSE)))
    from_txdb <- codon_mc(annotateCodonChange(whe, txdb, ref))
    cols <- c("TXNAME", "GENEID", "CONSEQUENCE", "REFCODON", "VARCODON",
              "REFAA", "VARAA", "CDSLOC", "PROTEINLOC")
    ## without this the loop would pass on two NULLs if a column were renamed
    expect_true(all(cols %in% colnames(from_gr)))
    expect_true(all(cols %in% colnames(from_txdb)))
    for (col in cols)
        expect_identical(from_txdb[[col]], from_gr[[col]], info = col)
})

test_that("annotateCodonChange stops on a TxDb with no coding sequences", {
    skip_if_not_installed("txdbmaker")
    txdb <- suppressWarnings(txdbmaker::makeTxDb(
        data.frame(tx_id = 1L, tx_name = "nc", tx_chrom = "seg1",
                   tx_strand = "+", tx_start = 1L, tx_end = 15L),
        data.frame(tx_id = 1L, exon_rank = 1L, exon_start = 1L,
                   exon_end = 15L),
        chrominfo = data.frame(chrom = "seg1", length = 15L,
                               is_circular = FALSE)))
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    expect_error(annotateCodonChange(codon_whe(1, "A", "G"), txdb, ref),
                 "no coding sequences")
})

test_that("annotateCodonChange rejects file paths and non-CDS records", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(1, "A", "G")
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), "ref.fa"),
                 "readDNAStringSet")
    expect_error(annotateCodonChange(whe, "genes.gff3", ref),
                 "rtracklayer::import")
    mixed <- cds_gr(c(1, 1), c(15, 15), ID = c("g", "c"), gene = "geneA",
                    type = c("gene", "CDS"))
    expect_error(annotateCodonChange(whe, mixed, ref), "type gene")
})

test_that("annotateCodonChange names missing columns and unused arguments", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(1, "A", "G")
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), ref,
                                     txCol = "protein_id"),
                 "protein_id not found.*Available columns: ID, gene")
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), ref,
                                     nameCol = "gene"),
                 "Unused argument")
    no_alt <- whe
    rr <- SummarizedExperiment::rowRanges(no_alt)
    S4Vectors::mcols(rr)$alt <- NULL
    SummarizedExperiment::rowRanges(no_alt) <- rr
    expect_error(annotateCodonChange(no_alt, cds_gr(1, 15), ref),
                 "'ref' and 'alt'")
})

test_that("annotateCodonChange reports transcript names with geneCol = NULL", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
                                  strand = "+", ID = "cds-A")
    mc <- codon_mc(annotateCodonChange(codon_whe(1, "A", "G"), cds, ref,
                                       geneCol = NULL))
    expect_identical(mc$GENEID, "cds-A")
})

test_that("the initiator rule applies only to the first base of a CDS", {
    skip_if_not_installed("txdbmaker")
    ## predictCoding() reports M for an alternative start codon only when
    ## the substitution is CDS position 1. At positions 2 and 3 of the same
    ## codon the standard table applies, so TTG reads as L.
    ref <- Biostrings::DNAStringSet(c(seg1 = "TTGCTGAAAGGGTAA"))
    whe <- codon_whe(c(1, 2, 3), c("T", "T", "G"), c("C", "C", "A"))
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref))
    expect_identical(mc$CDSLOC, 1:3)
    expect_identical(mc$REFCODON, rep("TTG", 3L))
    expect_identical(mc$VARCODON, c("CTG", "TCG", "TTA"))
    expect_identical(mc$REFAA, c("M", "L", "L"))
    expect_identical(mc$VARAA, c("M", "S", "L"))
    expect_identical(mc$CONSEQUENCE,
                     c("synonymous", "nonsynonymous", "synonymous"))
})

test_that("a codon with an ambiguous base is left unannotated", {
    skip_if_not_installed("txdbmaker")
    ## predictCoding() returns CONSEQUENCE "not translated" and empty
    ## amino acids for such a codon; counting that as synonymous would
    ## inflate the synonymous iSNVs of pN/pS.
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCNGAAAGGGTAA"))
    whe <- codon_whe(c(6, 7), c("G", "A"), c("A", "T"))
    expect_warning(out <- annotateCodonChange(whe, cds_gr(1, 15), ref),
                   "1 site\\(s\\) have a reference codon with an ambiguous")
    mc <- codon_mc(out)
    expect_identical(mc$CONSEQUENCE, c(NA, "nonsense"))
    expect_identical(mc$REFAA, c(NA, "K"))
    expect_identical(mc$VARAA, c(NA, "*"))
    expect_identical(mc$GENEID, c(NA, "geneA"))
})

test_that("annotateCodonChange shortens a minus-strand first part by phase", {
    skip_if_not_installed("txdbmaker")
    ## reverse complement of C + ATG CTG AAA GGG TAA; on the minus strand
    ## the phase base is the last one, at position 16
    ref <- Biostrings::DNAStringSet(c(seg1 = "TTACCCTTTCAGCATG"))
    whe <- codon_whe(c(15, 10, 16), c("T", "C", "G"), c("C", "T", "A"))
    cds <- cds_gr(1, 16, strand = "-", phase = 1L)
    mc <- codon_mc(annotateCodonChange(whe, cds, ref))
    expect_identical(mc$REFCODON[1:2], c("ATG", "CTG"))
    expect_identical(mc$VARCODON[1:2], c("GTG", "CTA"))
    expect_identical(mc$CONSEQUENCE[1:2],
                     c("nonsynonymous", "synonymous"))
    expect_true(is.na(mc$REFCODON[3]))
})

test_that("annotateCodonChange orders minus-strand parts by transcription", {
    skip_if_not_installed("txdbmaker")
    ## parts 12-17 then 1-9 on the minus strand join to ATG CTG AAA GGG TAA
    ref <- Biostrings::DNAStringSet(c(seg1 = "TTACCCTTTGGCAGCAT"))
    whe <- codon_whe(c(17, 9), c("T", "T"), c("C", "G"))
    cds <- cds_gr(c(12, 1), c(17, 9), strand = "-", ID = "tx1", gene = "g1",
                  phase = c(0L, 0L))
    mc <- codon_mc(annotateCodonChange(whe, cds, ref))
    expect_identical(mc$CDSLOC, c(1L, 7L))
    expect_identical(mc$REFCODON, c("ATG", "AAA"))
    expect_identical(mc$VARCODON, c("GTG", "CAA"))
    expect_identical(mc$PROTEINLOC, c(1L, 3L))
    ## part 2 follows 6 nucleotides, so its phase must be 0, not 1
    bad <- cds_gr(c(12, 1), c(17, 9), strand = "-", ID = "tx1", gene = "g1",
                  phase = c(0L, 1L))
    expect_error(annotateCodonChange(whe, bad, ref),
                 "shift the reading frame")
})

test_that("collapse = 'all' gives an empty element for an unannotated site", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(
        seg1 = paste0("ATGCTGAAAGGGTAA", "CCCCCCCCCCCCCCC")))
    whe <- codon_whe(c(1, 20), c("A", "C"), c("G", "T"))
    mc <- codon_mc(annotateCodonChange(whe, cds_gr(1, 15), ref,
                                       collapse = "all"))
    expect_identical(as.list(mc$CONSEQUENCE),
                     list("nonsynonymous", character()))
    expect_identical(as.list(mc$PROTEINLOC), list(1L, integer()))
})

test_that("ties in coding length go to the first transcript in 'cds'", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = codon_seq20))
    cds <- cds_gr(c(1, 1), c(60, 60), ID = c("txB", "txA"),
                  gene = c("gB", "gA"))
    whe <- codon_whe(4, "C", "A")
    mc <- codon_mc(suppressMessages(annotateCodonChange(whe, cds, ref)))
    expect_identical(mc$TXNAME, "txB")
    expect_identical(mc$GENEID, "gB")
})

test_that("a site past the end of its sequence gets one clear warning", {
    skip_if_not_installed("txdbmaker")
    ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
    whe <- codon_whe(c(1, 200), c("A", "A"), c("G", "G"))
    msgs <- character()
    out <- withCallingHandlers(annotateCodonChange(whe, cds_gr(1, 15), ref),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_length(msgs, 1L)
    expect_match(msgs, "1 site\\(s\\) lie beyond the end")
    expect_identical(codon_mc(out)$REFAA, c("M", NA))
})

test_that("annotateCodonChange explains a single DNAString or empty path", {
    ref <- Biostrings::DNAString("ATGCTGAAAGGGTAA")
    whe <- codon_whe(1, "A", "G")
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), ref),
                 "DNAStringSet\\(c\\(chr1 = ref\\)\\)")
    expect_error(annotateCodonChange(whe, cds_gr(1, 15), character()),
                 "reference\\.fa")
})

test_that("annotateCodonChange agrees with the removed in-house annotator", {
    skip_if_not_installed("txdbmaker")
    ## Every substitution at every position of a 65-nucleotide sequence,
    ## once against a plus-strand coding sequence and once against the
    ## same coding sequence on the minus strand of the reverse complement,
    ## compared with the reference implementation in
    ## helper-reference-codon.R.
    bases <- c("A", "C", "G", "T")
    plus_seq <- paste0(codon_seq20, "CCCCC")
    minus_seq <- as.character(Biostrings::reverseComplement(
        Biostrings::DNAString(plus_seq)))
    cases <- list(
        list(seq = plus_seq, from = 1L, to = 60L, strand = "+"),
        list(seq = minus_seq, from = 6L, to = 65L, strand = "-"))
    for (case in cases) {
        chars <- strsplit(case$seq, "")[[1L]]
        pos <- rep(seq_along(chars), each = 3L)
        refb <- chars[pos]
        altb <- unlist(lapply(chars, function(b) setdiff(bases, b)))
        whe <- codon_whe(pos, refb, altb)

        fa <- tempfile(fileext = ".fa")
        writeLines(c(">seg1", case$seq), fa)
        gff <- tempfile(fileext = ".gff3")
        writeLines(c("##gff-version 3",
            paste("seg1", ".", "CDS", case$from, case$to, ".", case$strand,
                  "0", "ID=tx1;gene=g1", sep = "\t")), gff)
        old <- refAnnotateCodonChange(whe, gff, fa)
        cds <- cds_gr(case$from, case$to, strand = case$strand, ID = "tx1",
                      gene = "g1")
        new <- codon_mc(annotateCodonChange(whe, cds,
                                            Biostrings::readDNAStringSet(fa)))

        info <- paste("strand", case$strand)
        expect_identical(new$REFCODON, old$REF_CODON, info = info)
        expect_identical(new$VARCODON, old$ALT_CODON, info = info)
        expect_identical(new$REFAA, old$REF_AA, info = info)
        expect_identical(new$GENEID, old$GFF_FEATURE, info = info)
        ## The two differ in two documented places: the loss of a stop
        ## codon is "Nonsense" in the old code and nonsynonymous in
        ## predictCoding(); and the first codon of a coding sequence is
        ## read as an initiator, so a change at that one base into CTG or
        ## TTG gives M rather than L. Two more differences need more than
        ## one coding sequence and are pinned in the next test and in the
        ## frameshift test above: overlapping records resolve to the last
        ## record in the old code and to the longest transcript in the new,
        ## and a nucleotide read twice at a frameshift keeps its second
        ## reading in the old code and its first in the new.
        stop_loss <- !is.na(old$REF_AA) & old$REF_AA == "*" &
            old$ALT_AA != "*"
        init <- !is.na(new$CDSLOC) & new$CDSLOC == 1L &
            new$VARCODON %in% c("CTG", "TTG")
        expect_identical(new$VARAA[!init], old$ALT_AA[!init], info = info)
        expect_identical(new$VARAA[init], rep("M", sum(init)), info = info)
        old_class <- tolower(old$AA_CLASS)
        old_class[stop_loss] <- "nonsynonymous"
        old_class[init] <- "synonymous"
        expect_identical(new$CONSEQUENCE, old_class, info = info)
        ## ATG reaches CTG and TTG only through its first base
        expect_identical(sum(init), 2L, info = info)
        ## TAA: 7 of its 9 substitutions remove the stop codon
        expect_identical(sum(stop_loss), 7L, info = info)
        ## the five bases outside the coding sequence are NA in both
        expect_identical(sum(is.na(new$REFAA)), 15L, info = info)
        expect_identical(is.na(new$REFAA), is.na(old$REF_AA), info = info)
    }
})

test_that("annotateCodonChange agrees with the reference annotator transcript by transcript where coding sequences overlap", {
    skip_if_not_installed("txdbmaker")
    ## tx1 (1-15) reads ATG GAT GCT GCT GCT; tx2 (5-13), one frame along,
    ## reads ATG CTG CTG. Every substitution at every base is compared.
    seq15 <- "ATGGATGCTGCTGCT"
    bases <- c("A", "C", "G", "T")
    chars <- strsplit(seq15, "")[[1L]]
    pos <- rep(seq_along(chars), each = 3L)
    whe <- codon_whe(pos, chars[pos],
                     unlist(lapply(chars, function(b) setdiff(bases, b))))
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", seq15), fa)
    ref <- Biostrings::readDNAStringSet(fa)
    gff_line <- function(from, to, id, gene)
        paste("seg1", ".", "CDS", from, to, ".", "+", "0",
              paste0("ID=", id, ";gene=", gene), sep = "\t")
    write_gff <- function(...) {
        f <- tempfile(fileext = ".gff3")
        writeLines(c("##gff-version 3", ...), f)
        f
    }
    cds <- cds_gr(c(1, 5), c(15, 13), ID = c("tx1", "tx2"),
                  gene = c("g1", "g2"))
    all <- codon_mc(annotateCodonChange(whe, cds, ref, collapse = "all"))
    tx_of <- as.list(all$TXNAME)
    for (tx in list(list(id = "tx1", gff = gff_line(1, 15, "tx1", "g1")),
                    list(id = "tx2", gff = gff_line(5, 13, "tx2", "g2")))) {
        old <- refAnnotateCodonChange(whe, write_gff(tx$gff), fa)
        k <- vapply(tx_of, function(t) match(tx$id, t), integer(1))
        reads <- !is.na(k)
        pick <- function(col) {
            l <- as.list(all[[col]])
            vapply(seq_along(k), function(i)
                if (reads[i]) l[[i]][k[i]] else NA_character_, "")
        }
        expect_identical(is.na(old$REF_CODON), !reads, info = tx$id)
        expect_identical(pick("REFCODON")[reads], old$REF_CODON[reads],
                         info = tx$id)
        expect_identical(pick("VARCODON")[reads], old$ALT_CODON[reads],
                         info = tx$id)
    }
    ## Where the two overlap they resolve differently, by design: the
    ## reference annotator keeps the last overlapping record, tx2, and
    ## collapse = "longest" keeps the longest transcript, tx1.
    old_both <- refAnnotateCodonChange(whe,
        write_gff(gff_line(1, 15, "tx1", "g1"), gff_line(5, 13, "tx2", "g2")),
        fa)
    new_longest <- codon_mc(suppressMessages(annotateCodonChange(whe, cds,
                                                                 ref)))
    overlap <- pos >= 5L & pos <= 13L
    expect_true(all(old_both$GFF_FEATURE[overlap] == "g2"))
    expect_true(all(new_longest$GENEID[overlap] == "g1"))
    expect_identical(old_both$GFF_FEATURE[!overlap],
                     new_longest$GENEID[!overlap])
})
