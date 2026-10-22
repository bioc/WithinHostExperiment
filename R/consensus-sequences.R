# R/consensus-sequences.R
# Per-sample consensus sequences, for comparisons between hosts

#' @include annotate-codon.R
#' @importFrom Biostrings DNAStringSet replaceAt extractAt reverseComplement
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges width start end
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom SummarizedExperiment assay assayNames rowRanges
NULL

#' Consensus sequence of each sample
#'
#' Builds, for every sample of a \code{\link{WithinHostExperiment}}, the
#' reference sequence with each substitution present at or above
#' \code{minFreq} written in, and returns the sequences as a
#' \code{DNAStringSet}. Its purpose is comparison between hosts: the
#' within-host ratio of \code{\link{pNpSWithinHost}} describes
#' polymorphism inside one host, whereas dN/dS between hosts is a
#' divergence between their consensus sequences, which the models of
#' \code{MSA2dist::dnastring2kaks()} estimate.
#'
#' @details
#' Only substitutions are written in: single-nucleotide changes and
#' multi-nucleotide changes whose reference and alternative alleles have
#' the same length. Insertions and deletions are not applied, with a
#' warning, so every consensus has the length of the reference and stays
#' aligned to it and to every other sample's consensus; that is what
#' \code{MSA2dist::dnastring2kaks()} needs. Where a deletion or
#' insertion is fixed in a host, the returned sequence therefore
#' differs from that host's true consensus.
#'
#' A site enters a sample's consensus when its \code{altFreq} is at
#' least \code{minFreq} and, with \code{usePassedOnly = TRUE}, it passed
#' QC. With \code{minFreq} at or below 0.5 two alleles can qualify at
#' one position; the more frequent is written in, the first in the
#' object on a tie, and a warning reports how many were set aside. A
#' site whose \code{ref} allele differs from the reference at that
#' position is never written in, with a warning, since it was called
#' against another sequence.
#'
#' \code{region} restricts the output to part of the sequence, such as
#' one coding sequence: its ranges are taken in order along the
#' sequence and joined, and the result is reverse-complemented when
#' they lie on the minus strand, so a coding sequence split into parts
#' comes back as one sequence read in its own direction. Positions are
#' those of the reference, since no length-changing edit is made.
#'
#' @param x A \code{\link{WithinHostExperiment}} with \code{ref} and
#'   \code{alt} columns in \code{mcols(rowRanges(x))}.
#' @param ref A \code{DNAStringSet} holding the reference sequence the
#'   variants were called against, for example from
#'   \code{Biostrings::readDNAStringSet()}. Sequence names are cut at the
#'   first space.
#' @param minFreq Numeric scalar in (0, 1]. The frequency at or above
#'   which a substitution is written in (default 0.5).
#' @param usePassedOnly Logical. If \code{TRUE} (default), only sites
#'   with \code{qcPass == TRUE} in that sample are written in.
#' @param seqname The sequence of \code{ref} to build; required when
#'   \code{ref} holds more than one sequence and \code{region} does not
#'   name one.
#' @param region Optional \code{GRanges} on \code{seqname}, all on one
#'   strand: the part of the sequence to return; see Details.
#' @param ... Not used; supplying anything is an error.
#'
#' @return A \code{DNAStringSet} with one sequence per sample, named by
#'   sample. Its \code{mcols()} hold \code{n_substitutions}, the number
#'   of substitutions written into the returned part of each sequence.
#'
#' @seealso \code{\link{pNpSWithinHost}} for the within-host ratio;
#'   \code{MSA2dist::dnastring2kaks()} for dN/dS between sequences.
#'
#' @name consensusSequences
#' @aliases consensusSequences,WithinHostExperiment,DNAStringSet-method
#'   consensusSequences,WithinHostExperiment,character-method
#'   consensusSequences,WithinHostExperiment,XString-method
#' @examples
#' ## A coding sequence of 15 codons and three hosts, each carrying one
#' ## synonymous and one nonsynonymous change at high frequency:
#' ## host1 6 G>A (CTG -> CTA) and 7 A>G (AAA -> GAA),
#' ## host2 12 G>A (GGG -> GGA) and 13 T>C (TTT -> CTT),
#' ## host3 18 C>A (CCC -> CCA) and 19 A>C (AAA -> CAA)
#' ref <- Biostrings::DNAStringSet(c(
#'     seg1 = "ATGCTGAAAGGGTTTCCCAAACTGATGGTTCTGAAAGGGTTTTAA"))
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(6, 7, 12, 13, 18, 19), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("G", "A", "G", "T", "C", "A")
#' S4Vectors::mcols(gr)$alt <- c("A", "G", "A", "C", "A", "C")
#' freqs <- cbind(host1 = c(0.95, 0.91, 0.02, 0.03, 0.01, 0.02),
#'                host2 = c(0.03, 0.02, 0.88, 0.93, 0.04, 0.01),
#'                host3 = c(0.01, 0.04, 0.02, 0.01, 0.97, 0.90))
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freqs),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = colnames(freqs),
#'                                    row.names = colnames(freqs)))
#' cons <- consensusSequences(whe, ref)
#' cons
#'
#' ## dN/dS between the hosts' consensus coding sequences: every pair
#' ## differs by two synonymous and two nonsynonymous changes
#' if (requireNamespace("MSA2dist", quietly = TRUE)) {
#'     kaks <- suppressMessages(MSA2dist::dnastring2kaks(cons, model = "Li"))
#'     kaks[, c("seq1", "seq2", "Ka", "Ks", "Ka/Ks")]
#' }
NULL

## Keep, among one sample's candidate substitutions, a disjoint set:
## the more frequent of two overlapping edits wins, the earlier on a tie.
.disjointByFreq <- function(starts, widths, freqs) {
    o <- order(-freqs, seq_along(freqs))
    kept <- integer()
    for (i in o) {
        clash <- length(kept) &&
            any(starts[kept] <= starts[i] + widths[i] - 1L &
                starts[i] <= starts[kept] + widths[kept] - 1L)
        if (!clash) kept <- c(kept, i)
    }
    sort(kept)
}

#' @rdname consensusSequences
#' @export
setMethod("consensusSequences", c("WithinHostExperiment", "DNAStringSet"),
    function(x, ref, minFreq = 0.5, usePassedOnly = TRUE, seqname = NULL,
             region = NULL, ...) {
        .stop_on_unused_dots("consensusSequences", ...)
        if (!is.numeric(minFreq) || length(minFreq) != 1L ||
            is.na(minFreq) || minFreq <= 0 || minFreq > 1)
            stop("'minFreq' must be a single number in (0, 1].")
        if (!is.logical(usePassedOnly) || length(usePassedOnly) != 1L ||
            is.na(usePassedOnly))
            stop("'usePassedOnly' must be TRUE or FALSE.")
        ref <- .prepareReference(ref)
        rr <- rowRanges(x)
        mc <- mcols(rr)
        if (!all(c("ref", "alt") %in% colnames(mc)))
            stop("mcols(rowRanges(x)) must contain 'ref' and 'alt' columns.")

        if (!is.null(region)) {
            if (!is(region, "GRanges") || length(region) == 0L)
                stop("'region' must be a non-empty GRanges.")
            region_seq <- unique(as.character(seqnames(region)))
            region_strand <- unique(as.character(BiocGenerics::strand(region)))
            if (length(region_seq) != 1L)
                stop("'region' must lie on one sequence.")
            if (length(region_strand) != 1L)
                stop("'region' must lie on one strand.")
            if (is.null(seqname)) seqname <- region_seq
            if (!identical(region_seq, seqname))
                stop("'region' lies on ", region_seq, ", not on '",
                     seqname, "'.")
        }
        if (is.null(seqname)) {
            if (length(ref) != 1L)
                stop("'ref' holds ", length(ref), " sequences; choose one ",
                     "with 'seqname': ", paste(names(ref), collapse = ", "),
                     ".")
            seqname <- names(ref)
        }
        if (!is.character(seqname) || length(seqname) != 1L ||
            !seqname %in% names(ref))
            stop("'seqname' must name one sequence of 'ref': ",
                 paste(names(ref), collapse = ", "), ".")
        .checkSeqlevels(rr, GRanges(names(ref),
                                    IRanges(1L, width(ref))), "ref")
        genome <- ref[[seqname]]
        len <- length(genome)
        if (!is.null(region) &&
            (any(start(region) < 1L) || any(end(region) > len)))
            stop("'region' extends beyond the ", len, " bases of '",
                 seqname, "'.")

        ref_allele <- toupper(as.character(mc$ref))
        alt_allele <- toupper(as.character(mc$alt))
        on_seq <- as.character(seqnames(rr)) == seqname
        pos <- start(rr)
        w <- nchar(ref_allele)
        known <- on_seq & !is.na(ref_allele) & !is.na(alt_allele)
        acgt <- known & grepl("^[ACGT]+$", ref_allele) &
            grepl("^[ACGT]+$", alt_allele)
        substitution <- acgt & nchar(alt_allele) == w &
            width(rr) == w & pos + w - 1L <= len
        base <- rep(NA_character_, length(rr))
        if (any(substitution))
            base[substitution] <- as.character(extractAt(genome,
                IRanges(pos[substitution], width = w[substitution])))
        mismatch <- substitution & base != ref_allele
        usable <- substitution & !mismatch

        freq <- assay(x, "altFreq")
        qc <- if (usePassedOnly && "qcPass" %in% assayNames(x))
            assay(x, "qcPass") else NULL
        present <- !is.na(freq) & freq >= minFreq
        if (!is.null(qc)) present <- present & !is.na(qc) & qc
        rows_present <- rowSums(present) > 0L

        n_mismatch <- sum(mismatch & rows_present)
        if (n_mismatch) {
            k <- which(mismatch & rows_present)[1L]
            warning(n_mismatch, " site(s) carry a 'ref' allele that differs ",
                    "from 'ref' (first: ", seqname, ":", pos[k], ", ",
                    ref_allele[k], " against ", base[k], ") and are not ",
                    "written in.", call. = FALSE)
        }
        skipped <- on_seq & rows_present & !substitution
        if (any(skipped))
            warning(sum(skipped), " call(s) at or above minFreq are ",
                    "insertions, deletions, alleles other than A, C, G or ",
                    "T, or run past the sequence end, and are not written ",
                    "in; each consensus keeps the reference length.",
                    call. = FALSE)

        sample_ids <- colnames(freq)
        if (is.null(sample_ids))
            sample_ids <- paste0("S", seq_len(ncol(freq)))
        region_ir <- if (is.null(region)) IRanges(1L, len) else
            IRanges(start(region), end(region))
        region_ir <- region_ir[order(start(region_ir))]
        minus <- !is.null(region) && region_strand == "-"

        set_aside <- 0L
        seqs <- vector("list", ncol(freq))
        n_sub <- integer(ncol(freq))
        for (j in seq_len(ncol(freq))) {
            cand <- which(usable & present[, j])
            if (length(cand) > 1L) {
                kept <- .disjointByFreq(pos[cand], w[cand], freq[cand, j])
                set_aside <- set_aside + length(cand) - length(kept)
                cand <- cand[kept]
            }
            edited <- if (length(cand)) {
                replaceAt(genome, IRanges(pos[cand], width = w[cand]),
                          DNAStringSet(alt_allele[cand]))
            } else genome
            part <- unlist(extractAt(edited, region_ir))
            if (minus) part <- reverseComplement(part)
            seqs[[j]] <- as.character(part)
            n_sub[j] <- if (length(cand)) {
                sum(IRanges::overlapsAny(
                    IRanges(pos[cand], width = w[cand]), region_ir))
            } else 0L
        }
        if (set_aside)
            warning(set_aside, " substitution(s) overlapped a more frequent ",
                    "one in the same sample and were set aside.",
                    call. = FALSE)
        out <- DNAStringSet(unlist(seqs))
        names(out) <- sample_ids
        mcols(out) <- S4Vectors::DataFrame(n_substitutions = n_sub)
        out
    })

#' @rdname consensusSequences
#' @export
setMethod("consensusSequences", c("WithinHostExperiment", "character"),
    function(x, ref, ...) .refPathError(ref))

#' @rdname consensusSequences
#' @export
setMethod("consensusSequences", c("WithinHostExperiment", "XString"),
    function(x, ref, ...) .xstringRefError(ref))
