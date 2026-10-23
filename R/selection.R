# R/selection.R
# Within-host pN/pS and selection coefficients

#' @include methods.R annotate-codon.R
#' @importFrom SummarizedExperiment assay assayNames rowRanges
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomicFeatures cdsBy transcripts extractTranscriptSeqs
#' @importFrom GenomicRanges findOverlaps reduce intersect
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand
#' @importFrom IRanges width start end
NULL

## Nei and Gojobori (1986) synonymous sites of every sense codon under the
## standard genetic code (Biostrings::GENETIC_CODE): at each codon
## position, the fraction of the three possible substitutions that keep
## the amino acid. Substitutions to a stop codon count as nonsynonymous;
## the nonsynonymous sites of a codon are 3 minus its synonymous sites.
.NG_SYN_SITES <- local({
    code <- Biostrings::GENETIC_CODE
    bases <- c("A", "C", "G", "T")
    sense <- names(code)[code != "*"]
    vapply(sense, function(codon) {
        aa <- code[[codon]]
        s <- 0
        for (p in seq_len(3L)) {
            for (b in setdiff(bases, substr(codon, p, p))) {
                mutant <- codon
                substr(mutant, p, p) <- b
                if (code[[mutant]] == aa) s <- s + 1 / 3
            }
        }
        s
    }, numeric(1))
})

## Whether each change from a reference to a variant codon keeps the amino
## acid under the standard genetic code, with no initiator rule; NA where
## the reference codon is a stop codon or either codon is not one of the
## 64. These are the codons and the code that .NG_SYN_SITES counts sites
## over, so a change to a stop codon is nonsynonymous.
.ngSynonymous <- function(refcodon, varcodon) {
    code <- Biostrings::GENETIC_CODE
    ref_aa <- unname(code[toupper(as.character(refcodon))])
    var_aa <- unname(code[toupper(as.character(varcodon))])
    syn <- ref_aa == var_aa
    syn[is.na(ref_aa) | ref_aa == "*" | is.na(var_aa)] <- NA
    syn
}

## Codons of every coding sequence in a TxDb, one row per transcript and
## codon, keyed by the genomic positions of the codon's three bases. Two
## transcripts reading the same bases in the same frame give the same key,
## so the key identifies a codon of the genome rather than of a transcript.
.txCodons <- function(txdb, ref) {
    cds_by_tx <- cdsBy(txdb, by = "tx")
    if (length(cds_by_tx) == 0L)
        stop("'cds' contains no coding sequences.")
    tx <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
    tx_id <- as.character(tx$tx_id)
    tx_name <- as.character(tx$tx_name)
    tx_name[is.na(tx_name)] <- tx_id[is.na(tx_name)]
    gene <- vapply(as.list(tx$gene_id), function(g)
        if (length(g)) as.character(g[1L]) else NA_character_, character(1))
    names(tx_name) <- names(gene) <- tx_id

    seqs <- as.character(extractTranscriptSeqs(ref, cds_by_tx))
    rows <- lapply(names(cds_by_tx), function(id) {
        parts <- cds_by_tx[[id]]
        minus <- as.character(strand(parts))[1L] == "-"
        pos <- unlist(lapply(seq_along(parts), function(k) {
            if (minus) end(parts)[k]:start(parts)[k]
            else start(parts)[k]:end(parts)[k]
        }), use.names = FALSE)
        n_codons <- length(pos) %/% 3L
        if (n_codons == 0L) return(NULL)
        first <- 3L * seq_len(n_codons) - 2L
        chrom <- as.character(seqnames(parts))[1L]
        data.frame(
            tx_id = id, tx_name = tx_name[[id]], gene = gene[[id]],
            key = paste(chrom, if (minus) "-" else "+", pos[first],
                        pos[first + 1L], pos[first + 2L], sep = ":"),
            codon = toupper(substring(seqs[[id]], first, first + 2L)),
            stringsAsFactors = FALSE)
    })
    codons <- do.call(rbind, rows)
    if (is.null(codons))
        codons <- data.frame(tx_id = character(), tx_name = character(),
                             gene = character(), key = character(),
                             codon = character(), stringsAsFactors = FALSE)
    ## stop codons and codons with ambiguous bases carry no sites
    codons$S <- unname(.NG_SYN_SITES[codons$codon])
    list(codons = codons[!is.na(codons$S), , drop = FALSE],
         cds_by_tx = cds_by_tx, tx_name = tx_name, gene = gene)
}

## Pairs of transcripts whose coding sequences share nucleotides, with the
## number of nucleotides and of codons they share.
.cdsOverlaps <- function(tc) {
    grl <- tc$cds_by_tx
    spans <- unlist(range(grl), use.names = FALSE)
    ids <- names(grl)
    hits <- findOverlaps(spans, spans, ignore.strand = TRUE)
    q <- S4Vectors::queryHits(hits)
    s <- S4Vectors::subjectHits(hits)
    keep <- q < s
    q <- q[keep]
    s <- s[keep]
    if (!length(q))
        return(data.frame(tx1 = character(), tx2 = character(),
                          shared_nt = integer(), shared_codons = integer(),
                          stringsAsFactors = FALSE))
    keys <- split(tc$codons$key, factor(tc$codons$tx_id, levels = ids))
    shared_nt <- vapply(seq_along(q), function(k) {
        a <- reduce(grl[[ids[q[k]]]], ignore.strand = TRUE)
        b <- reduce(grl[[ids[s[k]]]], ignore.strand = TRUE)
        sum(width(intersect(a, b, ignore.strand = TRUE)))
    }, integer(1))
    shared_codons <- vapply(seq_along(q), function(k)
        length(intersect(keys[[ids[q[k]]]], keys[[ids[s[k]]]])), integer(1))
    label <- function(id) {
        g <- tc$gene[id]
        ifelse(is.na(g) | g == tc$tx_name[id], tc$tx_name[id],
               paste0(tc$tx_name[id], " (", g, ")"))
    }
    out <- data.frame(tx1 = unname(label(ids[q])),
                      tx2 = unname(label(ids[s])),
                      shared_nt = shared_nt, shared_codons = shared_codons,
                      stringsAsFactors = FALSE)
    out[out$shared_nt > 0L, , drop = FALSE]
}

# ============================================================
# pNpSWithinHost -- pN/pS ratio per sample
# ============================================================

#' Within-host pN/pS ratio
#'
#' Computes, for each sample, the ratio of nonsynonymous to synonymous
#' polymorphism per site among its iSNVs. This is a within-host
#' polymorphism ratio, pN/pS: it counts variants segregating inside one
#' host, not substitutions fixed between two sequences, so it is not the
#' dN/dS of a between-sequence comparison and carries no multiple-hit
#' correction. Sites are counted with the method of Nei and Gojobori
#' (1986). For dN/dS between hosts or between consensus sequences, use
#' the models in the \pkg{MSA2dist} package.
#'
#' @details
#' With \eqn{n_S} synonymous and \eqn{n_N} nonsynonymous iSNVs present in
#' a sample, and \eqn{S} synonymous and \eqn{N} nonsynonymous sites in
#' the coding sequences,
#' \deqn{p_S = n_S / S, \quad p_N = n_N / N, \quad
#'   \mathrm{pN/pS} = p_N / p_S.}
#' Under the assumptions of the counting model (every substitution
#' equally likely, no selection) the expected ratio is 1. A ratio below
#' 1 is what purifying selection against amino acid changes would give,
#' but a within-host ratio rests on few iSNVs, is highly variable, and
#' should always be read with its counts. It is \code{NA} when no
#' synonymous iSNV is present.
#'
#' \strong{Classification.} The counts use only the codons and the
#' genetic code that the sites are counted over. This function calls
#' \code{\link{annotateCodonChange}} on \code{x} with the same \code{cds}
#' and \code{ref}, so annotation columns already present in \code{x} are
#' not used, and takes from it, for each transcript reading an iSNV, the
#' reference and variant codon. A reading counts only when its reference
#' codon is a sense codon: an iSNV in a stop codon, including the loss of
#' a stop codon, is counted in neither \eqn{n_S} nor \eqn{n_N}, since a
#' stop codon carries no sites. The reading is then classified by
#' translating both codons with the standard genetic code: synonymous
#' when the amino acid is kept, counting towards \eqn{n_S}, and
#' nonsynonymous otherwise, a change to a stop codon included, counting
#' towards \eqn{n_N}. The initiator rule is not applied, so \code{ATG} to
#' \code{CTG} or \code{TTG} at the first base of a coding sequence is
#' nonsynonymous here, as \code{ATG} carries no synonymous site, although
#' \code{annotateCodonChange} reports it as synonymous.
#'
#' The totals take each iSNV by its first counted reading, in the order
#' of \code{collapse = "longest"}: the longest transcript that reads it in
#' a sense codon, and where that transcript reads the nucleotide twice,
#' its first reading. Each gene's counts take it by the first counted
#' reading of that gene's own transcripts, so an iSNV in a codon two
#' genes share counts in both genes, as the codon's sites do, but once in
#' the total. Sites outside every coding sequence, insertions and
#' deletions, and codons holding an ambiguous base are not counted. An
#' iSNV is present in a sample when its frequency is above zero and, with
#' \code{usePassedOnly = TRUE}, it passed QC; it is counted once whatever
#' its frequency.
#'
#' \strong{Sites.} \eqn{S} and \eqn{N} are summed over the distinct
#' codons of the genome that some coding sequence reads: a codon is
#' identified by the genomic positions of its three bases, so where two
#' coding sequences read the same bases in the same frame (ORF1a inside
#' ORF1ab in SARS-CoV-2) the shared codons are counted once, while
#' coding sequences that overlap in different frames (ORF7a and ORF7b)
#' contribute their codons separately. Each sense codon carries three
#' sites, split between synonymous and nonsynonymous under the standard
#' genetic code; stop codons and codons with an ambiguous base carry
#' none, so \eqn{S + N} is three times the number of distinct sense
#' codons. When coding sequences of different transcripts overlap, a
#' warning names them with the nucleotides and codons they share. The
#' sites of a gene are summed in the same way over the codons of that
#' gene's own transcripts, so a codon shared by two genes is counted in
#' each of them but once in the total.
#'
#' @param x A \code{\link{WithinHostExperiment}} with \code{ref} and
#'   \code{alt} columns in \code{mcols(rowRanges(x))}.
#' @param cds The coding sequences, as for
#'   \code{\link{annotateCodonChange}}: a \code{TxDb}, a \code{GRanges} of
#'   CDS records, or a \code{GRangesList} with one element per transcript.
#' @param ref A \code{DNAStringSet} with the reference sequence of each
#'   seqlevel carrying coding sequences, for example from
#'   \code{Biostrings::readDNAStringSet()}.
#' @param usePassedOnly Logical. If \code{TRUE} (default), only iSNVs
#'   with \code{qcPass == TRUE} are counted.
#' @param txCol,geneCol As for \code{\link{annotateCodonChange}}: the
#'   transcript and gene columns of a \code{GRanges} or
#'   \code{GRangesList} \code{cds}.
#' @param ... Not used; supplying anything is an error.
#'
#' @return A \code{DataFrame} with one row per sample and gene carrying
#'   at least one counted iSNV (one row with \code{gene = NA} when a
#'   sample has none):
#'   \describe{
#'     \item{sample_id}{Sample identifier.}
#'     \item{nS, nN}{Synonymous and nonsynonymous iSNVs in the sample.}
#'     \item{S_sites, N_sites}{Synonymous and nonsynonymous sites over
#'       all coding sequences.}
#'     \item{pS, pN}{The proportions \eqn{n_S/S} and \eqn{n_N/N}.}
#'     \item{pNpS}{\eqn{p_N/p_S}; \code{NA} when \eqn{n_S = 0}.}
#'     \item{gene, gene_nS, gene_nN, gene_S_sites, gene_N_sites,
#'       gene_pNpS}{The same quantities for each gene.}
#'   }
#'
#' @references
#' Nei M, Gojobori T (1986). Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol Biol Evol} 3:418-426.
#' \doi{10.1093/oxfordjournals.molbev.a040410}
#'
#' @seealso \code{\link{annotateCodonChange}} for the classification of
#'   each iSNV.
#'
#' @name pNpSWithinHost
#' @aliases pNpSWithinHost,WithinHostExperiment,TxDb,DNAStringSet-method
#'   pNpSWithinHost,WithinHostExperiment,GRangesList,DNAStringSet-method
#'   pNpSWithinHost,WithinHostExperiment,GRanges,DNAStringSet-method
#'   pNpSWithinHost,WithinHostExperiment,ANY,character-method
#'   pNpSWithinHost,WithinHostExperiment,ANY,XString-method
#'   pNpSWithinHost,WithinHostExperiment,character,DNAStringSet-method
#' @examples
#' ## A 15-nucleotide coding sequence, ATG CTG AAA GGG TAA
#' ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
#' cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
#'     strand = "+", ID = "cds-A", gene = "geneA")
#' ## 1 A>G is ATG -> GTG (nonsynonymous); 6 G>A is CTG -> CTA (synonymous)
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(1, 6), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "G")
#' S4Vectors::mcols(gr)$alt <- c("G", "A")
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' if (requireNamespace("txdbmaker", quietly = TRUE)) {
#'     pNpSWithinHost(whe, cds, ref)
#' }
NULL

## Count iSNVs and sites; x is already annotated by annotateCodonChange().
.pNpSCount <- function(x, tc, usePassedOnly) {
    codons <- tc$codons
    if (nrow(codons) == 0L)
        stop("The coding sequences in 'cds' hold no complete sense codon ",
             "on the sequences in 'ref'.")
    total <- codons[!duplicated(codons$key), , drop = FALSE]
    S_total <- sum(total$S)
    N_total <- sum(3 - total$S)
    by_gene <- codons[!is.na(codons$gene), , drop = FALSE]
    by_gene <- by_gene[!duplicated(by_gene[, c("gene", "key")]), ,
                       drop = FALSE]
    gene_S <- tapply(by_gene$S, by_gene$gene, sum)
    gene_N <- tapply(3 - by_gene$S, by_gene$gene, sum)

    ## x was annotated with collapse = "all": one element per transcript
    ## and reading, the longest transcript's first reading first. A
    ## reading counts only at a sense codon, the codons that carry sites,
    ## and is classified by the standard code without the initiator rule,
    ## the code the sites are counted under
    mc <- mcols(rowRanges(x))
    n_rows <- lengths(mc$CONSEQUENCE)
    readings <- data.frame(
        site = rep(seq_along(n_rows), n_rows),
        gene = as.character(unlist(mc$GENEID, use.names = FALSE)),
        syn = .ngSynonymous(unlist(mc$REFCODON, use.names = FALSE),
                            unlist(mc$VARCODON, use.names = FALSE)),
        stringsAsFactors = FALSE)
    readings <- readings[!is.na(readings$syn), , drop = FALSE]
    ## the total takes a site by its first counted reading
    first <- readings[!duplicated(readings$site), , drop = FALSE]
    is_syn <- is_non <- logical(length(n_rows))
    is_syn[first$site] <- first$syn
    is_non[first$site] <- !first$syn
    ## for each gene a site counts by the first counted reading of that
    ## gene's transcripts, so a codon two genes share counts in both
    by_tx <- readings[!is.na(readings$gene) &
                          !duplicated(readings[, c("site", "gene")]), ,
                      drop = FALSE]
    tx_syn <- by_tx$syn
    tx_non <- !by_tx$syn

    freq_mat <- assay(x, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(x)) {
        assay(x, "qcPass")
    } else {
        NULL
    }
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids))
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    ratio <- function(nS, nN, S, N)
        if (nS == 0L) NA_real_ else (nN / N) / (nS / S)

    results <- lapply(seq_len(ncol(freq_mat)), function(j) {
        freq_j <- freq_mat[, j]
        if (!is.null(qc_mat)) freq_j[!qc_mat[, j]] <- NA
        present <- !is.na(freq_j) & freq_j > 0
        nS <- sum(present & is_syn)
        nN <- sum(present & is_non)
        base_row <- data.frame(
            sample_id = sample_ids[j], nS = nS, nN = nN,
            S_sites = S_total, N_sites = N_total,
            pS = nS / S_total, pN = nN / N_total,
            pNpS = ratio(nS, nN, S_total, N_total),
            stringsAsFactors = FALSE)
        counted <- present[by_tx$site] & (tx_syn | tx_non)
        genes <- sort(unique(by_tx$gene[counted]))
        if (length(genes) == 0L) {
            return(cbind(base_row, data.frame(
                gene = NA_character_, gene_nS = NA_integer_,
                gene_nN = NA_integer_, gene_S_sites = NA_real_,
                gene_N_sites = NA_real_, gene_pNpS = NA_real_,
                stringsAsFactors = FALSE)))
        }
        gene_rows <- lapply(genes, function(g) {
            in_gene <- counted & by_tx$gene == g
            gnS <- sum(in_gene & tx_syn)
            gnN <- sum(in_gene & tx_non)
            gS <- unname(gene_S[g])
            gN <- unname(gene_N[g])
            data.frame(gene = g, gene_nS = gnS, gene_nN = gnN,
                       gene_S_sites = gS, gene_N_sites = gN,
                       gene_pNpS = ratio(gnS, gnN, gS, gN),
                       stringsAsFactors = FALSE)
        })
        cbind(base_row[rep(1L, length(genes)), , drop = FALSE],
              do.call(rbind, gene_rows))
    })
    out <- do.call(rbind, results)
    rownames(out) <- NULL
    DataFrame(out)
}

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "TxDb", "DNAStringSet"),
    function(x, cds, ref, usePassedOnly = TRUE, ...) {
        .stop_on_unused_dots("pNpSWithinHost", ...)
        if (!is.logical(usePassedOnly) || length(usePassedOnly) != 1L ||
            is.na(usePassedOnly))
            stop("'usePassedOnly' must be TRUE or FALSE.")
        ref <- .prepareReference(ref)
        x <- annotateCodonChange(x, cds, ref, collapse = "all")
        tc <- .txCodons(cds, ref)
        ov <- .cdsOverlaps(tc)
        if (nrow(ov)) {
            shown <- utils::head(ov, 5L)
            warning("Coding sequences of different transcripts overlap; ",
                    "codons they share are counted once in the total ",
                    "sites:\n", paste0("  ", shown$tx1, " / ", shown$tx2,
                    ": ", shown$shared_nt, " nt, ", shown$shared_codons,
                    " codon(s) shared", collapse = "\n"),
                    if (nrow(ov) > 5L)
                        paste0("\n  and ", nrow(ov) - 5L, " more pair(s)"),
                    call. = FALSE)
        }
        .pNpSCount(x, tc, usePassedOnly)
    })

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "GRangesList", "DNAStringSet"),
    function(x, cds, ref, usePassedOnly = TRUE, geneCol = NULL, ...) {
        ref <- .prepareReference(ref)
        txdb <- .cdsListToTxDb(cds, geneCol, ref)
        pNpSWithinHost(x, txdb, ref, usePassedOnly = usePassedOnly, ...)
    })

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "GRanges", "DNAStringSet"),
    function(x, cds, ref, usePassedOnly = TRUE, txCol = "ID",
             geneCol = "gene", ...) {
        grl <- .cdsRangesToList(cds, txCol, geneCol)
        pNpSWithinHost(x, grl, ref, usePassedOnly = usePassedOnly,
                       geneCol = if (is.null(geneCol)) NULL else "gene",
                       ...)
    })

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "ANY", "character"),
    function(x, cds, ref, ...) .refPathError(ref))

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "ANY", "XString"),
    function(x, cds, ref, ...) .xstringRefError(ref))

#' @rdname pNpSWithinHost
#' @export
setMethod("pNpSWithinHost",
    c("WithinHostExperiment", "character", "DNAStringSet"),
    function(x, cds, ref, ...) .cdsPathError(cds))


# ============================================================
# estimateSelectionCoefficient
# ============================================================

#' Estimate selection coefficient from frequency trajectory
#'
#' Given a vector of allele frequencies observed over equally-spaced
#' time points, estimates the selection coefficient \eqn{s} using
#' the change in logit-transformed frequency per unit time:
#' \deqn{s = \frac{\Delta \mathrm{logit}(p)}{\Delta t}}
#' where \eqn{\mathrm{logit}(p) = \log(p / (1 - p))}.
#'
#' This assumes a simple deterministic selection model. Each
#' consecutive pair of time points yields one estimate of \eqn{s}.
#'
#' @param freqs Numeric vector. Allele frequencies at successive
#'   time points (values must be in (0, 1)).
#' @param times Numeric vector (optional). Time values for each
#'   observation. If \code{NULL}, integer time steps 1, 2, ... are
#'   used.
#' @param generation_time Numeric scalar. Time per generation
#'   (default: 1). Used to scale the selection coefficient.
#'
#' @return Numeric vector of selection coefficient estimates, one
#'   per consecutive time interval. Length is
#'   \code{length(freqs) - 1}.
#'
#' @references
#' Feder AF et al. (2016). More effective drugs lead to harder
#' selective sweeps in the evolution of drug resistance in HIV-1.
#' \emph{eLife} 5:e10670. \doi{10.7554/eLife.10670}
#'
#' @export
#' @examples
#' freqs <- c(0.05, 0.10, 0.25, 0.50)
#' estimateSelectionCoefficient(freqs)
#' estimateSelectionCoefficient(freqs, times = c(0, 3, 7, 14))
estimateSelectionCoefficient <- function(freqs, times = NULL,
                                         generation_time = 1) {
    if (!is.numeric(freqs) || length(freqs) < 2L)
        stop("'freqs' must be a numeric vector with at least 2 elements.")
    if (any(is.na(freqs)))
        stop("'freqs' must not contain NA values.")
    if (any(freqs <= 0 | freqs >= 1))
        stop("'freqs' must contain values strictly between 0 and 1.")
    if (!is.numeric(generation_time) || length(generation_time) != 1L ||
        generation_time <= 0)
        stop("'generation_time' must be a single positive number.")

    if (is.null(times)) {
        times <- seq_along(freqs)
    } else {
        if (!is.numeric(times) || length(times) != length(freqs))
            stop("'times' must be a numeric vector with the same ",
                 "length as 'freqs'.")
    }

    logit_p <- log(freqs / (1 - freqs))
    delta_logit <- diff(logit_p)
    delta_t <- diff(times) / generation_time

    if (any(delta_t <= 0))
        stop("'times' must contain strictly increasing values.")

    delta_logit / delta_t
}
