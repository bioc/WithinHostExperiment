# R/annotate-codon.R
# Codon and amino acid consequence of each SNV, from
# VariantAnnotation::predictCoding()

#' @include annotate-genes.R
#' @importFrom VariantAnnotation predictCoding
#' @importFrom GenomicFeatures cdsBy transcripts
#' @importFrom Biostrings DNAStringSet subseq
#' @importClassesFrom Biostrings DNAStringSet XString
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevelsInUse seqlengths seqlengths<- seqnames
#'   seqlevels seqlevels<-
#' @importFrom GenomicRanges GRanges granges
#' @importFrom IRanges IRanges CharacterList IntegerList width start end
#'   start<- end<-
#' @importFrom S4Vectors mcols mcols<- split
#' @importFrom SummarizedExperiment rowRanges
#' @importClassesFrom GenomicRanges GRanges GRangesList
NULL

## Reference sequences: names cut at the first space, as FASTA readers
## and samtools do, so that a header such as
## ">NC_045512.2 Severe acute respiratory syndrome..." gives the sequence
## name used in annotation files.
.prepareReference <- function(ref) {
    if (length(ref) == 0L)
        stop("'ref' contains no sequences.")
    nms <- names(ref)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))
        stop("Every sequence in 'ref' needs a name matching a seqlevel ",
             "of the variants.")
    names(ref) <- sub("[[:space:]].*$", "", nms)
    if (anyDuplicated(names(ref)))
        stop("Sequence names in 'ref' are not unique after cutting them ",
             "at the first space: ",
             paste(unique(names(ref)[duplicated(names(ref))]),
                   collapse = ", "))
    ref
}

## Record the length of each reference sequence in seqlengths(rr), for
## seqlevels that have none, and add the reference sequences that carry no
## call as seqlevels, since the genome length sums over every sequence. A
## recorded length that disagrees with the reference is an error: the
## variants were called against something else.
.fillSeqlengths <- function(rr, ref) {
    lens <- width(ref)
    names(lens) <- names(ref)
    extra <- setdiff(names(lens), seqlevels(rr))
    if (length(extra))
        seqlevels(rr) <- c(seqlevels(rr), extra)
    known <- seqlengths(rr)
    shared <- intersect(names(known), names(lens))
    clash <- shared[!is.na(known[shared]) & known[shared] != lens[shared]]
    if (length(clash))
        stop("Sequence length disagrees between the variants and 'ref' ",
             "for ", paste0(clash, " (", known[clash], " vs ",
                            lens[clash], ")", collapse = ", "), ".")
    fill <- shared[is.na(known[shared])]
    if (length(fill)) {
        known[fill] <- lens[fill]
        ## Recording a length makes a site past the end of its sequence
        ## out of bound. The caller warns about those in its own words.
        withCallingHandlers({
            seqlengths(rr) <- known
        }, warning = function(w) {
            if (grepl("out-of-bound range", conditionMessage(w),
                      fixed = TRUE))
                invokeRestart("muffleWarning")
        })
    }
    rr
}

.stopNoRef <- function(missing_names, ref) {
    stop("'ref' has no sequence named ",
         paste(missing_names, collapse = ", "),
         ", which carries coding sequences. Sequences in 'ref': ",
         paste(names(ref), collapse = ", "), ".", call. = FALSE)
}

## The three messages the path and single-sequence guards share, so that
## annotateCodonChange() and pNpSWithinHost() say the same thing.
.refPathError <- function(ref) {
    path <- if (length(ref)) ref[1L] else "reference.fa"
    stop("'ref' must be a DNAStringSet, not a file path. Read the ",
         "FASTA file first, for example ",
         "ref <- Biostrings::readDNAStringSet(\"", path, "\").",
         call. = FALSE)
}

.xstringRefError <- function(ref) {
    stop("'ref' must be a DNAStringSet with one named sequence per ",
         "seqlevel, not a single ", class(ref)[1L], ". Wrap it, for ",
         "example ref <- Biostrings::DNAStringSet(c(chr1 = ref)).",
         call. = FALSE)
}

.cdsPathError <- function(cds) {
    path <- if (length(cds)) cds[1L] else "annotation.gff3"
    stop("'cds' must be a GRanges, GRangesList or TxDb, not a file ",
         "path. Read a GFF3 or GTF file with rtracklayer::import() and ",
         "pass its CDS records, for example ",
         "features <- rtracklayer::import(\"", path, "\"); ",
         "cds <- features[features$type == \"CDS\"].", call. = FALSE)
}

#' Codon and amino acid consequence of each SNV
#'
#' Annotates every single-nucleotide substitution in a
#' \code{\link{WithinHostExperiment}} with the codon it changes and the
#' resulting amino acid, by passing the sites to
#' \code{VariantAnnotation::predictCoding()}. The coding sequences are
#' supplied either as a \code{TxDb} or as the CDS records of an
#' annotation in a \code{GRanges} or \code{GRangesList}; the reference
#' sequence is a \code{DNAStringSet}, as returned by
#' \code{Biostrings::readDNAStringSet()}.
#'
#' @details
#' \strong{Coding sequences.} \code{predictCoding()} reads coding
#' sequences from a \code{TxDb}. A \code{GRanges} of CDS records, for
#' example the rows of \code{rtracklayer::import()} with
#' \code{type == "CDS"}, is grouped into transcripts by the column
#' \code{txCol}: in GFF3 files from NCBI and Ensembl the CDS records of
#' one protein share an \code{ID}. Each transcript's parts are joined in
#' order of transcription, so a CDS split across a ribosomal frameshift,
#' such as SARS-CoV-2 ORF1ab, is translated as one sequence. The
#' \code{TxDb} is then built with \code{txdbmaker::makeTxDb()}, which is
#' why that route needs the \pkg{txdbmaker} package.
#'
#' Building a \code{TxDb} directly from a GFF3 file with
#' \code{txdbmaker::makeTxDbFromGFF()} is not equivalent for every
#' annotation. For the NCBI SARS-CoV-2 annotation (NC_045512.2) it drops
#' the ORF1ab gene, whose two CDS parts overlap by one nucleotide, with a
#' warning that its exon ranks could not be inferred; sites in ORF1a and
#' ORF1ab are then left unannotated. A \code{TxDb} gives no record of
#' what was dropped while it was built, so this function cannot detect
#' it. Pass the CDS records as a \code{GRanges} when in doubt.
#'
#' A GFF3 \code{phase} column, when present, is honoured: the first part
#' of each transcript is shortened by its phase, the number of bases
#' before the first complete codon, and the phase of every later part
#' must agree with the length of the parts before it. A later part whose
#' phase disagrees would shift the reading frame of the joined sequence,
#' so the function stops.
#'
#' \strong{Reference.} Sequence names in \code{ref} are cut at the first
#' space, so a FASTA header carrying a description still matches the
#' annotation. The function stops if the variant sites and the coding
#' sequences share no sequence name (see \code{\link{annotateGenes}}),
#' or if a sequence carrying coding sequences is missing from
#' \code{ref}. Sites whose \code{ref} allele differs from the reference
#' base at that position are left unannotated, with a warning, since
#' they were called against a different sequence. The length of each
#' reference sequence is recorded in \code{seqlengths(rowRanges(x))}.
#'
#' \strong{Sites annotated.} Only single-nucleotide substitutions with
#' alleles A, C, G or T are annotated; insertions, deletions and
#' multi-nucleotide substitutions get \code{NA}. Alleles are read on the
#' plus strand and complemented by \code{predictCoding()} for
#' transcripts on the minus strand.
#'
#' \strong{Classification.} \code{CONSEQUENCE} is taken from
#' \code{predictCoding()} unchanged: \code{"synonymous"} when the amino
#' acid is kept, \code{"nonsense"} when a stop codon is gained, and
#' \code{"nonsynonymous"} for every other change, including the loss of
#' a stop codon. A codon holding an ambiguous base is not translated;
#' such a site is left unannotated, with a warning.
#'
#' \code{predictCoding()} applies the initiator rule of
#' \code{Biostrings::translate()} only when the substitution is the
#' first base of a coding sequence: there an alternative start codon
#' (\code{CTG} or \code{TTG}) is reported as \code{M}, so
#' \code{ATG} to \code{CTG} at that one position is synonymous. A
#' change at the second or third base of the first codon is translated
#' with the standard table, giving \code{L} for \code{CTG}.
#'
#' \strong{Sites in several transcripts.} \code{predictCoding()} returns
#' one row per site and transcript, so a site where coding sequences
#' overlap (ORF1a and ORF1ab in SARS-CoV-2) has several. With
#' \code{collapse = "longest"} (the default) the row of the transcript
#' with the longest coding sequence is kept, ties going to the
#' transcript that appears first in \code{cds}, and a message reports
#' how many sites had
#' more than one transcript and at how many of them the consequence
#' differed between transcripts. With \code{collapse = "all"} every row
#' is kept: each column becomes a \code{CharacterList} or
#' \code{IntegerList} with one element per transcript.
#'
#' A nucleotide that one transcript reads twice, the nucleotide shared by
#' the two parts of a CDS split across a -1 ribosomal frameshift (SARS-CoV-2
#' ORF1ab at 13468), gives that transcript two rows, one for each reading.
#' With \code{collapse = "longest"} the first reading, at the lower CDS
#' position, is kept, and a message says how many sites this affected and
#' at how many the two readings differ in consequence;
#' \code{collapse = "all"} keeps both. \code{\link{pNpSWithinHost}} counts
#' such a site by its first reading.
#'
#' @param x A \code{\link{WithinHostExperiment}} with \code{ref} and
#'   \code{alt} columns in \code{mcols(rowRanges(x))}.
#' @param cds The coding sequences: a \code{TxDb}; a \code{GRanges} of
#'   CDS records; or a \code{GRangesList} with one element per
#'   transcript, holding its CDS parts.
#' @param ref A \code{DNAStringSet} with the reference sequence of each
#'   seqlevel carrying coding sequences, for example from
#'   \code{Biostrings::readDNAStringSet()}.
#' @param txCol For a \code{GRanges} \code{cds}: the column identifying
#'   the transcript each CDS record belongs to (default \code{"ID"}).
#' @param geneCol The column of gene labels reported as \code{GENEID}.
#'   For a \code{GRanges} it is a column of the CDS records,
#'   \code{mcols(cds)}, defaulting to \code{"gene"}; for a
#'   \code{GRangesList} it is a column of the list's own metadata, one
#'   value per transcript, defaulting to \code{NULL}. \code{NULL} reports
#'   the transcript name instead of a gene label.
#' @param collapse Either \code{"longest"} (default) or \code{"all"};
#'   see Details.
#' @param ... Not used; supplying anything is an error.
#'
#' @return \code{x}, with these columns added to
#'   \code{mcols(rowRanges(x))} (\code{NA}, or an empty element, for
#'   sites not annotated):
#'   \describe{
#'     \item{TXNAME}{Transcript the row describes.}
#'     \item{GENEID}{Gene of that transcript.}
#'     \item{CONSEQUENCE}{\code{"synonymous"}, \code{"nonsynonymous"} or
#'       \code{"nonsense"}; see Details. A site whose reference codon
#'       holds an ambiguous base is \code{NA} in every column.}
#'     \item{CDSLOC}{Position of the site within the coding sequence.}
#'     \item{PROTEINLOC}{Position of the affected codon in the protein.}
#'     \item{REFCODON, VARCODON}{Reference and variant codon, on the
#'       transcript's strand.}
#'     \item{REFAA, VARAA}{Reference and variant amino acid, \code{"*"}
#'       for a stop codon.}
#'   }
#'
#' @references
#' Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M
#' (2014). VariantAnnotation: a Bioconductor package for exploration
#' and annotation of genetic variants. \emph{Bioinformatics}
#' 30:2076-2078. \doi{10.1093/bioinformatics/btu168}
#'
#' @seealso \code{\link{annotateGenes}} for gene labels without codons;
#'   \code{VariantAnnotation::predictCoding()}.
#'
#' @name annotateCodonChange
#' @aliases annotateCodonChange,WithinHostExperiment,TxDb,DNAStringSet-method
#'   annotateCodonChange,WithinHostExperiment,GRangesList,DNAStringSet-method
#'   annotateCodonChange,WithinHostExperiment,GRanges,DNAStringSet-method
#'   annotateCodonChange,WithinHostExperiment,ANY,character-method
#'   annotateCodonChange,WithinHostExperiment,ANY,XString-method
#'   annotateCodonChange,WithinHostExperiment,character,DNAStringSet-method
#' @examples
#' ## A 15-nucleotide coding sequence, ATG CTG AAA GGG TAA
#' ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
#' cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
#'     strand = "+", ID = "cds-A", gene = "geneA")
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(1, 6, 13), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "G", "T")
#' S4Vectors::mcols(gr)$alt <- c("G", "A", "C")
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' if (requireNamespace("txdbmaker", quietly = TRUE)) {
#'     whe <- annotateCodonChange(whe, cds, ref)
#'     mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))
#'     mc[, c("GENEID", "REFCODON", "VARCODON", "REFAA", "VARAA",
#'            "CONSEQUENCE")]
#' }
NULL

## Group CDS records into one GRangesList element per transcript, keeping
## the phase column (if any) for the frame checks in .cdsListToTxDb() and
## the gene label in mcols() of the list.
.cdsRangesToList <- function(cds, txCol, geneCol) {
    if (length(cds) == 0L)
        stop("'cds' contains no ranges.")
    cols <- colnames(mcols(cds))
    if ("type" %in% cols) {
        other <- setdiff(unique(as.character(mcols(cds)$type)), "CDS")
        if (length(other))
            stop("'cds' holds records of type ",
                 paste(other, collapse = ", "), " as well as CDS. Pass ",
                 "the CDS records only, for example ",
                 "cds[cds$type == \"CDS\"].")
    }
    absent <- setdiff(c(txCol, geneCol), cols)
    if (length(absent))
        stop("Column(s) ", paste(absent, collapse = ", "),
             " not found in mcols(cds). Available columns: ",
             if (length(cols)) paste(cols, collapse = ", ") else "(none)")
    tx <- mcols(cds)[[txCol]]
    if (!is.atomic(tx))
        stop("Column '", txCol, "' must hold one transcript identifier ",
             "per CDS record; it is a ", class(tx)[1L], ".")
    tx <- as.character(tx)
    if (anyNA(tx))
        stop(sum(is.na(tx)), " CDS record(s) have no '", txCol, "'.")
    f <- factor(tx, levels = unique(tx))
    parts <- granges(cds)
    if ("phase" %in% cols)
        mcols(parts)$phase <- mcols(cds)$phase
    grl <- split(parts, f)
    if (!is.null(geneCol)) {
        g <- mcols(cds)[[geneCol]]
        if (!is.atomic(g))
            stop("Column '", geneCol, "' must hold one label per CDS ",
                 "record; it is a ", class(g)[1L], ".")
        per_tx <- split(as.character(g), f)
        gene <- vapply(per_tx, function(z) {
            u <- unique(z[!is.na(z)])
            if (length(u) == 1L) u else NA_character_
        }, character(1))
        bad <- names(gene)[is.na(gene)]
        if (length(bad))
            stop("Transcript(s) ", paste(utils::head(bad, 5L),
                 collapse = ", "), " do not carry exactly one '", geneCol,
                 "' value across their CDS records.")
        mcols(grl)$gene <- unname(gene)
    }
    grl
}

## Build a TxDb from CDS parts grouped by transcript. Parts are ranked in
## order of transcription; overlapping parts (as at a -1 ribosomal
## frameshift) are kept as they are. A GFF3 phase on the first part
## shortens it to the first complete codon.
.cdsListToTxDb <- function(grl, geneCol, ref) {
    if (!requireNamespace("txdbmaker", quietly = TRUE))
        stop("Annotating from CDS ranges builds a TxDb with the txdbmaker ",
             "package, which is not installed. Install it with ",
             "BiocManager::install(\"txdbmaker\"), or pass a TxDb as ",
             "'cds'.", call. = FALSE)
    if (length(grl) == 0L)
        stop("'cds' contains no transcripts.")
    tx_names <- names(grl)
    if (is.null(tx_names) || anyNA(tx_names) || any(!nzchar(tx_names)) ||
        anyDuplicated(tx_names))
        stop("'cds' must have unique names, one per transcript.")
    n_parts <- lengths(grl)
    if (any(n_parts == 0L))
        stop("Transcript(s) with no CDS range: ",
             paste(tx_names[n_parts == 0L], collapse = ", "))
    gene <- if (is.null(geneCol)) {
        tx_names
    } else {
        if (!geneCol %in% colnames(mcols(grl)))
            stop("Column '", geneCol, "' not found in mcols(cds). ",
                 "Available columns: ",
                 if (ncol(mcols(grl))) paste(colnames(mcols(grl)),
                                             collapse = ", ")
                 else "(none)", ". Use geneCol = NULL to report ",
                 "transcript names.")
        as.character(mcols(grl)[[geneCol]])
    }
    if (anyNA(gene))
        stop("Transcript(s) ", paste(tx_names[is.na(gene)], collapse = ", "),
             " have no '", geneCol, "'.")

    parts <- unlist(grl, use.names = FALSE)
    tx_index <- rep(seq_along(grl), n_parts)
    chrom <- as.character(seqnames(parts))
    strand_chr <- as.character(strand(parts))
    if (any(strand_chr == "*"))
        stop("Every CDS range needs a strand; transcript(s) ",
             paste(unique(tx_names[tx_index[strand_chr == "*"]]),
                   collapse = ", "), " have ranges without one.")
    mixed <- vapply(split(paste(chrom, strand_chr), tx_index),
                    function(z) length(unique(z)) > 1L, logical(1))
    if (any(mixed))
        stop("Transcript(s) ", paste(tx_names[mixed], collapse = ", "),
             " have CDS ranges on more than one sequence or strand.")
    no_ref <- setdiff(unique(chrom), names(ref))
    if (length(no_ref))
        .stopNoRef(no_ref, ref)

    o <- order(tx_index,
               ifelse(strand_chr == "-", -start(parts), start(parts)))
    parts <- parts[o]
    tx_index <- tx_index[o]
    chrom <- chrom[o]
    strand_chr <- strand_chr[o]
    rank <- stats::ave(seq_along(tx_index), tx_index, FUN = seq_along)

    if ("phase" %in% colnames(mcols(parts))) {
        ## a GFF3 phase is 0, 1, 2 or missing; rtracklayer gives an
        ## integer column, other readers may give text
        raw_phase <- as.character(mcols(parts)$phase)
        known <- !is.na(raw_phase) & raw_phase != "."
        if (any(known & !raw_phase %in% c("0", "1", "2")))
            stop("GFF3 phase must be 0, 1 or 2.")
        phase <- ifelse(known, match(raw_phase, c("0", "1", "2")) - 1L,
                        NA_integer_)
        first_phase <- rank == 1L & !is.na(phase) & phase > 0L
        if (any(first_phase)) {
            k <- which(first_phase)
            if (any(width(parts)[k] <= phase[k]))
                stop("Transcript(s) ",
                     paste(unique(tx_names[tx_index[k[width(parts)[k] <=
                                                      phase[k]]]]),
                           collapse = ", "),
                     " have a first CDS part no longer than its phase.")
            plus <- k[strand_chr[k] == "+"]
            minus <- k[strand_chr[k] == "-"]
            start(parts)[plus] <- start(parts)[plus] + phase[plus]
            end(parts)[minus] <- end(parts)[minus] - phase[minus]
        }
        w <- width(parts)
        before <- stats::ave(w, tx_index,
                             FUN = function(v) c(0L, cumsum(v)[-length(v)]))
        expected <- (3L - before %% 3L) %% 3L
        later_bad <- rank > 1L & !is.na(phase) & phase != expected
        if (any(later_bad))
            stop("For transcript(s) ",
                 paste(unique(tx_names[tx_index[later_bad]]),
                       collapse = ", "),
                 ", the phase of a later CDS part does not follow from the ",
                 "length of the parts before it, so joining the parts would ",
                 "shift the reading frame.")
    }

    first <- !duplicated(tx_index)
    transcripts_df <- data.frame(
        tx_id = seq_along(grl), tx_name = tx_names,
        tx_chrom = chrom[first], tx_strand = strand_chr[first],
        tx_start = as.integer(tapply(start(parts), tx_index, min)),
        tx_end = as.integer(tapply(end(parts), tx_index, max)),
        stringsAsFactors = FALSE)
    splicings <- data.frame(
        tx_id = tx_index, exon_rank = rank,
        exon_start = start(parts), exon_end = end(parts),
        cds_start = start(parts), cds_end = end(parts))
    genes <- data.frame(tx_id = seq_along(grl), gene_id = gene,
                        stringsAsFactors = FALSE)
    chrominfo <- data.frame(chrom = names(ref), length = width(ref),
                            is_circular = FALSE, stringsAsFactors = FALSE)
    withCallingHandlers(
        txdbmaker::makeTxDb(transcripts_df, splicings, genes,
                            chrominfo = chrominfo),
        warning = function(w) {
            if (grepl("genome version information is not available",
                      conditionMessage(w), fixed = TRUE))
                invokeRestart("muffleWarning")
        })
}

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "TxDb", "DNAStringSet"),
    function(x, cds, ref, collapse = c("longest", "all"), ...) {
        .stop_on_unused_dots("annotateCodonChange", ...)
        collapse <- match.arg(collapse)
        ref <- .prepareReference(ref)
        rr <- rowRanges(x)
        mc <- mcols(rr)
        if (!all(c("ref", "alt") %in% colnames(mc)))
            stop("mcols(rowRanges(x)) must contain 'ref' and 'alt' columns.")

        cds_by_tx <- cdsBy(cds, by = "tx")
        if (length(cds_by_tx) == 0L)
            stop("'cds' contains no coding sequences.")
        cds_parts <- unlist(cds_by_tx, use.names = FALSE)
        .checkSeqlevels(rr, cds_parts, "cds")
        no_ref <- setdiff(seqlevelsInUse(cds_parts), names(ref))
        if (length(no_ref))
            .stopNoRef(no_ref, ref)
        n <- length(rr)
        ref_allele <- toupper(as.character(mc$ref))
        alt_allele <- toupper(as.character(mc$alt))
        chrom <- as.character(seqnames(rr))
        pos <- start(rr)
        snv <- width(rr) == 1L & !is.na(ref_allele) & !is.na(alt_allele) &
            grepl("^[ACGT]$", ref_allele) & grepl("^[ACGT]$", alt_allele)
        len <- width(ref)[match(chrom, names(ref))]
        beyond <- snv & !is.na(len) & pos > len
        if (any(beyond))
            warning(sum(beyond), " site(s) lie beyond the end of their ",
                    "sequence in 'ref' and are left unannotated.",
                    call. = FALSE)
        on_ref <- snv & !is.na(len) & !beyond
        base <- rep(NA_character_, n)
        if (any(on_ref))
            base[on_ref] <- as.character(subseq(
                ref[chrom[on_ref]], start = pos[on_ref], width = 1L))
        mismatch <- on_ref & base != ref_allele
        if (any(mismatch)) {
            k <- which(mismatch)[1L]
            warning(sum(mismatch), " site(s) carry a 'ref' allele that ",
                    "differs from the base in 'ref' at that position ",
                    "(first: ", chrom[k], ":", pos[k], ", ", ref_allele[k],
                    " against ", base[k], ") and are left unannotated. ",
                    "Check that 'ref' is the sequence the variants were ",
                    "called against.", call. = FALSE)
        }
        rr <- .fillSeqlengths(rr, ref)
        coding_seq <- chrom %in% seqlevelsInUse(cds_parts)
        use <- which(on_ref & !mismatch & coding_seq)

        cols <- list(TXNAME = character(), GENEID = character(),
                     CONSEQUENCE = character(), CDSLOC = integer(),
                     PROTEINLOC = integer(), REFCODON = character(),
                     VARCODON = character(), REFAA = character(),
                     VARAA = character())
        site <- integer()
        if (length(use)) {
            query <- GRanges(chrom[use], IRanges(pos[use], width = 1L))
            pc <- withCallingHandlers(
                predictCoding(query, cds, seqSource = ref,
                              varAllele = DNAStringSet(alt_allele[use])),
                warning = function(w) {
                    if (grepl("were not translated", conditionMessage(w),
                              fixed = TRUE))
                        invokeRestart("muffleWarning")
                })
            ## A codon holding N, ., + or - comes back as "not translated"
            ## with empty amino acids. Such a row carries no consequence,
            ## so the site is left unannotated rather than counted.
            if (length(pc)) {
                blank <- as.character(pc$CONSEQUENCE) == "not translated" |
                    !nzchar(as.character(pc$REFAA)) |
                    !nzchar(as.character(pc$VARAA))
                if (any(blank)) {
                    warning(length(unique(pc$QUERYID[blank])), " site(s) ",
                            "have a reference codon with an ambiguous base ",
                            "and are left unannotated.", call. = FALSE)
                    pc <- pc[!blank]
                }
            }
            if (length(pc)) {
                site <- use[pc$QUERYID]
                txid <- as.character(pc$TXID)
                tx_tab <- transcripts(cds, columns = c("tx_id", "tx_name"))
                tx_name <- as.character(tx_tab$tx_name)
                names(tx_name) <- as.character(tx_tab$tx_id)
                txname <- tx_name[txid]
                txname[is.na(txname)] <- txid[is.na(txname)]
                chr <- function(v) unname(as.character(v))
                cols <- list(
                    TXNAME = unname(txname),
                    GENEID = chr(pc$GENEID),
                    CONSEQUENCE = chr(pc$CONSEQUENCE),
                    CDSLOC = unname(start(pc$CDSLOC)),
                    PROTEINLOC = unname(as.integer(min(pc$PROTEINLOC))),
                    REFCODON = chr(pc$REFCODON),
                    VARCODON = chr(pc$VARCODON),
                    REFAA = chr(pc$REFAA),
                    VARAA = chr(pc$VARAA))
                tx_len <- sum(width(cds_by_tx))[txid]
                ## the first row of a site is the one "longest" keeps: the
                ## longest transcript and, where that transcript reads the
                ## nucleotide twice, its first reading
                o <- order(site, -tx_len, as.integer(txid), cols$CDSLOC)
                site <- site[o]
                cols <- lapply(cols, `[`, o)
            }
        }

        if (collapse == "longest") {
            ## one row per site and transcript: its first reading
            tx_key <- paste(site, cols$TXNAME, sep = "\r")
            first_tx <- !duplicated(tx_key)
            multi <- tabulate(site[first_tx], nbins = n) > 1L
            if (any(multi)) {
                by_site <- split(cols$CONSEQUENCE[first_tx],
                                 factor(site[first_tx],
                                        levels = seq_len(n)))[multi]
                differ <- vapply(by_site, function(z) length(unique(z)) > 1L,
                                 logical(1))
                message(sum(multi), " site(s) fall in more than one ",
                        "transcript; the longest was kept. The consequence ",
                        "differs between transcripts at ", sum(differ),
                        " of them. Use collapse = \"all\" to keep every ",
                        "transcript.")
            }
            keep <- !duplicated(site)
            ## a nucleotide the kept transcript reads twice
            kept_rows <- tx_key %in% tx_key[keep]
            readings <- split(cols$CONSEQUENCE[kept_rows],
                              factor(site[kept_rows], levels = seq_len(n)))
            twice <- lengths(readings) > 1L
            if (any(twice)) {
                differ <- vapply(readings[twice],
                                 function(z) length(unique(z)) > 1L,
                                 logical(1))
                message(sum(twice), " site(s) lie on a nucleotide that the ",
                        "kept transcript reads twice, as at a -1 ribosomal ",
                        "frameshift; the first reading was kept. The ",
                        "consequence differs between the two readings at ",
                        sum(differ), " of them. Use collapse = \"all\" to ",
                        "keep both.")
            }
            out <- lapply(cols, function(v) {
                filled <- rep(v[0L][NA_integer_], n)
                filled[site[keep]] <- v[keep]
                filled
            })
        } else {
            f <- factor(site, levels = seq_len(n))
            out <- lapply(cols, function(v) {
                l <- if (is.integer(v)) IntegerList(split(v, f))
                     else CharacterList(split(v, f))
                names(l) <- NULL
                l
            })
        }
        for (nm in names(out)) mcols(rr)[[nm]] <- out[[nm]]
        SummarizedExperiment::rowRanges(x) <- rr
        x
    })

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "GRangesList", "DNAStringSet"),
    function(x, cds, ref, geneCol = NULL, collapse = c("longest", "all"),
             ...) {
        ref <- .prepareReference(ref)
        txdb <- .cdsListToTxDb(cds, geneCol, ref)
        annotateCodonChange(x, txdb, ref, collapse = collapse, ...)
    })

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "GRanges", "DNAStringSet"),
    function(x, cds, ref, txCol = "ID", geneCol = "gene",
             collapse = c("longest", "all"), ...) {
        grl <- .cdsRangesToList(cds, txCol, geneCol)
        annotateCodonChange(x, grl, ref,
                            geneCol = if (is.null(geneCol)) NULL else "gene",
                            collapse = collapse, ...)
    })

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "ANY", "character"),
    function(x, cds, ref, ...) .refPathError(ref))

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "ANY", "XString"),
    function(x, cds, ref, ...) .xstringRefError(ref))

#' @rdname annotateCodonChange
#' @export
setMethod("annotateCodonChange",
    c("WithinHostExperiment", "character", "DNAStringSet"),
    function(x, cds, ref, ...) .cdsPathError(cds))
