# R/annotate-genes.R
# Label iSNV sites with the genes they fall in

#' @include methods.R
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb seqlevelsInUse
#' @importFrom IRanges CharacterList width
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
#' @importFrom SummarizedExperiment rowRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicFeatures TxDb
NULL

## Stop when the variant sites and a set of ranges share no sequence
## name, and warn about variant sequences the ranges do not cover.
## Without this check a naming mismatch (for example NC_045512.2 against
## MN908947.3 for the same SARS-CoV-2 reference) gives an annotation
## column that is NA everywhere, or a QC step that flags nothing, and no
## error. 'effect' says what the ranges would have done to a site.
.checkSeqlevels <- function(variants, annotation, what,
                            effect = "annotated") {
    v <- seqlevelsInUse(variants)
    a <- seqlevelsInUse(annotation)
    if (length(v) == 0L) return(invisible(TRUE))
    if (length(intersect(v, a)) == 0L) {
        stop("The variant sites and '", what, "' share no sequence ",
             "names, so no site could be ", effect, ".\n",
             "  variant seqlevels: ", paste(v, collapse = ", "), "\n",
             "  '", what, "' seqlevels: ",
             if (length(a)) paste(a, collapse = ", ") else "(none)", "\n",
             "Rename one of them to match, for example with ",
             "GenomeInfoDb::renameSeqlevels().", call. = FALSE)
    }
    uncovered <- setdiff(v, a)
    if (length(uncovered)) {
        warning("'", what, "' has no ranges on variant seqlevel(s) ",
                paste(uncovered, collapse = ", "),
                "; sites there are not ", effect, ".", call. = FALSE)
    }
    invisible(TRUE)
}

#' Label variant sites with the genes they overlap
#'
#' Adds a \code{GENE} column to \code{mcols(rowRanges(x))} giving, for
#' each variant site, the label of the gene range it falls in. Genes are
#' supplied as a \code{GRanges}, which is how
#' \code{rtracklayer::import()} returns a GFF3 or GTF file, or as a
#' \code{TxDb}, from which \code{GenomicFeatures::genes()} extracts them.
#' The overlap itself is \code{GenomicRanges::findOverlaps()} with strand
#' ignored, since variant sites carry no strand.
#'
#' Before overlapping, the sequence names in use on both sides are
#' compared. If they share none, the function stops and prints both sets:
#' the same reference is often named differently by different sources
#' (for SARS-CoV-2, \code{NC_045512.2} in RefSeq and \code{MN908947.3} in
#' GenBank), and an overlap across mismatched names would otherwise
#' return a column that is \code{NA} everywhere. If they share some but
#' not all, a warning names the variant sequences that have no gene range.
#'
#' A site can fall in more than one gene where genes overlap.
#' \code{collapse = "longest"} keeps the label of the longest overlapping
#' range (the first in \code{genes} on a tie) and reports with a message
#' how many sites carried more than one distinct label;
#' \code{collapse = "all"} keeps every distinct label in a
#' \code{CharacterList}, one element per site, in the order of
#' \code{genes}. Several ranges with the same label, such as the two CDS
#' parts of a gene split by a ribosomal frameshift, count as one gene.
#'
#' Ranges whose label is missing are dropped with a warning; if none is
#' left, the function stops. For a \code{TxDb},
#' \code{GenomicFeatures::genes()} leaves out, with a message, genes whose
#' transcripts lie on more than one strand or sequence, so sites in them
#' are not labelled.
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param genes A \code{GRanges} of gene ranges, or a \code{TxDb}.
#' @param nameCol Character scalar. The column of \code{mcols(genes)}
#'   holding the label to copy. Defaults to \code{"Name"}, the GFF3
#'   attribute that \code{rtracklayer::import()} returns as a column; for
#'   a \code{TxDb} it defaults to \code{"gene_id"}. An unknown column is an
#'   error that lists the columns available.
#' @param collapse Either \code{"longest"} (default) or \code{"all"}; see
#'   Details.
#' @param ... Not used; supplying anything is an error.
#'
#' @return \code{x} with a \code{GENE} column in
#'   \code{mcols(rowRanges(x))}: a character vector for
#'   \code{collapse = "longest"}, \code{NA} for sites in no gene, or a
#'   \code{CharacterList} for \code{collapse = "all"}, with an empty
#'   element for sites in no gene.
#'
#' @seealso \code{\link{annotateCodonChange}} for the codon and amino
#'   acid consequence of each substitution.
#'
#' @name annotateGenes
#' @aliases annotateGenes,WithinHostExperiment,GRanges-method
#'   annotateGenes,WithinHostExperiment,TxDb-method
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#'
#' ## Gene ranges built directly
#' genes <- GenomicRanges::GRanges("test_segment",
#'     IRanges::IRanges(c(1, 501), c(500, 1000)), strand = "+",
#'     Name = c("geneA", "geneB"))
#' whe <- annotateGenes(whe, genes)
#' table(S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$GENE,
#'     useNA = "ifany")
#'
#' ## The same genes read from a GFF3 file
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'     gff <- tempfile(fileext = ".gff3")
#'     writeLines(c("##gff-version 3",
#'         "test_segment\t.\tgene\t1\t500\t.\t+\t.\tID=gA;Name=geneA",
#'         "test_segment\t.\tgene\t501\t1000\t.\t+\t.\tID=gB;Name=geneB"),
#'         gff)
#'     features <- rtracklayer::import(gff)
#'     whe <- annotateGenes(whe, features[features$type == "gene"])
#' }
#'
#' ## Sequence names that do not match are an error, not a column of NA
#' renamed <- GenomeInfoDb::renameSeqlevels(genes,
#'     c(test_segment = "segment_1"))
#' try(annotateGenes(whe, renamed))
NULL

#' @rdname annotateGenes
#' @export
setMethod("annotateGenes", c("WithinHostExperiment", "GRanges"),
    function(x, genes, nameCol = "Name", collapse = c("longest", "all"),
             ...) {
        .stop_on_unused_dots("annotateGenes", ...)
        collapse <- match.arg(collapse)
        if (!is.character(nameCol) || length(nameCol) != 1L ||
            is.na(nameCol))
            stop("'nameCol' must be a single column name.")
        if (length(genes) == 0L)
            stop("'genes' contains no ranges.")
        available <- colnames(mcols(genes))
        if (!nameCol %in% available)
            stop("Column '", nameCol, "' not found in mcols(genes). ",
                 "Available columns: ",
                 if (length(available)) paste(available, collapse = ", ")
                 else "(none)")
        labels <- mcols(genes)[[nameCol]]
        if (is(labels, "Rle"))
            labels <- S4Vectors::decode(labels)
        if (!is.atomic(labels))
            stop("Column '", nameCol, "' must hold one label per range; ",
                 "it is a ", class(labels)[1L], ".")
        labels <- as.character(labels)
        unlabelled <- is.na(labels)
        if (any(unlabelled)) {
            if (all(unlabelled))
                stop("No range in 'genes' has a '", nameCol, "' label.",
                     call. = FALSE)
            warning(sum(unlabelled), " range(s) in 'genes' have no '",
                    nameCol, "' and are ignored.", call. = FALSE)
            genes <- genes[!unlabelled]
            labels <- labels[!unlabelled]
        }

        rr <- rowRanges(x)
        .checkSeqlevels(rr, genes, "genes")

        hits <- findOverlaps(rr, genes, ignore.strand = TRUE)
        q <- queryHits(hits)
        s <- subjectHits(hits)
        ## sites carrying more than one distinct label
        distinct <- !duplicated(data.frame(q, labels[s]))
        n_multi <- sum(tabulate(q[distinct], nbins = length(rr)) > 1L)
        if (collapse == "all") {
            o <- order(q, s)
            gene <- unique(CharacterList(split(labels[s[o]],
                factor(q[o], levels = seq_along(rr)))))
            names(gene) <- NULL
        } else {
            o <- order(q, -width(genes)[s], s)
            q <- q[o]
            s <- s[o]
            keep <- !duplicated(q)
            if (n_multi > 0L)
                message(n_multi, " site(s) fall in more than one gene; ",
                        "the longest was kept. Use collapse = \"all\" ",
                        "to keep every label.")
            gene <- rep(NA_character_, length(rr))
            gene[q[keep]] <- labels[s[keep]]
        }
        mcols(rr)$GENE <- gene
        SummarizedExperiment::rowRanges(x) <- rr
        x
    })

#' @rdname annotateGenes
#' @export
setMethod("annotateGenes", c("WithinHostExperiment", "TxDb"),
    function(x, genes, nameCol = "gene_id", collapse = c("longest", "all"),
             ...) {
        annotateGenes(x, GenomicFeatures::genes(genes), nameCol = nameCol,
                      collapse = collapse, ...)
    })
