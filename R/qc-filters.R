# R/qc-filters.R
# Quality control: flag-based filtering (never deletes rows)

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames colData
#'   rowRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors DataFrame queryHits
NULL

# ============================================================
# flagISNV -- main QC function
# ============================================================

#' Flag iSNV sites based on quality criteria
#'
#' Applies quality control filters to a \code{\link{WithinHostExperiment}}.
#' Instead of removing variants, this function sets the corresponding
#' entries in the \code{qcPass} assay to \code{FALSE}. Each filtering
#' step is recorded in the \code{qcLog} for full reproducibility.
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @param filter An \code{\link{ISNVFilter}} object (default:
#'   \code{ISNVFilter()}).
#' @param ... Must be empty. Any further argument is an error, so that a
#'   misspelled argument name cannot be silently ignored.
#'
#' @return A \code{\link{WithinHostExperiment}} with updated
#'   \code{qcPass} assay and \code{qcLog}.
#'
#' @details
#' Filtering steps are applied in order:
#' \enumerate{
#'   \item Minimum read depth (\code{minDepth})
#'   \item Minimum alternative allele frequency (\code{minFreq})
#'   \item Maximum alternative allele frequency (\code{maxFreq})
#'   \item Minimum alternative allele read count (\code{minAltReads})
#'   \item Strand bias (\code{maxStrandBias}) -- only if
#'     \code{fwdAltCount} and \code{revAltCount} assays exist
#'   \item Replicate concordance (\code{replicateConc}) -- only if the
#'     filter holds a \code{\link{ConcordanceRule}}
#' }
#'
#' The filter is checked with \code{validObject()} first, so a threshold
#' or a rule changed with \code{@<-} into one its class refuses is an
#' error rather than a step applied.
#'
#' Base and mapping quality are not filtered here: the object holds no
#' per-site quality, so those thresholds belong in the variant caller.
#'
#' The replicate step compares every pair of samples that share a value
#' of the \code{replicate_group} column of \code{colData} and flags a site
#' that fails the rule in every sample of the group. It compares the calls
#' as made, before the thresholds above, so it flags what
#' \code{\link{flagReplicateDiscordance}} flags with the same settings,
#' and it records the same two \code{qcLog} rows. A rule with no
#' \code{replicate_group} column, or with no group of two or more samples,
#' is an error rather than a step that checks nothing. Samples outside
#' any such group are named in a warning and are filtered by the other
#' thresholds only.
#'
#' @seealso \code{\link{ISNVFilter}} for filter configuration,
#'   \code{\link{ConcordanceRule}} for the replicate rule,
#'   \code{\link{passedISNV}} for extracting passed variants,
#'   \code{\link{qcSummary}} for QC report,
#'   \code{\link{flagReplicateDiscordance}} for replicate QC.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
#' qcLog(whe)
#' passedISNV(whe)
#' @rdname flagISNV
#' @aliases flagISNV,WithinHostExperiment-method
setMethod("flagISNV", "WithinHostExperiment",
    function(x, filter = ISNVFilter(), ...) {
    whe <- x
    .stop_on_unused_dots("flagISNV", ...)
    if (!is(filter, "ISNVFilter"))
        stop("'filter' must be an ISNVFilter object. ",
             "Create one with ISNVFilter().")
    ## a slot changed with @<- after construction, or a filter saved by an
    ## earlier version, is checked here rather than applied
    validObject(filter)

    if (!"qcPass" %in% assayNames(whe)) {
        stop("WithinHostExperiment must have a 'qcPass' assay")
    }

    ## A replicate rule with nothing to compare must not pass unnoticed
    rule <- slot(filter, "replicateConc")
    if (!is.null(rule))
        .checkReplicateGroups(colData(whe), colnames(whe))

    qc <- assay(whe, "qcPass")

    ## ---- Step 1: Depth filter ----
    if ("totalDepth" %in% assayNames(whe)) {
        depth_mat <- assay(whe, "totalDepth")
        min_depth <- slot(filter, "minDepth")
        fail <- !is.na(depth_mat) & depth_mat < min_depth & qc
        n_flag <- sum(fail, na.rm = TRUE)
        qc[fail] <- FALSE
        whe <- .addQcStep(whe, "depth_filter", "minDepth",
                          min_depth, n_flag,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Step 2: Min frequency filter ----
    if ("altFreq" %in% assayNames(whe)) {
        freq_mat <- assay(whe, "altFreq")
        min_freq <- slot(filter, "minFreq")
        fail <- !is.na(freq_mat) & freq_mat < min_freq & qc
        n_flag <- sum(fail, na.rm = TRUE)
        qc[fail] <- FALSE
        whe <- .addQcStep(whe, "min_freq_filter", "minFreq",
                          min_freq, n_flag,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Step 3: Max frequency filter ----
    if ("altFreq" %in% assayNames(whe)) {
        freq_mat <- assay(whe, "altFreq")
        max_freq <- slot(filter, "maxFreq")
        fail <- !is.na(freq_mat) & freq_mat > max_freq & qc
        n_flag <- sum(fail, na.rm = TRUE)
        qc[fail] <- FALSE
        whe <- .addQcStep(whe, "max_freq_filter", "maxFreq",
                          max_freq, n_flag,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Step 4: Min alt reads filter ----
    if ("altCount" %in% assayNames(whe)) {
        alt_mat <- assay(whe, "altCount")
        min_alt <- slot(filter, "minAltReads")
        fail <- !is.na(alt_mat) & alt_mat < min_alt & qc
        n_flag <- sum(fail, na.rm = TRUE)
        qc[fail] <- FALSE
        whe <- .addQcStep(whe, "alt_reads_filter", "minAltReads",
                          min_alt, n_flag,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Step 5: Strand bias filter (only if strand assays exist) ----
    ## Uses fold-difference max(fwd/rev, rev/fwd) as a simple strand
    ## bias metric. For amplicon data with inherent strand coverage
    ## imbalance, flagStrandBiasSCI() normalises by the global strand
    ## coverage ratio instead; it needs fwdTotalCount/revTotalCount
    ## assays. Mostefai et al. (2024) NAR Genomics & Bioinformatics
    ## 6:lqae145 use a per-strand binomial likelihood.
    if (all(c("fwdAltCount", "revAltCount") %in% assayNames(whe))) {
        fwd <- assay(whe, "fwdAltCount")
        rev_counts <- assay(whe, "revAltCount")
        ## Strand bias = max(fwd/rev, rev/fwd); avoid div by zero
        fwd_safe <- pmax(fwd, 1L)
        rev_safe <- pmax(rev_counts, 1L)
        sb <- pmax(fwd_safe / rev_safe, rev_safe / fwd_safe)
        max_sb <- slot(filter, "maxStrandBias")
        fail <- !is.na(sb) & sb > max_sb & qc
        n_flag <- sum(fail, na.rm = TRUE)
        qc[fail] <- FALSE
        whe <- .addQcStep(whe, "strand_bias_filter", "maxStrandBias",
                          max_sb, n_flag,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Step 6: Replicate concordance (only with a ConcordanceRule) ----
    if (!is.null(rule)) {
        tol <- slot(rule, "freqTolerance")
        both <- slot(rule, "requireBothDetected")
        res <- .flagDiscordantReplicates(colData(whe)$replicate_group,
                                         assay(whe, "altFreq"), qc,
                                         tol, both)
        qc <- res$qc
        whe <- .addQcStep(whe, "replicate_discordance", "freqTolerance",
                          tol, res$n_flagged, sum(qc, na.rm = TRUE))
        whe <- .addQcStep(whe, "replicate_discordance",
                          "requireBothDetected", both, 0L,
                          sum(qc, na.rm = TRUE))
    }

    ## ---- Write back qcPass and record filter ----
    assay(whe, "qcPass") <- qc
    qcFilters(whe) <- c(qcFilters(whe), list(filter))
    whe
})


# ============================================================
# qcSummary
# ============================================================

#' Generate a QC summary report
#'
#' Produces per-sample summary statistics from a
#' \code{\link{WithinHostExperiment}}.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#'
#' @return A \code{DataFrame} with one row per sample and columns:
#'   \describe{
#'     \item{sample_id}{Sample identifier.}
#'     \item{total_sites}{The sample's calls: sites where its frequency is
#'       above 0. A recorded 0 is a site assessed without the allele and
#'       is not counted.}
#'     \item{passed_sites, flagged_sites}{The calls whose \code{qcPass}
#'       is \code{TRUE} and \code{FALSE}.}
#'     \item{pass_rate}{\code{passed_sites / total_sites}, \code{NA}
#'       without calls.}
#'     \item{median_depth, median_freq}{The median read depth and the
#'       median frequency of the calls.}
#'   }
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
#' qcSummary(whe)
qcSummary <- function(whe) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    qc_mat <- if ("qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        matrix(TRUE, nrow = nrow(whe), ncol = ncol(whe))
    }

    freq_mat <- if ("altFreq" %in% assayNames(whe)) {
        assay(whe, "altFreq")
    } else {
        NULL
    }

    depth_mat <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")
    } else {
        NULL
    }

    sample_ids <- colnames(qc_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(qc_mat)))
    }

    n_samples <- ncol(qc_mat)

    ## A call is a frequency above 0; a recorded 0 is a site assessed
    ## without the allele and is not counted
    is_call <- if (is.null(freq_mat)) {
        matrix(TRUE, nrow = nrow(qc_mat), ncol = n_samples)
    } else {
        .carries(freq_mat)
    }

    DataFrame(
        sample_id = sample_ids,
        total_sites = vapply(seq_len(n_samples), function(j) {
            sum(is_call[, j])
        }, integer(1)),
        passed_sites = vapply(seq_len(n_samples), function(j) {
            sum(qc_mat[, j] & is_call[, j], na.rm = TRUE)
        }, integer(1)),
        flagged_sites = vapply(seq_len(n_samples), function(j) {
            sum(!qc_mat[, j] & is_call[, j], na.rm = TRUE)
        }, integer(1)),
        pass_rate = vapply(seq_len(n_samples), function(j) {
            total <- sum(is_call[, j])
            if (total == 0L) return(NA_real_)
            sum(qc_mat[, j] & is_call[, j], na.rm = TRUE) / total
        }, numeric(1)),
        median_depth = vapply(seq_len(n_samples), function(j) {
            if (is.null(depth_mat)) return(NA_real_)
            stats::median(depth_mat[is_call[, j], j], na.rm = TRUE)
        }, numeric(1)),
        median_freq = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) return(NA_real_)
            stats::median(freq_mat[is_call[, j], j], na.rm = TRUE)
        }, numeric(1))
    )
}


# ============================================================
# flagPrimerSites -- amplicon primer region QC
# ============================================================

#' Flag iSNVs in primer binding regions
#'
#' For amplicon-based sequencing data, marks variants that fall within
#' primer binding sites in the \code{qcPass} assay.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param primers \code{GRanges}. Primer binding site coordinates.
#'   Can be loaded from a BED file using
#'   \code{rtracklayer::import.bed()}. Their sequence names must match
#'   those of the sites: when the two share none, the function stops and
#'   prints both, rather than report that no site lies in a primer.
#' @param action Character. \code{"flag"} (default) sets qcPass to
#'   FALSE for affected sites, preserving all data for audit.
#'   \code{"remove"} is provided as an explicit convenience for
#'   users who require destructive filtering after inspection.
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' primers <- GRanges("test_segment", IRanges::IRanges(
#'     start = c(40, 880), width = c(25, 30)))
#' whe <- flagPrimerSites(whe, primers)
#' qcLog(whe)
flagPrimerSites <- function(whe, primers, action = c("flag", "remove")) {
    action <- match.arg(action)
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!is(primers, "GRanges"))
        stop("'primers' must be a GRanges object.")

    ## A naming mismatch between primers and sites would flag nothing
    ## and look like a clean result; an empty primer set is not one
    if (length(primers))
        .checkSeqlevels(rowRanges(whe), primers, "primers",
                        effect = "flagged")

    ## Find variant sites overlapping primers
    hits <- findOverlaps(rowRanges(whe), primers)
    affected_rows <- unique(queryHits(hits))

    if (length(affected_rows) == 0L) {
        message("No variant sites overlap primer regions")
        whe <- .addQcStep(whe, "primer_filter", "n_primer_regions",
                          length(primers), 0L,
                          sum(assay(whe, "qcPass"), na.rm = TRUE))
        return(whe)
    }

    if (action == "flag") {
        qc <- assay(whe, "qcPass")
        n_before <- sum(qc, na.rm = TRUE)
        qc[affected_rows, ] <- FALSE
        assay(whe, "qcPass") <- qc
        n_after <- sum(qc, na.rm = TRUE)
        whe <- .addQcStep(whe, "primer_filter", "n_primer_regions",
                          length(primers),
                          n_before - n_after, n_after)
    } else {
        n_before <- nrow(whe)
        whe <- whe[-affected_rows, ]
        whe <- .addQcStep(whe, "primer_filter_remove",
                          "n_primer_regions",
                          length(primers),
                          n_before - nrow(whe), nrow(whe))
    }
    whe
}
