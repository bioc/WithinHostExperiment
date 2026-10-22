# R/qc-filters.R
# Quality control: flag-based filtering (never deletes rows)

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames rowRanges
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
#' }
#'
#' Base and mapping quality are not filtered here: the object holds no
#' per-site quality, so those thresholds belong in the variant caller.
#' The \code{replicateConc} slot of the filter is recorded for
#' provenance but not applied here; \code{TRUE} raises a warning. Use
#' \code{\link{flagReplicateDiscordance}} for replicate concordance.
#'
#' @seealso \code{\link{ISNVFilter}} for filter configuration,
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

    if (!"qcPass" %in% assayNames(whe)) {
        stop("WithinHostExperiment must have a 'qcPass' assay")
    }

    ## A threshold this function cannot apply must not pass unnoticed
    if (isTRUE(slot(filter, "replicateConc"))) {
        warning("replicateConc is recorded but not applied by ",
                "flagISNV(); use flagReplicateDiscordance().",
                call. = FALSE)
    }

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
#' @return A \code{DataFrame} with columns: sample_id, total_sites,
#'   passed_sites, flagged_sites, pass_rate, median_depth, median_freq.
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

    DataFrame(
        sample_id = sample_ids,
        total_sites = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) return(nrow(qc_mat))
            sum(!is.na(freq_mat[, j]))
        }, integer(1)),
        passed_sites = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) return(sum(qc_mat[, j], na.rm = TRUE))
            sum(qc_mat[, j] & !is.na(freq_mat[, j]), na.rm = TRUE)
        }, integer(1)),
        flagged_sites = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) {
                return(sum(!qc_mat[, j], na.rm = TRUE))
            }
            sum(!qc_mat[, j] & !is.na(freq_mat[, j]), na.rm = TRUE)
        }, integer(1)),
        pass_rate = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) {
                total <- nrow(qc_mat)
            } else {
                total <- sum(!is.na(freq_mat[, j]))
            }
            if (total == 0L) return(NA_real_)
            if (is.null(freq_mat)) {
                sum(qc_mat[, j], na.rm = TRUE) / total
            } else {
                sum(qc_mat[, j] & !is.na(freq_mat[, j]),
                    na.rm = TRUE) / total
            }
        }, numeric(1)),
        median_depth = vapply(seq_len(n_samples), function(j) {
            if (is.null(depth_mat)) return(NA_real_)
            stats::median(depth_mat[, j], na.rm = TRUE)
        }, numeric(1)),
        median_freq = vapply(seq_len(n_samples), function(j) {
            if (is.null(freq_mat)) return(NA_real_)
            stats::median(freq_mat[, j], na.rm = TRUE)
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
#'   \code{rtracklayer::import.bed()}.
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
