# R/AllClasses.R
# S4 class definitions for WithinHostExperiment
# Collate order: loaded first

#' @importFrom methods setClass setValidity new is validObject slot
#' @importFrom S4Vectors DataFrame
NULL

# ============================================================
# ISNVFilter -- QC filter configuration class
# ============================================================

#' ISNVFilter: Quality control filter for iSNV analysis
#'
#' An S4 class holding parameters for filtering intra-host single
#' nucleotide variants. Used by \code{\link{flagISNV}}.
#'
#' @slot minDepth Integer scalar. Minimum read depth at iSNV site
#'   (default: 100).
#' @slot minFreq Numeric scalar. Minimum alternative allele frequency
#'   (default: 0.03).
#' @slot maxFreq Numeric scalar. Maximum alternative allele frequency;
#'   variants above this are excluded as near-fixed (default: 0.50).
#' @slot minAltReads Integer scalar. Minimum reads supporting the
#'   alternative allele (default: 10).
#' @slot maxStrandBias Numeric scalar. Maximum strand bias as
#'   fold-difference between strands (default: 10.0).
#' @slot replicateConc Logical scalar. Whether replicate concordance
#'   is required (default: FALSE). Recorded for provenance only;
#'   apply it with \code{\link{flagReplicateDiscordance}}.
#'
#' @return An object of class \code{ISNVFilter}. Objects are created with
#'   the constructor \code{\link{ISNVFilter}}; the \code{show} method prints
#'   the filter thresholds and returns \code{NULL} invisibly.
#'
#' @details
#' Default values are based on commonly used thresholds in the
#' within-host pathogen diversity literature:
#' \describe{
#'   \item{\code{minDepth = 100}}{Minimum coverage for variant
#'     detection. This is a permissive floor: Grubaugh et al. (2019,
#'     Genome Biology) measured iSNVs above 3\% accurately from at
#'     least 1,000 RNA copies sequenced to at least 400x.}
#'   \item{\code{minFreq = 0.03}}{3\% minor allele frequency
#'     threshold, the lower limit of accurate iSNV measurement in
#'     Grubaugh et al. (2019); Popa et al. (2020, Science
#'     Translational Medicine) used 2\%.}
#'   \item{\code{maxFreq = 0.50}}{Variants above 50\% are considered
#'     consensus-level changes, not intra-host variants.}
#'   \item{\code{minAltReads = 10}}{Minimum supporting reads to
#'     reduce false positives from sequencing errors.}
#'   \item{\code{maxStrandBias = 10}}{Maximum fold-difference between
#'     forward and reverse strand alt allele counts. Excessive strand
#'     bias suggests systematic error.}
#' }
#'
#' @references
#' Roder AE et al. (2023). Optimized quantification of intra-host
#' viral diversity in SARS-CoV-2 and influenza virus sequence data.
#' \emph{mBio} 14:e01046-23.
#' \doi{10.1128/mbio.01046-23}
#'
#' McCrone JT, Lauring AS (2016). Measurements of intrahost viral
#' diversity are extremely sensitive to systematic errors in variant
#' calling. \emph{J Virol} 90:6884-6895.
#' \doi{10.1128/JVI.00667-16}
#'
#' Grubaugh ND et al. (2019). An amplicon-based sequencing framework
#' for accurately measuring intrahost virus diversity using PrimalSeq
#' and iVar. \emph{Genome Biology} 20:8.
#' \doi{10.1186/s13059-018-1618-7}
#'
#' Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host
#' variation by leveraging large-scale sequencing data.
#' \emph{NAR Genomics and Bioinformatics} 6:lqae145.
#' \doi{10.1093/nargab/lqae145}
#'
#' @seealso \code{\link{ISNVFilter}} for the constructor,
#'   \code{\link{flagISNV}} for applying filters.
#'
#' @exportClass ISNVFilter
#' @examples
#' # Default filter
#' ISNVFilter()
#'
#' # Stringent filter
#' ISNVFilter(minDepth = 500L, minFreq = 0.05)
.ISNVFilter <- setClass("ISNVFilter",
    slots = list(
        minDepth      = "integer",
        minFreq       = "numeric",
        maxFreq       = "numeric",
        minAltReads   = "integer",
        maxStrandBias = "numeric",
        replicateConc = "logical"
    ),
    prototype = list(
        minDepth      = 100L,
        minFreq       = 0.03,
        maxFreq       = 0.50,
        minAltReads   = 10L,
        maxStrandBias = 10.0,
        replicateConc = FALSE
    )
)

setValidity("ISNVFilter", function(object) {
    msg <- character()

    # minDepth
    v <- slot(object, "minDepth")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'minDepth' must be a single non-NA integer")
    } else if (v < 0L) {
        msg <- c(msg, sprintf("'minDepth' must be non-negative, got %d", v))
    }

    # minFreq
    v <- slot(object, "minFreq")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'minFreq' must be a single non-NA numeric")
    } else if (v < 0 || v > 1) {
        msg <- c(msg, sprintf(
            "'minFreq' must be between 0 and 1, got %.4f", v))
    }

    # maxFreq
    v <- slot(object, "maxFreq")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'maxFreq' must be a single non-NA numeric")
    } else if (v > 1) {
        msg <- c(msg, sprintf("'maxFreq' must be <= 1, got %.4f", v))
    }

    # cross-check
    mn <- slot(object, "minFreq")
    mx <- slot(object, "maxFreq")
    if (length(mn) == 1L && length(mx) == 1L &&
        !is.na(mn) && !is.na(mx) && mx < mn) {
        msg <- c(msg, sprintf(
            "'maxFreq' (%.4f) must be >= 'minFreq' (%.4f)", mx, mn))
    }

    # minAltReads
    v <- slot(object, "minAltReads")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'minAltReads' must be a single non-NA integer")
    } else if (v < 0L) {
        msg <- c(msg, sprintf(
            "'minAltReads' must be non-negative, got %d", v))
    }

    # maxStrandBias
    v <- slot(object, "maxStrandBias")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'maxStrandBias' must be a single non-NA numeric")
    } else if (v <= 0) {
        msg <- c(msg, sprintf(
            "'maxStrandBias' must be positive, got %.4f", v))
    }

    # replicateConc
    v <- slot(object, "replicateConc")
    if (length(v) != 1L || is.na(v)) {
        msg <- c(msg, "'replicateConc' must be a single non-NA logical")
    }

    if (length(msg) == 0L) TRUE else msg
})

# ============================================================
# WithinHostExperiment -- main container class
# ============================================================

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
NULL

#' WithinHostExperiment: Container for within-host pathogen variant data
#'
#' Extends \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' to hold intra-host single nucleotide variant (iSNV) data with
#' integrated quality control audit logging.
#'
#' @section Standard assays:
#' \describe{
#'   \item{\code{altCount}}{Integer matrix. Alternative allele read counts.}
#'   \item{\code{refCount}}{Integer matrix. Reference allele read counts.}
#'   \item{\code{totalDepth}}{Integer matrix. Total read depth.}
#'   \item{\code{altFreq}}{Numeric matrix. Alt allele frequencies (0-1).}
#'   \item{\code{qcPass}}{Logical matrix. TRUE = passed QC, FALSE = flagged.
#'     Initialized to all TRUE; updated by \code{\link{flagISNV}}.}
#' }
#'
#' @section Expected colData columns:
#' \describe{
#'   \item{\code{sample_id}}{Character. Unique sample identifier (required).}
#'   \item{\code{host_id}}{Character. Host/patient identifier (optional).}
#'   \item{\code{replicate_group}}{Character. Technical replicate group
#'     identifier (optional).}
#'   \item{\code{replicate_id}}{Character. Replicate number within group
#'     (optional).}
#'   \item{\code{role}}{Character. Sample role: "donor", "recipient",
#'     or "independent" (optional).}
#'   \item{\code{pair_id}}{Character. Transmission pair identifier (optional).}
#' }
#'
#' @slot qcLog \code{DataFrame}. Audit trail of QC filtering steps.
#' @slot qcFilters \code{list}. History of \code{\link{ISNVFilter}} objects
#'   that have been applied.
#'
#' @seealso
#'   \code{\link{readWithinHost}} and \code{\link{readWithinHostTable}}
#'   for constructing from files.
#'   \code{\link{flagISNV}} for quality control.
#'   \code{\link{passedISNV}} for extracting QC-passed variants.
#'   \code{\link{transmissionPairs}} for extracting donor-recipient pairs.
#'
#' @examples
#' data(example_whe)
#' example_whe
#'
#' @exportClass WithinHostExperiment
.WithinHostExperiment <- setClass("WithinHostExperiment",
    contains = "RangedSummarizedExperiment",
    slots = list(
        qcLog     = "DataFrame",
        qcFilters = "list"
    ),
    prototype = list(
        qcLog = DataFrame(
            step      = character(),
            parameter = character(),
            value     = character(),
            n_flagged = integer(),
            n_passed  = integer(),
            timestamp = character()
        ),
        qcFilters = list()
    )
)

setValidity("WithinHostExperiment", function(object) {
    msg <- character()

    ## Validate qcLog structure
    log <- object@qcLog
    expected_cols <- c("step", "parameter", "value",
                       "n_flagged", "n_passed", "timestamp")
    if (!all(expected_cols %in% colnames(log))) {
        missing <- setdiff(expected_cols, colnames(log))
        msg <- c(msg, paste0("'qcLog' missing columns: ",
                             paste(missing, collapse = ", ")))
    }

    ## Validate qcFilters is a list of ISNVFilter objects (if non-empty)
    filters <- object@qcFilters
    if (length(filters) > 0L) {
        are_filters <- vapply(filters, is, logical(1), "ISNVFilter")
        if (!all(are_filters)) {
            msg <- c(msg,
                "'qcFilters' must contain only ISNVFilter objects")
        }
    }

    ## Validate required assays
    anames <- SummarizedExperiment::assayNames(object)
    if (length(anames) > 0L && !"altFreq" %in% anames) {
        msg <- c(msg, "'altFreq' assay is required")
    }

    ## Validate qcPass assay is logical if present
    if ("qcPass" %in% anames) {
        qc <- SummarizedExperiment::assay(object, "qcPass")
        if (!is.logical(qc)) {
            msg <- c(msg, "'qcPass' assay must be a logical matrix")
        }
    }

    if (length(msg) == 0L) TRUE else msg
})
