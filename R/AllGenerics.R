# R/AllGenerics.R
# Generic function definitions
#
# Design principle: every analytical operation is a generic.
# This allows downstream packages to specialize methods for
# their own container classes (e.g., a LongitudinalWHE or
# MultiPathogenWHE), following the same pattern as
# SummarizedExperiment's assay/colData/rowRanges generics.

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

# ---- WithinHostExperiment accessors ----

#' Access the QC audit log
#'
#' Returns the quality control log recording all filtering steps
#' applied to a \code{\link{WithinHostExperiment}}.
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @return A \code{DataFrame} with columns: step, parameter, value,
#'   n_flagged, n_passed, timestamp.
#'
#' @export
#' @rdname WHE-accessors
#' @examples
#' gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3, width = 1))
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(runif(3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' qcLog(whe)
setGeneric("qcLog", function(x) standardGeneric("qcLog"))

#' Access the list of applied QC filters
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @return A list of \code{\link{ISNVFilter}} objects.
#'
#' @export
#' @rdname WHE-accessors
setGeneric("qcFilters", function(x) standardGeneric("qcFilters"))

#' Extract QC-passed variants
#'
#' Returns a subset of the \code{\link{WithinHostExperiment}} containing
#' only variant sites where at least one sample passed QC.
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @return A \code{\link{WithinHostExperiment}} (subset of rows).
#'
#' @export
#' @rdname WHE-accessors
#' @examples
#' gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3, width = 1))
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(runif(3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' passedISNV(whe)
setGeneric("passedISNV", function(x) standardGeneric("passedISNV"))

#' Extract transmission pairs from colData
#'
#' Identifies donor-recipient pairs based on \code{role} and
#' \code{pair_id} columns in \code{colData}.
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @return A \code{DataFrame} with columns: donor (sample_id),
#'   recipient (sample_id), pair_id.
#'
#' @export
#' @rdname WHE-accessors
setGeneric("transmissionPairs",
    function(x) standardGeneric("transmissionPairs"))

#' Set the QC audit log
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @param value A \code{DataFrame} replacement value.
#' @return The modified \code{\link{WithinHostExperiment}}.
#'
#' @export
#' @rdname WHE-accessors
setGeneric("qcLog<-", function(x, value) standardGeneric("qcLog<-"))

#' Set the list of applied QC filters
#'
#' @param x A \code{\link{WithinHostExperiment}} object.
#' @param value A list replacement value.
#' @return The modified \code{\link{WithinHostExperiment}}.
#'
#' @export
#' @rdname WHE-accessors
setGeneric("qcFilters<-", function(x, value) standardGeneric("qcFilters<-"))

# ---- Analytical generics ----
# These are extension points: downstream packages can define
# methods for their own container classes.

#' Flag iSNV sites based on quality criteria
#'
#' Generic for applying quality control filters to within-host
#' variant data. The default method operates on
#' \code{\link{WithinHostExperiment}} objects. Other packages
#' can define methods for alternative containers.
#'
#' @param x A within-host variant container.
#' @param filter An \code{\link{ISNVFilter}} object.
#' @param ... Additional arguments for methods.
#' @return The filtered object (same class as input).
#'
#' @export
#' @rdname flagISNV
setGeneric("flagISNV", function(x, filter = ISNVFilter(), ...)
    standardGeneric("flagISNV"))

#' Compute within-host diversity indices
#'
#' Generic for computing diversity metrics from within-host
#' variant data. Methods should return a \code{DataFrame} with
#' one row per sample.
#'
#' @param x A within-host variant container.
#' @param ... Additional arguments (e.g., \code{genomeLength},
#'   \code{indices}).
#' @return A \code{DataFrame} of diversity metrics.
#'
#' @export
#' @rdname calcDiversity
setGeneric("calcDiversity", function(x, ...)
    standardGeneric("calcDiversity"))

#' @export
#' @rdname annotateGenes
setGeneric("annotateGenes", function(x, genes, ...)
    standardGeneric("annotateGenes"))

#' @export
#' @rdname annotateCodonChange
setGeneric("annotateCodonChange", function(x, cds, ref, ...)
    standardGeneric("annotateCodonChange"))

#' @export
#' @rdname pNpSWithinHost
setGeneric("pNpSWithinHost", function(x, cds, ref, ...)
    standardGeneric("pNpSWithinHost"))

#' @export
#' @rdname consensusSequences
setGeneric("consensusSequences", function(x, ref, ...)
    standardGeneric("consensusSequences"))
