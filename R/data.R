# R/data.R
# Documentation for example datasets

#' Example WithinHostExperiment dataset
#'
#' A small synthetic \code{\link{WithinHostExperiment}} with 3 samples
#' (donor, recipient, independent) containing iSNV data at positions 33 to
#' 945 of a synthetic genome segment. The donor and recipient form a
#' transmission pair (pair_id = "pair_1").
#'
#' @format A \code{\link{WithinHostExperiment}} with 28 variant sites
#'   and 3 samples. Assays: altFreq, totalDepth, refCount, altCount,
#'   qcPass. The length of the segment, 1000 bases, is recorded in
#'   \code{seqlengths()} from \code{test_ref.fa}.
#'
#' @source Synthetic. Built by \code{inst/scripts/save_example_data.R} from
#'   \code{test_donor.vcf}, \code{test_recipient.vcf} and
#'   \code{test_independent.vcf} in \code{inst/extdata}, read with
#'   \code{\link{readWithinHost}}.
#'
#' @docType data
#' @name example_whe
#' @usage data(example_whe)
#'
#' @examples
#' data(example_whe)
#' example_whe
#' transmissionPairs(example_whe)
NULL
