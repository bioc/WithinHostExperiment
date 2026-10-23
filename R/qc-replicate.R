# R/qc-replicate.R
# Replicate-aware QC

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames colData
#' @importFrom S4Vectors DataFrame
NULL

#' Flag iSNVs with replicate discordance
#'
#' Identifies variant sites where technical replicates disagree and
#' flags them in the \code{qcPass} assay.
#'
#' @param whe A \code{\linkS4class{WithinHostExperiment}} with
#'   \code{replicate_group} and \code{replicate_id} columns in
#'   \code{colData}.
#' @param freqTolerance Numeric scalar. Maximum allowable difference
#'   in \code{altFreq} between replicates (default: 0.05).
#' @param requireBothDetected Logical scalar. If \code{TRUE} (default),
#'   sites detected in only one replicate are also flagged.
#'
#' @return A \code{\linkS4class{WithinHostExperiment}} with updated
#'   \code{qcPass} and \code{qcLog}.
#'
#' @details
#' For each replicate group, the function compares \code{altFreq}
#' between replicates and flags sites where the difference exceeds
#' \code{freqTolerance}. Sites are flagged in ALL samples of the
#' replicate group.
#'
#' Roder et al. (2023) reported that the choice of variant caller and
#' the use of replicate sequencing had the largest effects on
#' single-nucleotide variant discovery. An
#' alternative approach is replicate intersection (requiring
#' detection in both replicates); use \code{requireBothDetected =
#' TRUE} with \code{freqTolerance = Inf} to achieve this.
#'
#' @references
#' Roder AE et al. (2023). Optimized quantification of intra-host
#' viral diversity in SARS-CoV-2 and influenza virus sequence data.
#' \emph{mBio} 14:e01046-23.
#' \doi{10.1128/mbio.01046-23}
#'
#' @seealso \code{\link{flagISNV}} for standard QC filtering.
#'
#' @export
#' @examples
#' rep1 <- system.file("extdata", "test_donor_rep1.vcf",
#'     package = "WithinHostExperiment")
#' rep2 <- system.file("extdata", "test_donor_rep2.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(rep1, rep2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("rep1", "rep2"),
#'         replicate_group = c("donor_reps", "donor_reps"),
#'         replicate_id = c("1", "2")),
#'     caller = "ivar")
#' whe <- flagReplicateDiscordance(whe, freqTolerance = 0.05)
#' qcLog(whe)
flagReplicateDiscordance <- function(whe, freqTolerance = 0.05,
                                     requireBothDetected = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!is.numeric(freqTolerance) || length(freqTolerance) != 1L ||
        freqTolerance <= 0)
        stop("'freqTolerance' must be a single positive number.")

    cd <- colData(whe)

    if (!"replicate_group" %in% colnames(cd)) {
        warning("No 'replicate_group' column in colData; ",
                "skipping replicate discordance check.")
        return(whe)
    }

    freq_mat <- assay(whe, "altFreq")
    qc <- assay(whe, "qcPass")

    groups <- unique(as.character(cd$replicate_group))
    groups <- groups[!is.na(groups)]

    if (length(groups) == 0L) {
        warning("No non-NA replicate_group values found; ",
                "skipping replicate discordance check.")
        return(whe)
    }

    res <- .flagDiscordantReplicates(cd$replicate_group, freq_mat, qc,
                                     freqTolerance, requireBothDetected)
    qc <- res$qc

    assay(whe, "qcPass") <- qc
    whe <- .addQcStep(whe, "replicate_discordance", "freqTolerance",
                      freqTolerance, res$n_flagged,
                      sum(qc, na.rm = TRUE))
    whe <- .addQcStep(whe, "replicate_discordance",
                      "requireBothDetected",
                      requireBothDetected, 0L,
                      sum(qc, na.rm = TRUE))
    whe
}

## Compares every pair of samples sharing a value of 'group' and sets
## qcPass to FALSE, in every sample of the group, at each site that fails
## for any pair. A call is a non-missing frequency. Returns the updated
## qcPass matrix and the number of cells newly flagged. Shared by
## flagReplicateDiscordance() and the ConcordanceRule step of flagISNV().
.flagDiscordantReplicates <- function(group, freq_mat, qc, freqTolerance,
                                      requireBothDetected) {
    group <- as.character(group)
    groups <- unique(group[!is.na(group)])
    total_flagged <- 0L

    for (grp in groups) {
        grp_idx <- which(group == grp & !is.na(group))
        if (length(grp_idx) < 2L) next

        n_reps <- length(grp_idx)
        for (i in seq_len(n_reps - 1L)) {
            for (j in (i + 1L):n_reps) {
                ci <- grp_idx[i]
                cj <- grp_idx[j]

                freq_i <- freq_mat[, ci]
                freq_j <- freq_mat[, cj]

                both_present <- !is.na(freq_i) & !is.na(freq_j)
                freq_diff <- abs(freq_i - freq_j)
                discordant_freq <- both_present &
                    freq_diff > freqTolerance

                if (requireBothDetected) {
                    one_only <- xor(is.na(freq_i), is.na(freq_j))
                    discordant <- discordant_freq | one_only
                } else {
                    discordant <- discordant_freq
                }

                if (any(discordant)) {
                    for (k in grp_idx) {
                        newly_flagged <- discordant & qc[, k]
                        qc[newly_flagged, k] <- FALSE
                        total_flagged <- total_flagged +
                            sum(newly_flagged)
                    }
                }
            }
        }
    }

    list(qc = qc, n_flagged = total_flagged)
}
