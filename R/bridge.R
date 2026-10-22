# R/bridge.R
# Bridge to downstream tools (ViralBottleneck, etc.)

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors DataFrame mcols
NULL

#' Export donor-recipient frequency table
#'
#' For each transmission pair, creates a table of variant frequencies
#' in the donor and recipient.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with \code{role}
#'   and \code{pair_id} in \code{colData}.
#' @param pairId Character (optional). Specific pair to export.
#'   If NULL, exports all pairs.
#' @param usePassedOnly Logical. Use only QC-passed variants
#'   (default: TRUE).
#'
#' @return A \code{DataFrame} with columns: chrom, position, ref, alt,
#'   donor_freq, recipient_freq, shared (logical).
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' exportPairFrequencies(whe)
exportPairFrequencies <- function(whe, pairId = NULL,
                                  usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    pairs <- transmissionPairs(whe)
    if (nrow(pairs) == 0L) stop("No transmission pairs defined in colData")

    if (!is.null(pairId)) {
        pairs <- pairs[pairs$pair_id == pairId, , drop = FALSE]
        if (nrow(pairs) == 0L) stop("Pair '", pairId, "' not found")
    }

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if ("qcPass" %in% assayNames(whe) && usePassedOnly) {
        assay(whe, "qcPass")
    } else {
        NULL
    }

    sample_ids <- as.character(colData(whe)$sample_id)
    rr <- rowRanges(whe)
    ref_col <- if ("ref" %in% colnames(mcols(rr))) {
        mcols(rr)$ref
    } else {
        rep(NA_character_, length(rr))
    }
    alt_col <- if ("alt" %in% colnames(mcols(rr))) {
        mcols(rr)$alt
    } else {
        rep(NA_character_, length(rr))
    }

    result_list <- lapply(seq_len(nrow(pairs)), function(p) {
        d_col <- match(as.character(pairs$donor[p]), sample_ids)
        r_col <- match(as.character(pairs$recipient[p]), sample_ids)

        d_freq <- freq_mat[, d_col]
        r_freq <- freq_mat[, r_col]

        ## Apply QC mask
        if (!is.null(qc_mat)) {
            d_freq[!qc_mat[, d_col]] <- NA
            r_freq[!qc_mat[, r_col]] <- NA
        }

        ## Keep rows where at least one has data
        has_data <- !is.na(d_freq) | !is.na(r_freq)

        DataFrame(
            chrom          = as.character(seqnames(rr))[has_data],
            position       = start(rr)[has_data],
            ref            = ref_col[has_data],
            alt            = alt_col[has_data],
            donor_freq     = d_freq[has_data],
            recipient_freq = r_freq[has_data],
            shared         = !is.na(d_freq[has_data]) &
                !is.na(r_freq[has_data]),
            pair_id        = as.character(pairs$pair_id[p])
        )
    })

    do.call(rbind, result_list)
}

#' Export data for ViralBottleneck package
#'
#' Converts a \code{\link{WithinHostExperiment}} to the input format
#' expected by the \pkg{ViralBottleneck} package.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param pairId Character (optional). Specific pair.
#' @param threshold Numeric. Variant calling threshold (default: 0.03).
#' @param usePassedOnly Logical (default: TRUE).
#'
#' @return A \code{data.frame} with columns: pos, donor_freq,
#'   recipient_freq -- a simplified format for beta-binomial MLE.
#'   For the full ViralBottleneck input format (per-base read
#'   counts), prepare data externally using the read count
#'   assays (\code{altCount}, \code{refCount}).
#'
#' @details
#' For comprehensive bottleneck analysis with multiple methods
#' (presence-absence, KL divergence, binomial, beta-binomial
#' approximate/exact, and Wright-Fisher), use the
#' \pkg{ViralBottleneck} package (Zheng et al. 2025,
#' \emph{Virus Evolution} 11:veaf071). This function prepares
#' a simplified input suitable for the beta-binomial methods.
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' asViralBottleneckInput(whe, pairId = "p1")
asViralBottleneckInput <- function(whe, pairId = NULL, threshold = 0.03,
                                   usePassedOnly = TRUE) {
    pf <- exportPairFrequencies(whe, pairId = pairId,
                                usePassedOnly = usePassedOnly)

    ## Filter to donor variants above threshold
    donor_above <- !is.na(pf$donor_freq) & pf$donor_freq >= threshold
    pf <- pf[donor_above, , drop = FALSE]

    ## Replace NA recipient freq with 0
    recip_freq <- pf$recipient_freq
    recip_freq[is.na(recip_freq)] <- 0

    data.frame(
        pos            = pf$position,
        donor_freq     = pf$donor_freq,
        recipient_freq = recip_freq,
        stringsAsFactors = FALSE
    )
}

#' Calculate shared variant statistics
#'
#' Counts, for one transmission pair, the variant sites present in both
#' hosts, in the donor alone and in the recipient alone, using the calls
#' that passed QC.
#'
#' @details
#' A donor variant absent from the recipient means one thing when the
#' recipient was sequenced deeply enough at that site to have shown it,
#' and another when it was not. The donor-only count is therefore
#' split by the recipient's read depth at the site, from the
#' \code{totalDepth} assay: \code{n_donor_only_adequate_depth} where
#' that depth is at least \code{minDepth},
#' \code{n_donor_only_low_depth} where it is lower, and
#' \code{n_donor_only_unknown_depth} where it is missing. Depth is
#' missing at every site the recipient's caller did not report when the
#' object was read from VCF files, which list only the sites with a
#' call, so with such input most donor-only sites fall in the last
#' group; a depth matrix from a coverage file, or a table listing every
#' assessed site, is needed to tell the first two apart. Only a
#' donor-only site at adequate depth is evidence of loss in
#' transmission.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param pairId Character. Which pair to analyse.
#' @param minDepth Integer. The recipient depth at or above which an
#'   absent donor variant counts as absent at adequate depth. Defaults
#'   to the strictest \code{minDepth} of the filters applied with
#'   \code{\link{flagISNV}}, or to 100 when none was applied.
#'
#' @return A named list: \code{n_donor_only} (the three groups below
#'   added together), \code{n_donor_only_adequate_depth},
#'   \code{n_donor_only_low_depth}, \code{n_donor_only_unknown_depth},
#'   \code{n_recipient_only}, \code{n_shared}, \code{n_total},
#'   \code{sharing_rate} (shared sites over all sites) and
#'   \code{min_depth}, the threshold used.
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' calcSharedVariants(whe, pairId = "p1")
calcSharedVariants <- function(whe, pairId, minDepth = NULL) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (is.null(minDepth)) {
        filters <- qcFilters(whe)
        minDepth <- if (length(filters)) {
            max(vapply(filters, function(f) as.integer(slot(f, "minDepth")),
                       integer(1)))
        } else {
            slot(ISNVFilter(), "minDepth")
        }
    }
    if (!is.numeric(minDepth) || length(minDepth) != 1L || is.na(minDepth) ||
        minDepth < 0)
        stop("'minDepth' must be a single non-negative number.")
    pf <- exportPairFrequencies(whe, pairId = pairId)

    d_only <- !is.na(pf$donor_freq) & is.na(pf$recipient_freq)
    r_only <- is.na(pf$donor_freq) & !is.na(pf$recipient_freq)
    shared <- sum(pf$shared)
    total  <- nrow(pf)

    ## the recipient's depth at the donor-only sites
    pairs <- transmissionPairs(whe)
    recipient <- as.character(pairs$recipient[pairs$pair_id == pairId][1L])
    r_col <- match(recipient, as.character(colData(whe)$sample_id))
    rr <- rowRanges(whe)
    key_all <- paste(as.character(seqnames(rr)), start(rr),
                     if ("ref" %in% colnames(mcols(rr))) mcols(rr)$ref else NA,
                     if ("alt" %in% colnames(mcols(rr))) mcols(rr)$alt else NA)
    key_pf <- paste(pf$chrom, pf$position, pf$ref, pf$alt)
    r_depth <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")[match(key_pf, key_all), r_col]
    } else {
        rep(NA_integer_, nrow(pf))
    }
    adequate <- d_only & !is.na(r_depth) & r_depth >= minDepth
    low      <- d_only & !is.na(r_depth) & r_depth < minDepth
    unknown  <- d_only & is.na(r_depth)

    list(
        n_donor_only                = sum(d_only),
        n_donor_only_adequate_depth = sum(adequate),
        n_donor_only_low_depth      = sum(low),
        n_donor_only_unknown_depth  = sum(unknown),
        n_recipient_only            = sum(r_only),
        n_shared                    = shared,
        n_total                     = total,
        sharing_rate     = if (total > 0L) shared / total else NA_real_,
        min_depth        = as.integer(minDepth)
    )
}
