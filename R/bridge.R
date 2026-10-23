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
#'   donor_freq, recipient_freq, shared (logical) and pair_id. A host
#'   carries the allele at a site when its frequency there is known and
#'   above 0; \code{shared} is \code{TRUE} where both hosts carry it.
#'   A recorded frequency of 0, such as a site that a table lists
#'   without alternative reads, is an assessed site without the allele:
#'   it is kept as 0 in its frequency column, and is not shared. A
#'   frequency is \code{NA} where the host has no call, or its call
#'   failed QC. Sites that neither host carries are left out.
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
    pf <- .pairTable(whe, pairId, usePassedOnly)
    pf[, setdiff(colnames(pf), c("row_index", "recipient_col")),
       drop = FALSE]
}

## The table of exportPairFrequencies(), with the row of the object and
## the recipient's column that each line comes from
.pairTable <- function(whe, pairId, usePassedOnly) {
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

        ## A host carries the allele where its frequency is known and
        ## above 0; a recorded 0 is an assessed site without the allele
        d_pres <- .carries(d_freq)
        r_pres <- .carries(r_freq)
        has_data <- d_pres | r_pres

        DataFrame(
            chrom          = as.character(seqnames(rr))[has_data],
            position       = start(rr)[has_data],
            ref            = ref_col[has_data],
            alt            = alt_col[has_data],
            donor_freq     = d_freq[has_data],
            recipient_freq = r_freq[has_data],
            shared         = (d_pres & r_pres)[has_data],
            pair_id        = rep(as.character(pairs$pair_id[p]),
                                 sum(has_data)),
            row_index      = unname(which(has_data)),
            recipient_col  = rep(r_col, sum(has_data))
        )
    })

    do.call(rbind, result_list)
}

## TRUE where a frequency is known and above 0
.carries <- function(freq) !is.na(freq) & freq > 0

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
#' One row is written for each donor variant at or above
#' \code{threshold}. Where the recipient has no call, its frequency is
#' written as 0, which an estimator reads as a variant lost in
#' transmission. A recipient call at or above \code{threshold} that QC
#' removed is not such an absence: the call may have been removed for
#' lying above \code{maxFreq}, as when the variant rose towards fixation
#' in the recipient. With \code{usePassedOnly = TRUE}, a donor variant at
#' such a site is therefore left out, and a message gives the number of
#' sites. A removed call below \code{threshold} is an absence by the
#' estimators' own definition, and is written as 0.
#' (\code{n_donor_only_recipient_failed_qc} in
#' \code{\link{calcSharedVariants}} counts every removed call above 0,
#' whatever its frequency.) \code{\link{quickBottleneck}} and
#' \code{\link{exactBottleneck}} use the same donor variants.
#'
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
    .bottleneckInput(whe, pairId, threshold, usePassedOnly)$table
}

## The donor variants a bottleneck estimate uses: the table of
## asViralBottleneckInput(), the lines of .pairTable() it comes from,
## and the number of sites left out because QC removed the recipient's
## call there
.bottleneckInput <- function(whe, pairId, threshold,
                             usePassedOnly = TRUE) {
    pf <- .pairTable(whe, pairId, usePassedOnly)

    ## Filter to donor variants above threshold
    donor_above <- !is.na(pf$donor_freq) & pf$donor_freq >= threshold
    pf <- pf[donor_above, , drop = FALSE]

    ## A recipient call at or above threshold that QC removed is not an
    ## absence; one below threshold is, by the estimators' definition
    failed <- logical(nrow(pf))
    if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        r_raw <- assay(whe, "altFreq")[cbind(pf$row_index,
                                             pf$recipient_col)]
        failed <- is.na(pf$recipient_freq) & !is.na(r_raw) &
            r_raw >= threshold
    }
    n_failed <- sum(failed)
    if (n_failed > 0L)
        message(n_failed, " donor variant site(s) where QC removed a ",
                "recipient call at or above threshold were left out: ",
                "such a site is not an absence.")
    pf <- pf[!failed, , drop = FALSE]

    ## Elsewhere a recipient without a call is written as 0
    recip_freq <- pf$recipient_freq
    recip_freq[is.na(recip_freq)] <- 0

    list(
        table = data.frame(
            pos            = pf$position,
            donor_freq     = pf$donor_freq,
            recipient_freq = recip_freq,
            stringsAsFactors = FALSE
        ),
        sites = pf,
        n_failed_qc = n_failed
    )
}

#' Calculate shared variant statistics
#'
#' Counts, for one transmission pair, the variant sites present in both
#' hosts, in the donor alone and in the recipient alone, using the calls
#' that passed QC.
#'
#' @details
#' A host carries a variant at a site when its frequency there is known
#' and above 0, as in \code{\link{exportPairFrequencies}}; a recorded
#' frequency of 0 is a site that was assessed and found without the
#' variant.
#'
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
#' assessed site with its depth and a frequency of 0 where the recipient
#' had no alternative reads (\code{\link{readWithinHostTable}} with
#' \code{format = "generic"} and \code{aaf = 0}), is needed to tell the
#' first two apart. A fourth group,
#' \code{n_donor_only_recipient_failed_qc}, holds the sites where the
#' recipient did make a call, at a frequency above 0, but QC removed
#' it, for example a call
#' above \code{maxFreq} in a recipient where the variant rose towards
#' fixation: such a site is not absent, so it is kept out of the three
#' depth groups. Only a donor-only site at adequate depth is evidence of
#' loss in transmission.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param pairId Character. Which pair to analyse.
#' @param minDepth Integer. The recipient depth at or above which an
#'   absent donor variant counts as absent at adequate depth. Defaults
#'   to the strictest \code{minDepth} of the filters applied with
#'   \code{\link{flagISNV}}, or to 100 when none was applied.
#'
#' @return A named list: \code{n_donor_only} (the four groups below
#'   added together), \code{n_donor_only_adequate_depth},
#'   \code{n_donor_only_low_depth}, \code{n_donor_only_unknown_depth},
#'   \code{n_donor_only_recipient_failed_qc},
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

    d_pres <- .carries(pf$donor_freq)
    r_pres <- .carries(pf$recipient_freq)
    d_only <- d_pres & !r_pres
    r_only <- r_pres & !d_pres
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
    ## a recipient call at the site that QC removed is not an absence; a
    ## recorded frequency of 0 is not a call
    r_call <- assay(whe, "altFreq")[match(key_pf, key_all), r_col]
    failed   <- d_only & .carries(r_call)
    adequate <- d_only & !failed & !is.na(r_depth) & r_depth >= minDepth
    low      <- d_only & !failed & !is.na(r_depth) & r_depth < minDepth
    unknown  <- d_only & !failed & is.na(r_depth)

    list(
        n_donor_only                = sum(d_only),
        n_donor_only_adequate_depth = sum(adequate),
        n_donor_only_low_depth      = sum(low),
        n_donor_only_unknown_depth  = sum(unknown),
        n_donor_only_recipient_failed_qc = sum(failed),
        n_recipient_only            = sum(r_only),
        n_shared                    = shared,
        n_total                     = total,
        sharing_rate     = if (total > 0L) shared / total else NA_real_,
        min_depth        = as.integer(minDepth)
    )
}
