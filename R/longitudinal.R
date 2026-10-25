# R/longitudinal.R
# Longitudinal within-host variant tracking

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors DataFrame mcols
NULL

## Avoid R CMD check notes for ggplot2 .data pronoun
utils::globalVariables(".data")

# ============================================================
# trackFrequency -- extract frequency trajectories
# ============================================================

#' Track variant frequency trajectories over time
#'
#' Extracts per-variant allele frequency trajectories across
#' longitudinal time points within each host. Returns a long-format
#' \code{DataFrame} suitable for downstream analysis or plotting.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with longitudinal
#'   samples from the same host.
#' @param hostCol Character scalar. Column name in \code{colData}
#'   containing the host/patient identifier (default:
#'   \code{"host_id"}).
#' @param timeCol Character scalar. Column name in \code{colData}
#'   containing the time point (default: \code{"timepoint"}). It must
#'   hold numbers, dates, date-times, time differences or a factor whose
#'   levels are in time order; text is read as numbers when every value
#'   is one, and is otherwise an error, as in
#'   \code{\link{flagTemporalInconsistency}}. Samples whose time point is
#'   \code{NA} are left out, with a warning.
#' @param usePassedOnly Logical. If TRUE (default), only use
#'   variants where \code{qcPass == TRUE}. With \code{FALSE} every call
#'   is returned and \code{qc_passed} marks those that failed.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{host_id}{Host identifier.}
#'     \item{variant_key}{Unique variant identifier
#'       (chrom:pos:alt).}
#'     \item{timepoint}{Time point value.}
#'     \item{frequency}{Alternative allele frequency at that
#'       time point.}
#'     \item{qc_passed}{Logical. Whether the call passed QC at that
#'       time point; an \code{NA} \code{qcPass} counts as not passed
#'       (\code{TRUE} if no QC assay present).}
#'   }
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C")
#' S4Vectors::mcols(gr)$alt <- c("T", "G")
#' freq <- matrix(c(0.05, 0.10, 0.15, 0.20), nrow = 2, ncol = 2)
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freq),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("S1_t1", "S1_t2"),
#'         host_id = c("H1", "H1"),
#'         timepoint = c(1, 2)))
#' trackFrequency(whe)
trackFrequency <- function(whe,
                           hostCol = "host_id",
                           timeCol = "timepoint",
                           usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    cd <- colData(whe)
    if (!hostCol %in% colnames(cd))
        stop("Column '", hostCol, "' not found in colData. ",
             "Available columns: ",
             paste(colnames(cd), collapse = ", "))
    if (!timeCol %in% colnames(cd))
        stop("Column '", timeCol, "' not found in colData. ",
             "Available columns: ",
             paste(colnames(cd), collapse = ", "))

    freq_mat <- assay(whe, "altFreq")
    ## Each cell's QC state is read whatever usePassedOnly is, so that with
    ## FALSE the failing calls are kept and marked
    passed <- if ("qcPass" %in% assayNames(whe)) {
        .passedQC(assay(whe, "qcPass"))
    } else {
        matrix(TRUE, nrow = nrow(freq_mat), ncol = ncol(freq_mat))
    }

    rr <- rowRanges(whe)
    mc <- mcols(rr)
    chr_vals <- as.character(seqnames(rr))
    pos_vals <- start(rr)
    alt_vals <- if ("alt" %in% colnames(mc)) {
        as.character(mc$alt)
    } else {
        rep(NA_character_, length(rr))
    }
    variant_keys <- paste(chr_vals, pos_vals, alt_vals, sep = ":")

    host_ids <- as.character(cd[[hostCol]])
    timepoints <- cd[[timeCol]]
    n_sites <- nrow(freq_mat)
    n_samples <- ncol(freq_mat)

    ## A trajectory needs a time order: text without one is refused, and a
    ## sample without a timepoint is left out
    no_time <- is.na(.timeOrder(timepoints, timeCol))
    if (any(no_time)) {
        warning(sum(no_time), " sample(s) with NA in '", timeCol,
                "' are left out.")
        freq_mat[, no_time] <- NA
    }

    if (usePassedOnly) freq_mat[!passed] <- NA
    entries <- which(!is.na(freq_mat), arr.ind = TRUE)

    if (nrow(entries) == 0L) {
        return(DataFrame(
            host_id = character(),
            variant_key = character(),
            timepoint = numeric(),
            frequency = numeric(),
            qc_passed = logical()))
    }

    DataFrame(
        host_id = host_ids[entries[, "col"]],
        variant_key = variant_keys[entries[, "row"]],
        timepoint = timepoints[entries[, "col"]],
        frequency = freq_mat[entries],
        qc_passed = passed[entries])
}


# ============================================================
# detectEmergingVariants -- find frequency-crossing variants
# ============================================================

#' Detect emerging variants across time points
#'
#' Identifies variants whose allele frequency crosses a threshold
#' between consecutive time points within each host. A variant is
#' considered emerging if it was below \code{freqThreshold} at one
#' time point and increased by at least \code{minIncrease} at the
#' next. Replicate samples taken at one time point are combined first:
#' the variant's frequency there is the mean of the replicates' calls,
#' so replicates never form an interval of their own and the result
#' does not depend on the order of the columns. Intervals join the
#' consecutive time points at which the variant has a call in the
#' host; a time point without a call for it is passed over, since no
#' call is not a frequency of 0, so a variant first called at a later
#' time point is not reported.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with longitudinal
#'   samples.
#' @param hostCol Character scalar. Column name in \code{colData}
#'   for host identifier (default: \code{"host_id"}). Samples whose
#'   host is \code{NA} are left out, with a warning.
#' @param timeCol Character scalar. Column name in \code{colData}
#'   for time point (default: \code{"timepoint"}), read as in
#'   \code{\link{trackFrequency}}: it is put in time order, and samples
#'   whose time point is \code{NA} are left out, with a warning.
#' @param freqThreshold Numeric scalar. Frequency threshold for
#'   emergence (default: 0.05).
#' @param minIncrease Numeric scalar. Minimum frequency increase
#'   between consecutive time points (default: 0.05).
#' @param usePassedOnly Logical. If TRUE (default), an interval is
#'   reported only when the later time point has a call that passed QC,
#'   and its frequency is the mean of those calls. The earlier time
#'   point is the mean of its passing calls too, or of all its calls
#'   when none passed: \code{\link{flagISNV}} fails a frequency below
#'   the filter's \code{minFreq}, and such a frequency is the low start
#'   that emergence is measured from. \code{qc_from} says which.
#'   Intervals are formed on every call before QC is applied, so a time
#'   point that failed QC, for example above \code{maxFreq}, is never
#'   passed over. With \code{FALSE} every call is used.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{host_id}{Host identifier.}
#'     \item{variant_key}{Variant identifier (chrom:pos:alt).}
#'     \item{time_from}{Earlier time point.}
#'     \item{time_to}{Later time point.}
#'     \item{freq_from}{Frequency at the earlier time point.}
#'     \item{freq_to}{Frequency at the later time point.}
#'     \item{delta_freq}{Change in frequency.}
#'     \item{qc_from}{Logical. Whether a call at the earlier time point
#'       passed QC (\code{TRUE} without a \code{qcPass} assay). A
#'       \code{FALSE} start failed QC, which a frequency below
#'       \code{minFreq} does by design; one that failed for low depth
#'       may be set aside.}
#'   }
#'   One row is one interval, between two distinct time points. Rows
#'   are ordered by host, then by variant, each in order of first
#'   appearance in the object, then by time point.
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C")
#' S4Vectors::mcols(gr)$alt <- c("T", "G")
#' freq <- matrix(c(0.02, 0.10, 0.12, 0.15), nrow = 2, ncol = 2)
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freq),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("S1_t1", "S1_t2"),
#'         host_id = c("H1", "H1"),
#'         timepoint = c(1, 2)))
#' detectEmergingVariants(whe)
detectEmergingVariants <- function(whe,
                                   hostCol = "host_id",
                                   timeCol = "timepoint",
                                   freqThreshold = 0.05,
                                   minIncrease = 0.05,
                                   usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!is.numeric(freqThreshold) || length(freqThreshold) != 1L ||
        freqThreshold < 0 || freqThreshold > 1)
        stop("'freqThreshold' must be a single number between 0 and 1.")
    if (!is.numeric(minIncrease) || length(minIncrease) != 1L ||
        minIncrease < 0)
        stop("'minIncrease' must be a single non-negative number.")

    if (!isTRUE(usePassedOnly) && !isFALSE(usePassedOnly))
        stop("'usePassedOnly' must be TRUE or FALSE.")

    ## Intervals are formed on every call, so that a timepoint that failed
    ## QC is never passed over; QC is applied to the end points below
    traj <- trackFrequency(whe, hostCol = hostCol,
                           timeCol = timeCol,
                           usePassedOnly = FALSE)
    empty <- DataFrame(
        host_id = character(),
        variant_key = character(),
        time_from = numeric(),
        time_to = numeric(),
        freq_from = numeric(),
        freq_to = numeric(),
        delta_freq = numeric(),
        qc_from = logical())

    host <- as.character(traj$host_id)
    keep <- !is.na(host)
    if (!all(keep))
        warning(sum(is.na(colData(whe)[[hostCol]])), " sample(s) with ",
                "NA in '", hostCol, "' are left out.")
    if (!any(keep)) return(empty)
    host <- host[keep]
    key <- as.character(traj$variant_key)[keep]
    freq <- traj$frequency[keep]
    time <- traj$timepoint[keep]
    pass_freq <- ifelse(traj$qc_passed[keep], freq, NA_real_)

    ## One point per host, variant and timepoint: replicate samples taken
    ## at one timepoint are combined by the mean of their calls, as
    ## flagTemporalInconsistency() counts them as one timepoint, so that
    ## replicates never form an interval and column order does not matter
    host_rank <- match(host, unique(host))
    pair <- paste(host_rank, key, sep = ":")
    pair_rank <- match(pair, unique(pair))
    time_num <- .timeOrder(time, timeCol)
    point <- paste(pair_rank, time_num)
    point <- factor(point, levels = unique(point))
    first <- !duplicated(point)
    ## the mean of all the calls at a point, and of those that passed QC
    freq <- unname(vapply(split(freq, point), mean, numeric(1)))
    pass_freq <- unname(vapply(split(pass_freq, point), function(f)
        if (all(is.na(f))) NA_real_ else mean(f, na.rm = TRUE), numeric(1)))
    host_rank <- host_rank[first]
    pair_rank <- pair_rank[first]
    time_num <- time_num[first]
    host <- host[first]
    key <- key[first]
    time <- time[first]

    ## Order once: hosts and, within each host, variants in order of
    ## first appearance, then by timepoint in time order, as
    ## flagTemporalInconsistency() reads it
    o <- order(host_rank, pair_rank, time_num)
    host <- host[o]
    key <- key[o]
    freq <- freq[o]
    pass_freq <- pass_freq[o]
    time <- time[o]
    has_pass <- !is.na(pass_freq)

    ## Each point and the next form an interval when they share host and
    ## variant. With usePassedOnly the later point must have a call that
    ## passed QC, and is the mean of those; the earlier point is the mean
    ## of its passing calls too, or of all of them when none passed, since
    ## a frequency below minFreq, which flagISNV() fails, is the low start
    ## that emergence is measured from
    n <- length(o)
    i <- which(host[-n] == host[-1L] & key[-n] == key[-1L])
    if (usePassedOnly) {
        f_from <- ifelse(has_pass[i], pass_freq[i], freq[i])
        f_to <- pass_freq[i + 1L]
    } else {
        f_from <- freq[i]
        f_to <- freq[i + 1L]
    }
    increase <- f_to - f_from
    ok <- !is.na(f_to) & f_from < freqThreshold & increase >= minIncrease
    hit <- i[ok]
    if (!length(hit)) return(empty)

    DataFrame(
        host_id = host[hit],
        variant_key = key[hit],
        time_from = time[hit],
        time_to = time[hit + 1L],
        freq_from = f_from[ok],
        freq_to = f_to[ok],
        delta_freq = increase[ok],
        qc_from = has_pass[hit])
}


# ============================================================
# plotFrequencyTrajectory -- spaghetti plot
# ============================================================

#' Plot allele frequency trajectories over time
#'
#' Creates a ggplot2 spaghetti plot of allele frequency over time
#' for each variant within a specified host. Each line represents
#' one variant, coloured with the Okabe-Ito palette (Wong 2011);
#' with more than eight variants the colours repeat.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with longitudinal
#'   samples.
#' @param hostId Character scalar. Host identifier to plot.
#' @param hostCol Character scalar. Column name in \code{colData}
#'   for host identifier (default: \code{"host_id"}).
#' @param timeCol Character scalar. Column name in \code{colData}
#'   for time point (default: \code{"timepoint"}), read as in
#'   \code{\link{trackFrequency}}. Numbers given as text are drawn in
#'   numeric order.
#' @param usePassedOnly Logical. If TRUE (default), only use
#'   QC-passed variants.
#'
#' @return A ggplot object.
#'
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     gr <- GenomicRanges::GRanges("seg1",
#'         IRanges::IRanges(c(100, 200), width = 1))
#'     S4Vectors::mcols(gr)$ref <- c("A", "C")
#'     S4Vectors::mcols(gr)$alt <- c("T", "G")
#'     freq <- matrix(c(0.05, 0.10, 0.15, 0.20), nrow = 2, ncol = 2)
#'     whe <- WithinHostExperiment(
#'         assays = list(altFreq = freq),
#'         rowRanges = gr,
#'         colData = S4Vectors::DataFrame(
#'             sample_id = c("S1_t1", "S1_t2"),
#'             host_id = c("H1", "H1"),
#'             timepoint = c(1, 2)))
#'     plotFrequencyTrajectory(whe, hostId = "H1")
#' }
plotFrequencyTrajectory <- function(whe, hostId,
                                    hostCol = "host_id",
                                    timeCol = "timepoint",
                                    usePassedOnly = TRUE) {
    .requireGgplot2()

    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (missing(hostId) || !is.character(hostId) || length(hostId) != 1L)
        stop("'hostId' must be a single character string.")

    traj <- trackFrequency(whe, hostCol = hostCol,
                           timeCol = timeCol,
                           usePassedOnly = usePassedOnly)

    host_traj <- traj[as.character(traj$host_id) == hostId, ]
    if (nrow(host_traj) == 0L) {
        message("No trajectory data for host '", hostId, "'")
        return(ggplot2::ggplot())
    }

    plot_df <- as.data.frame(host_traj)
    ## numbers given as text would be drawn in text order, "10" before "2"
    if (is.character(plot_df$timepoint)) {
        tp <- plot_df$timepoint
        plot_df$timepoint <- factor(tp, levels = unique(
            tp[order(.timeOrder(tp, timeCol))]))
    }

    ggplot2::ggplot(plot_df,
        ggplot2::aes(
            x = .data$timepoint,
            y = .data$frequency,
            color = .data$variant_key,
            group = .data$variant_key)) +
        ggplot2::geom_line(linewidth = 0.6, alpha = 0.7) +
        ggplot2::geom_point(size = 1.5, alpha = 0.9) +
        ggplot2::scale_colour_manual(
            values = rep(unlist(.whe_pal, use.names = FALSE),
                         length.out = length(unique(plot_df$variant_key)))) +
        ggplot2::scale_y_continuous(
            name = "Alternative allele frequency",
            limits = c(0, NA)) +
        ggplot2::labs(
            x = "Time point",
            color = "Variant") +
        ggplot2::ggtitle(paste("Host:", hostId)) +
        .theme_pub(base_size = 9) +
        ggplot2::theme(
            legend.position = "right",
            legend.title = ggplot2::element_text(
                size = 8, face = "bold"))
}
