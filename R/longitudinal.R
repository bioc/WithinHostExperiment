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
#'   variants where \code{qcPass == TRUE}.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{host_id}{Host identifier.}
#'     \item{variant_key}{Unique variant identifier
#'       (chrom:pos:alt).}
#'     \item{timepoint}{Time point value.}
#'     \item{frequency}{Alternative allele frequency at that
#'       time point.}
#'     \item{qc_passed}{Logical. Whether the variant passed QC
#'       at that time point (\code{TRUE} if no QC assay present).}
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
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
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

    if (!is.null(qc_mat)) freq_mat[!.passedQC(qc_mat)] <- NA
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
        qc_passed = if (is.null(qc_mat)) {
            rep(TRUE, nrow(entries))
        } else {
            qc_mat[entries]
        })
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
#' next.
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
#' @param usePassedOnly Logical. If TRUE (default), only use
#'   QC-passed variants.
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
#'   }
#'   One row is one interval. Rows are ordered by host, then by
#'   variant, each in order of first appearance in the object, then
#'   by time point; samples tied at one time point stay in column
#'   order.
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

    traj <- trackFrequency(whe, hostCol = hostCol,
                           timeCol = timeCol,
                           usePassedOnly = usePassedOnly)
    empty <- DataFrame(
        host_id = character(),
        variant_key = character(),
        time_from = numeric(),
        time_to = numeric(),
        freq_from = numeric(),
        freq_to = numeric(),
        delta_freq = numeric())

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

    ## Order once: hosts and, within each host, variants in order of
    ## first appearance, then by timepoint in time order, as
    ## flagTemporalInconsistency() reads it. order() is stable, so
    ## samples tied at one timepoint keep their column order.
    host_rank <- match(host, unique(host))
    pair <- paste(host_rank, key, sep = ":")
    o <- order(host_rank, match(pair, unique(pair)),
               .timeOrder(time, timeCol))
    host <- host[o]
    key <- key[o]
    freq <- freq[o]
    time <- time[o]

    ## Each row and the next form an interval when they share host and
    ## variant
    n <- length(o)
    i <- which(host[-n] == host[-1L] & key[-n] == key[-1L])
    increase <- freq[i + 1L] - freq[i]
    hit <- i[freq[i] < freqThreshold & increase >= minIncrease]
    if (!length(hit)) return(empty)

    DataFrame(
        host_id = host[hit],
        variant_key = key[hit],
        time_from = time[hit],
        time_to = time[hit + 1L],
        freq_from = freq[hit],
        freq_to = freq[hit + 1L],
        delta_freq = freq[hit + 1L] - freq[hit])
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
