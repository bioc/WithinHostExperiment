# R/visualization.R
# Publication-quality visualisation (Nature/Science style)

#' @include methods.R
#' @include bridge.R
#' @importFrom SummarizedExperiment assay assay<- assayNames colData rowRanges
#' @importFrom GenomicRanges start
#' @importFrom S4Vectors DataFrame
NULL

## Avoid R CMD check notes for ggplot2 .data pronoun
utils::globalVariables(".data")

#' @keywords internal
.requireGgplot2 <- function() {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Install with: install.packages('ggplot2')")
    }
}

## Colorblind-safe palette: Wong B (2011) Nat Methods 8:441
## Exact RGB values from the paper
.whe_pal <- list(
    blue      = "#0072B2",
    vermillon = "#D55E00",
    green     = "#009E73",
    skyblue   = "#56B4E9",
    orange    = "#E69F00",
    purple    = "#CC79A7",
    black     = "#000000",
    gray      = "#999999"
)

## An inside legend is placed with legend.position.inside from ggplot2 3.5.0;
## earlier versions take the coordinates in legend.position.
#' @keywords internal
.legend_inside <- function(x, y) {
    if (utils::packageVersion("ggplot2") >= "3.5.0") {
        ggplot2::theme(legend.position = "inside",
                       legend.position.inside = c(x, y))
    } else {
        ggplot2::theme(legend.position = c(x, y))
    }
}

## Nature/Science base theme -- Arial/Helvetica, no gridlines,
## white background, minimal chrome
#' @keywords internal
.theme_pub <- function(base_size = 8) {
    ggplot2::theme_classic(base_size = base_size,
                           base_family = "sans") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = base_size + 1, face = "bold",
                hjust = 0, margin = ggplot2::margin(b = 4)),
            plot.subtitle = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(
                size = base_size, color = "black"),
            axis.text = ggplot2::element_text(
                size = base_size - 1, color = "black"),
            axis.line = ggplot2::element_line(
                linewidth = 0.4, color = "black"),
            axis.ticks = ggplot2::element_line(
                linewidth = 0.3, color = "black"),
            axis.ticks.length = ggplot2::unit(1.5, "pt"),
            legend.title = ggplot2::element_text(
                size = base_size - 1, face = "bold"),
            legend.text = ggplot2::element_text(
                size = base_size - 1),
            legend.key.size = ggplot2::unit(8, "pt"),
            legend.background = ggplot2::element_blank(),
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            panel.background = ggplot2::element_rect(
                fill = "white", color = NA),
            panel.grid = ggplot2::element_blank(),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(
                size = base_size, face = "bold"),
            plot.background = ggplot2::element_rect(
                fill = "white", color = NA),
            plot.margin = ggplot2::margin(6, 8, 6, 6)
        )
}

#' Plot iSNV frequency spectrum
#'
#' Creates a publication-quality histogram or density plot of
#' alternative allele frequencies, optionally faceted by sample
#' or metadata group. Only calls are drawn: a recorded frequency of 0
#' is a site assessed without the allele and is left out. A site at
#' frequency 1 is drawn unless QC has removed it.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param type Character. \code{"histogram"} (default) or
#'   \code{"density"}.
#' @param facetBy Character (optional). Column name in
#'   \code{colData} to facet by.
#' @param usePassedOnly Logical. Only show QC-passed variants
#'   (default: TRUE).
#' @param threshold Numeric scalar (optional). Variant calling
#'   threshold to display as reference line. Default \code{NULL}
#'   (no line shown). Set to e.g. 0.03 for a 3\% threshold.
#' @param freqRange Numeric vector of length 2 (optional). Frequency
#'   range for the x-axis as proportions, e.g. \code{c(0, 0.5)}.
#'   Default \code{NULL} uses the data range.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     vcf <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(vcf,
#'         colData = S4Vectors::DataFrame(sample_id = "d"),
#'         caller = "ivar")
#'     plotFrequencySpectrum(whe, threshold = 0.03)
#' }
plotFrequencySpectrum <- function(whe,
                                  type = c("histogram", "density"),
                                  facetBy = NULL,
                                  usePassedOnly = TRUE,
                                  threshold = NULL,
                                  freqRange = NULL) {
    .requireGgplot2()
    type <- match.arg(type)

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

    df_list <- lapply(seq_len(ncol(freq_mat)), function(j) {
        freq <- freq_mat[, j]
        if (!is.null(qc_mat)) freq[!.passedQC(qc_mat[, j])] <- NA
        ## a recorded 0 is a site assessed without the allele, not an iSNV
        valid <- .carries(freq) & is.finite(freq)
        if (sum(valid) == 0L) return(NULL)
        data.frame(sample_id = sample_ids[j],
                   altFreq = freq[valid],
                   stringsAsFactors = FALSE)
    })
    plot_df <- do.call(rbind, Filter(Negate(is.null), df_list))

    if (is.null(plot_df) || nrow(plot_df) == 0L) {
        message("No data to plot after QC filtering")
        return(ggplot2::ggplot())
    }

    n_variants <- nrow(plot_df)
    n_samples <- length(unique(plot_df$sample_id))

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$altFreq))

    if (type == "histogram") {
        if (n_samples > 1 && n_samples <= 7) {
            p <- p +
                ggplot2::geom_histogram(
                    ggplot2::aes(fill = .data$sample_id),
                    bins = 20, color = "white",
                    linewidth = 0.2, alpha = 0.55,
                    position = "identity") +
                ggplot2::scale_fill_manual(
                    values = c(.whe_pal$blue,
                               .whe_pal$vermillon,
                               .whe_pal$orange,
                               .whe_pal$green,
                               .whe_pal$purple,
                               .whe_pal$skyblue,
                               .whe_pal$gray)[
                                   seq_len(n_samples)],
                    name = NULL)
        } else if (n_samples > 7) {
            ## Too many samples for individual colours
            p <- p +
                ggplot2::geom_histogram(
                    bins = 25, fill = .whe_pal$blue,
                    color = "white", linewidth = 0.2,
                    alpha = 0.8)
        } else {
            p <- p +
                ggplot2::geom_histogram(
                    bins = 20, fill = .whe_pal$blue,
                    color = "white", linewidth = 0.2,
                    alpha = 0.8)
        }
        p <- p + ggplot2::labs(y = "Number of iSNVs")
    } else {
        p <- p +
            ggplot2::geom_density(
                fill = .whe_pal$skyblue, alpha = 0.4,
                color = .whe_pal$blue, linewidth = 0.5) +
            ggplot2::labs(y = "Density")
    }

    ## Optional threshold reference line
    if (!is.null(threshold) && is.numeric(threshold) &&
        length(threshold) == 1L) {
        p <- p +
            ggplot2::geom_vline(
                xintercept = threshold, linetype = "dashed",
                color = .whe_pal$gray, linewidth = 0.4) +
            ggplot2::annotate(
                "text", x = threshold + 0.02, y = Inf,
                vjust = 1.5,
                label = paste0(threshold * 100, "%"),
                color = .whe_pal$gray,
                size = 2.5, fontface = "plain")
    }

    ## Determine x-axis range from data or user input
    all_freq <- plot_df$altFreq
    if (!is.null(freqRange) && length(freqRange) == 2L) {
        x_lim <- freqRange
    } else {
        x_lim <- c(max(0, min(all_freq) - 0.02),
                    min(1, max(all_freq) + 0.02))
    }
    x_breaks <- pretty(x_lim, n = 5)

    p <- p +
        ggplot2::scale_x_continuous(
            name = "Alternative allele frequency (%)",
            breaks = x_breaks,
            labels = function(x) as.integer(x * 100)) +
        ## Clip the view rather than dropping data, which would warn
        ggplot2::coord_cartesian(xlim = x_lim) +
        ggplot2::ggtitle(NULL) +
        .theme_pub(base_size = 9) +
        .legend_inside(0.85, 0.90) +
        ggplot2::theme(legend.justification = c(0.5, 1))

    if (!is.null(facetBy) &&
        facetBy %in% colnames(colData(whe))) {
        cd <- colData(whe)
        facet_map <- stats::setNames(
            as.character(cd[[facetBy]]),
            as.character(cd$sample_id))
        plot_df[[facetBy]] <- facet_map[plot_df$sample_id]
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = .data$altFreq)) +
            ggplot2::geom_histogram(
                bins = 20, fill = .whe_pal$blue,
                color = "white", linewidth = 0.2,
                alpha = 0.8) +
            ggplot2::scale_x_continuous(
                name = "Alternative allele frequency (%)",
                labels = function(x) as.integer(x * 100)) +
            ggplot2::coord_cartesian(xlim = x_lim) +
            ggplot2::labs(y = "Number of iSNVs") +
            ggplot2::facet_wrap(
                stats::as.formula(paste("~", facetBy)),
                scales = "free_y") +
            .theme_pub(base_size = 9)
    }

    p
}


#' Plot donor-recipient frequency scatter
#'
#' Creates a publication-quality scatter plot comparing variant
#' frequencies between a donor and recipient, highlighting shared
#' and unique variants with distinct shapes and colours.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with
#'   transmission pairs.
#' @param pairId Character. Which pair to plot.
#' @param threshold Numeric. Variant calling threshold shown as
#'   reference lines (default: 0.03).
#'
#' @return A ggplot object.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     vcf1 <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     vcf2 <- system.file("extdata", "test_recipient.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(c(vcf1, vcf2),
#'         colData = S4Vectors::DataFrame(
#'             sample_id = c("donor", "recipient"),
#'             role = c("donor", "recipient"),
#'             pair_id = c("p1", "p1")),
#'         caller = "ivar")
#'     plotPairScatter(whe, pairId = "p1")
#' }
plotPairScatter <- function(whe, pairId, threshold = 0.03) {
    .requireGgplot2()

    pf <- exportPairFrequencies(whe, pairId = pairId)
    if (nrow(pf) == 0L) {
        message("Pair '", pairId, "' has no site where either host ",
                "carries an allele; nothing to plot")
        return(ggplot2::ggplot())
    }

    d_freq <- ifelse(is.na(pf$donor_freq), 0, pf$donor_freq)
    r_freq <- ifelse(is.na(pf$recipient_freq), 0,
                     pf$recipient_freq)
    ## a host carries the allele where its frequency is known and above
    ## 0, as in exportPairFrequencies(); every row it returns is carried
    ## by one host at least, so "Other" is a fallback that should not occur
    d_pres <- .carries(pf$donor_freq)
    r_pres <- .carries(pf$recipient_freq)
    category <- ifelse(
        pf$shared, "Shared",
        ifelse(d_pres & !r_pres, "Donor only",
               ifelse(r_pres & !d_pres, "Recipient only", "Other")))

    plot_df <- data.frame(
        donor = d_freq, recipient = r_freq,
        category = factor(
            category,
            levels = c("Shared", "Donor only",
                       "Recipient only", "Other")),
        stringsAsFactors = FALSE
    )

    n_shared <- sum(plot_df$category == "Shared")
    n_donor <- sum(plot_df$category == "Donor only")
    n_recip <- sum(plot_df$category == "Recipient only")

    cat_colors <- c(
        "Shared"         = .whe_pal$green,
        "Donor only"     = .whe_pal$blue,
        "Recipient only" = .whe_pal$vermillon,
        "Other"          = .whe_pal$gray
    )
    cat_shapes <- c(
        "Shared" = 16, "Donor only" = 17,
        "Recipient only" = 15, "Other" = 1
    )

    ## Build legend labels with counts
    cat_labels <- c(
        "Shared"         = paste0("Shared (n = ", n_shared, ")"),
        "Donor only"     = paste0("Donor only (n = ",
                                  n_donor, ")"),
        "Recipient only" = paste0("Recipient only (n = ",
                                  n_recip, ")"),
        "Other"          = "Other"
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(
        x = .data$donor, y = .data$recipient,
        color = .data$category,
        shape = .data$category)) +
        ## Identity line
        ggplot2::geom_abline(
            slope = 1, intercept = 0,
            linetype = "dashed", color = .whe_pal$gray,
            linewidth = 0.3) +
        ## Threshold reference lines
        ggplot2::geom_vline(
            xintercept = threshold, linetype = "dotted",
            color = "#CCCCCC", linewidth = 0.3) +
        ggplot2::geom_hline(
            yintercept = threshold, linetype = "dotted",
            color = "#CCCCCC", linewidth = 0.3) +
        ## Points: non-shared with jitter, shared on top
        ggplot2::geom_point(
            data = plot_df[plot_df$category != "Shared", ],
            size = 2, alpha = 0.7, stroke = 0.3,
            position = ggplot2::position_jitter(
                width = 0.005, height = 0.005,
                seed = 42)) +
        ggplot2::geom_point(
            data = plot_df[plot_df$category == "Shared", ],
            size = 2.5, alpha = 0.9, stroke = 0.3) +
        ## Scales
        ggplot2::scale_color_manual(
            values = cat_colors, labels = cat_labels) +
        ggplot2::scale_shape_manual(
            values = cat_shapes, labels = cat_labels) +
        ggplot2::scale_x_continuous(
            name = "Donor frequency (%)",
            labels = function(x) as.integer(x * 100)) +
        ggplot2::scale_y_continuous(
            name = "Recipient frequency (%)",
            labels = function(x) as.integer(x * 100)) +
        ## Clip the view rather than dropping jittered points
        ggplot2::coord_equal(
            xlim = c(0, max(c(d_freq, r_freq)) * 1.1 + 0.02),
            ylim = c(0, max(c(d_freq, r_freq)) * 1.1 + 0.02)) +
        ggplot2::ggtitle(NULL) +
        .theme_pub(base_size = 9) +
        .legend_inside(0.98, 0.98) +
        ggplot2::theme(
            legend.justification = c(1, 1),
            legend.title = ggplot2::element_blank(),
            legend.key.size = ggplot2::unit(10, "pt"),
            legend.spacing.y = ggplot2::unit(1, "pt"))
}


#' QC dashboard plot
#'
#' Creates a publication-quality four-panel dashboard showing
#' frequency distribution, read depth, QC filtering waterfall,
#' and per-sample pass/fail summary. Panels a and d draw and count
#' calls, frequencies above 0, as \code{\link{qcSummary}} does: a
#' recorded 0 is a site assessed without the allele.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#'
#' @return A patchwork composite ggplot object.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'     vcf <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(vcf,
#'         colData = S4Vectors::DataFrame(sample_id = "d"),
#'         caller = "ivar")
#'     whe <- flagISNV(whe, ISNVFilter())
#'     plotQCDashboard(whe)
#' }
plotQCDashboard <- function(whe) {
    .requireGgplot2()
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required. ",
             "Install with: install.packages('patchwork')")
    }

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if ("qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    depth_mat <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

    bs <- 8

    ## A call is a frequency above 0; a recorded 0 is a site assessed
    ## without the allele and is neither drawn nor counted as an iSNV
    is_call <- .carries(freq_mat)

    ## ---- Panel A: Frequency distribution by sample ----
    ## a sample without calls has no bars
    freq_df <- do.call(rbind, Filter(Negate(is.null),
        lapply(seq_len(ncol(freq_mat)), function(j) {
            vals <- freq_mat[is_call[, j], j]
            if (!length(vals)) return(NULL)
            data.frame(sample = sample_ids[j], freq = vals,
                       stringsAsFactors = FALSE)
        })))

    n_samp <- length(unique(freq_df$sample))

    if (is.null(freq_df)) {
        pA <- ggplot2::ggplot() +
            ggplot2::labs(title = expression(bold("a"))) +
            ggplot2::theme_void()
    } else if (n_samp <= 7) {
        all_pal <- c(.whe_pal$blue, .whe_pal$vermillon,
                     .whe_pal$orange, .whe_pal$green,
                     .whe_pal$purple, .whe_pal$skyblue,
                     .whe_pal$gray)
        pA <- ggplot2::ggplot(freq_df,
            ggplot2::aes(x = .data$freq,
                         fill = .data$sample)) +
            ggplot2::geom_histogram(
                bins = 20, color = "white",
                linewidth = 0.15, alpha = 0.55,
                position = "identity") +
            ggplot2::scale_fill_manual(
                values = all_pal[seq_len(n_samp)],
                name = NULL)
    } else {
        pA <- ggplot2::ggplot(freq_df,
            ggplot2::aes(x = .data$freq)) +
            ggplot2::geom_histogram(
                bins = 25, fill = .whe_pal$blue,
                color = "white", linewidth = 0.15,
                alpha = 0.8)
    }

    if (!is.null(freq_df)) {
        all_freq_vals <- freq_df$freq
        x_lim_d <- c(max(0, min(all_freq_vals) - 0.02),
                     min(1, max(all_freq_vals) + 0.02))
        pA <- pA +
            ggplot2::scale_x_continuous(
                name = "Alternative allele frequency (%)",
                labels = function(x) as.integer(x * 100)) +
            ggplot2::coord_cartesian(xlim = x_lim_d) +
            ggplot2::labs(y = "Number of iSNVs",
                          title = expression(bold("a"))) +
            .theme_pub(base_size = bs) +
            .legend_inside(0.82, 0.88) +
            ggplot2::theme(legend.key.size = ggplot2::unit(7, "pt"))
    }

    ## ---- Panel B: Depth distribution ----
    if (!is.null(depth_mat)) {
        all_depth <- as.numeric(depth_mat[!is.na(depth_mat)])
        med_d <- stats::median(all_depth)
        pB <- ggplot2::ggplot(
            data.frame(d = all_depth),
            ggplot2::aes(x = .data$d)) +
            ggplot2::geom_histogram(
                bins = 20, fill = .whe_pal$orange,
                color = "white", linewidth = 0.15,
                alpha = 0.8) +
            ggplot2::geom_vline(
                xintercept = med_d,
                linetype = "dashed", color = "black",
                linewidth = 0.3) +
            ggplot2::annotate(
                "text", x = med_d, y = Inf,
                vjust = 1.5, hjust = -0.08,
                label = paste0("median = ",
                    format(round(med_d), big.mark = ",")),
                size = 2.2, color = "black") +
            ggplot2::labs(
                x = "Read depth",
                y = "Count",
                title = expression(bold("b"))) +
            .theme_pub(base_size = bs)
    } else {
        pB <- ggplot2::ggplot() +
            ggplot2::labs(title = expression(bold("b"))) +
            ggplot2::theme_void()
    }

    ## ---- Panel C: QC waterfall ----
    log_df <- as.data.frame(qcLog(whe))
    if (nrow(log_df) > 0) {
        first_total <- log_df$n_passed[1] + log_df$n_flagged[1]

        ## Human-readable labels with parameter values
        step_labels <- vapply(seq_len(nrow(log_df)), function(i) {
            p <- log_df$parameter[i]
            v <- log_df$value[i]
            switch(p,
                minDepth = paste0("Depth \u2265 ",
                    format(as.numeric(v), big.mark = ",")),
                minFreq = paste0("Freq \u2265 ",
                    as.numeric(v) * 100, "%"),
                maxFreq = paste0("Freq \u2264 ",
                    as.numeric(v) * 100, "%"),
                minAltReads = paste0("Alt reads \u2265 ", v),
                maxStrandBias = paste0("Strand bias \u2264 ",
                    v),
                paste0(p, "=", v))
        }, character(1))

        steps <- data.frame(
            step = c("Input", step_labels, "Passed"),
            count = c(first_total, log_df$n_passed,
                      log_df$n_passed[nrow(log_df)]),
            stringsAsFactors = FALSE
        )
        steps <- steps[!duplicated(steps$step,
                                    fromLast = TRUE), ]
        steps$step <- factor(steps$step, levels = steps$step)

        ## Colour: input=skyblue, filters=orange, passed=green
        bar_fills <- c(
            .whe_pal$skyblue,
            rep(.whe_pal$orange, nrow(steps) - 2),
            .whe_pal$green)

        pC <- ggplot2::ggplot(steps,
            ggplot2::aes(x = .data$step,
                         y = .data$count)) +
            ggplot2::geom_col(
                fill = bar_fills, width = 0.6) +
            ggplot2::geom_text(
                ggplot2::aes(label = .data$count),
                vjust = -0.4, size = 2.2,
                color = "black") +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(
                    mult = c(0, 0.12))) +
            ggplot2::labs(
                x = NULL,
                y = "Variant\u00d7sample entries",
                title = expression(bold("c"))) +
            .theme_pub(base_size = bs) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    angle = 35, hjust = 1, size = 6.5))
    } else {
        pC <- ggplot2::ggplot() +
            ggplot2::labs(title = expression(bold("c"))) +
            ggplot2::theme_void()
    }

    ## ---- Panel D: Per-sample pass/flagged ----
    n_total <- vapply(seq_len(ncol(freq_mat)), function(j) {
        sum(is_call[, j])
    }, integer(1))
    n_pass <- if (!is.null(qc_mat)) {
        vapply(seq_len(ncol(qc_mat)), function(j) {
            sum(.passedQC(qc_mat[, j]) & is_call[, j])
        }, integer(1))
    } else {
        n_total
    }

    n_samp_d <- length(sample_ids)

    if (n_samp_d <= 15) {
        ## Few samples: show individual bars
        bar_df <- data.frame(
            sample = rep(sample_ids, 2),
            status = factor(
                rep(c("Passed", "Flagged"),
                    each = n_samp_d),
                levels = c("Passed", "Flagged")),
            count = c(n_pass, n_total - n_pass),
            stringsAsFactors = FALSE)

        pD <- ggplot2::ggplot(bar_df,
            ggplot2::aes(x = .data$sample,
                         y = .data$count,
                         fill = .data$status)) +
            ggplot2::geom_col(width = 0.55) +
            ggplot2::geom_text(
                data = data.frame(
                    sample = sample_ids,
                    y = n_total,
                    label = paste0(n_pass, "/", n_total)),
                ggplot2::aes(x = .data$sample,
                             y = .data$y,
                             label = .data$label),
                inherit.aes = FALSE,
                vjust = -0.4, size = 2.2, color = "black") +
            ggplot2::scale_fill_manual(
                values = c(Passed = .whe_pal$green,
                           Flagged = .whe_pal$vermillon),
                name = NULL) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(
                    mult = c(0, 0.15))) +
            ggplot2::labs(x = NULL, y = "Number of iSNVs",
                          title = expression(bold("d"))) +
            .theme_pub(base_size = bs) +
            .legend_inside(0.85, 0.88) +
            ggplot2::theme(legend.key.size = ggplot2::unit(7, "pt"))
    } else {
        ## Many samples: histogram of iSNVs per sample
        site_df <- data.frame(n = n_total,
                              stringsAsFactors = FALSE)
        pass_rate <- sum(n_pass) / max(sum(n_total), 1)

        pD <- ggplot2::ggplot(site_df,
            ggplot2::aes(x = .data$n)) +
            ggplot2::geom_histogram(
                bins = max(max(n_total), 5),
                fill = .whe_pal$green,
                color = "white", linewidth = 0.15,
                alpha = 0.8) +
            ggplot2::annotate(
                "text", x = Inf, y = Inf,
                hjust = 1.1, vjust = 1.5,
                label = paste0(n_samp_d, " samples\n",
                    round(pass_rate * 100, 1), "% pass rate"),
                size = 2.2, color = "black") +
            ggplot2::labs(
                x = "iSNVs per sample",
                y = "Number of samples",
                title = expression(bold("d"))) +
            .theme_pub(base_size = bs)
    }

    ## ---- Compose ----
    (pA + pB) / (pC + pD)
}
