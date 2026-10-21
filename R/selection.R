# R/selection.R
# Within-host dN/dS and selection coefficients

#' @include methods.R annotate-codon.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
NULL

## Standard genetic code (NCBI translation table 1)
.CODON_TABLE <- c(
    TTT="F",TTC="F",TTA="L",TTG="L", CTT="L",CTC="L",CTA="L",CTG="L",
    ATT="I",ATC="I",ATA="I",ATG="M", GTT="V",GTC="V",GTA="V",GTG="V",
    TCT="S",TCC="S",TCA="S",TCG="S", CCT="P",CCC="P",CCA="P",CCG="P",
    ACT="T",ACC="T",ACA="T",ACG="T", GCT="A",GCC="A",GCA="A",GCG="A",
    TAT="Y",TAC="Y",TAA="*",TAG="*", CAT="H",CAC="H",CAA="Q",CAG="Q",
    AAT="N",AAC="N",AAA="K",AAG="K", GAT="D",GAC="D",GAA="E",GAG="E",
    TGT="C",TGC="C",TGA="*",TGG="W", CGT="R",CGC="R",CGA="R",CGG="R",
    AGT="S",AGC="S",AGA="R",AGG="R", GGT="G",GGC="G",GGA="G",GGG="G"
)

#' @keywords internal
.read_fasta_simple <- function(path) {
    lines <- readLines(path)
    header_idx <- grep("^>", lines)
    if (length(header_idx) == 0L)
        stop("No FASTA header lines found in: ", path)
    seqs <- list()
    for (i in seq_along(header_idx)) {
        name <- sub("^>\\s*", "", lines[header_idx[i]])
        name <- sub("\\s+.*", "", name)
        start <- header_idx[i] + 1L
        end <- if (i < length(header_idx)) {
            header_idx[i + 1L] - 1L
        } else {
            length(lines)
        }
        seq_lines <- if (start <= end) lines[start:end] else character()
        ## Whitespace inside a record would shift every later coordinate
        seqs[[name]] <- gsub("[[:space:]]", "",
                             paste(seq_lines, collapse = ""))
    }
    seqs
}

#' @keywords internal
.read_gff_cds <- function(gff, nameCol = "gene") {
    lines <- readLines(gff)
    lines <- lines[!startsWith(lines, "#") & nchar(lines) > 0]
    fields <- strsplit(lines, "\t")
    fields <- fields[vapply(fields, length, integer(1)) >= 9L]
    if (length(fields) == 0L)
        stop("No valid GFF lines found.")

    gff_df <- data.frame(
        seqid  = vapply(fields, `[`, character(1), 1L),
        type   = vapply(fields, `[`, character(1), 3L),
        start  = as.integer(vapply(fields, `[`, character(1), 4L)),
        end    = as.integer(vapply(fields, `[`, character(1), 5L)),
        strand = vapply(fields, `[`, character(1), 7L),
        phase  = vapply(fields, `[`, character(1), 8L),
        attrs  = vapply(fields, `[`, character(1), 9L),
        stringsAsFactors = FALSE
    )
    cds <- gff_df[gff_df$type == "CDS", , drop = FALSE]
    if (nrow(cds) == 0L)
        stop("No CDS features found in GFF file. Available types: ",
             paste(unique(gff_df$type), collapse = ", "))

    cds$phase <- as.integer(ifelse(cds$phase == ".", "0", cds$phase))

    parse_attr <- function(attr_str, key) {
        pairs <- strsplit(attr_str, ";")[[1L]]
        for (p in pairs) {
            kv <- strsplit(trimws(p), "=")[[1L]]
            if (length(kv) == 2L && kv[1L] == key) return(kv[2L])
        }
        NA_character_
    }
    cds$gene <- vapply(cds$attrs, parse_attr, character(1),
                       key = nameCol, USE.NAMES = FALSE)
    if (all(is.na(cds$gene)) && nameCol != "Name")
        cds$gene <- vapply(cds$attrs, parse_attr, character(1),
                           key = "Name", USE.NAMES = FALSE)
    if (all(is.na(cds$gene)))
        cds$gene <- vapply(cds$attrs, parse_attr, character(1),
                           key = "ID", USE.NAMES = FALSE)
    rownames(cds) <- NULL
    cds
}

## Nei and Gojobori (1986) synonymous sites of every sense codon under the
## standard genetic code: at each codon position, the fraction of the three
## possible substitutions that keep the amino acid. Substitutions to a stop
## codon count as nonsynonymous; the nonsynonymous sites are 3 minus this.
.NG_SYN_SITES <- local({
    bases <- c("A", "C", "G", "T")
    sense <- names(.CODON_TABLE)[.CODON_TABLE != "*"]
    vapply(sense, function(codon) {
        aa <- .CODON_TABLE[[codon]]
        s <- 0
        for (p in seq_len(3L)) {
            for (b in setdiff(bases, substr(codon, p, p))) {
                mutant <- codon
                substr(mutant, p, p) <- b
                if (.CODON_TABLE[[mutant]] == aa) s <- s + 1 / 3
            }
        }
        s
    }, numeric(1))
})

#' @keywords internal
.revcomp_vec <- function(x) {
    vapply(strsplit(x, "", fixed = TRUE), function(ch) {
        paste(rev(chartr("ACGTacgt", "TGCAtgca", ch)), collapse = "")
    }, character(1))
}

## Synonymous and nonsynonymous site totals per gene over the complete
## codons of its CDS features. Codons shared by overlapping CDS features of
## the same gene (for example ORF1a within ORF1ab) are counted once; stop
## codons and codons with ambiguous bases are skipped.
#' @keywords internal
.cds_site_counts <- function(cds, ref_seqs) {
    codon_rows <- lapply(seq_len(nrow(cds)), function(i) {
        if (is.na(cds$gene[i]) || !cds$seqid[i] %in% names(ref_seqs))
            return(NULL)
        usable <- cds$end[i] - cds$start[i] + 1L - cds$phase[i]
        n_codons <- usable %/% 3L
        if (n_codons < 1L) return(NULL)
        offsets <- 3L * (seq_len(n_codons) - 1L)
        low <- if (cds$strand[i] == "-") {
            cds$end[i] - cds$phase[i] - offsets - 2L
        } else {
            cds$start[i] + cds$phase[i] + offsets
        }
        data.frame(gene = cds$gene[i], seqid = cds$seqid[i],
                   strand = cds$strand[i], low = low,
                   stringsAsFactors = FALSE)
    })
    codons <- unique(do.call(rbind, codon_rows))
    if (is.null(codons) || nrow(codons) == 0L) {
        return(data.frame(gene = character(), S_sites = numeric(),
                          N_sites = numeric(), stringsAsFactors = FALSE))
    }
    codon_seq <- character(nrow(codons))
    for (sid in unique(codons$seqid)) {
        rows <- which(codons$seqid == sid)
        codon_seq[rows] <- toupper(substring(ref_seqs[[sid]],
            codons$low[rows], codons$low[rows] + 2L))
    }
    minus <- codons$strand == "-"
    codon_seq[minus] <- .revcomp_vec(codon_seq[minus])
    syn <- unname(.NG_SYN_SITES[codon_seq])
    keep <- !is.na(syn)
    s_by_gene <- tapply(syn[keep], codons$gene[keep], sum)
    n_by_gene <- tapply(3 - syn[keep], codons$gene[keep], sum)
    data.frame(gene = names(s_by_gene), S_sites = as.numeric(s_by_gene),
               N_sites = as.numeric(n_by_gene[names(s_by_gene)]),
               stringsAsFactors = FALSE)
}

# ============================================================
# dndsWithinHost -- dN/dS ratio per sample
# ============================================================

#' Within-host dN/dS ratio
#'
#' Estimates the ratio of nonsynonymous to synonymous divergence among
#' the iSNVs of each sample with the counting method of Nei and Gojobori
#' (1986). The numbers of synonymous (\eqn{S}) and nonsynonymous
#' (\eqn{N}) sites are counted over the complete codons of every CDS in
#' \code{gff}, using the reference sequence in \code{refFasta} and the
#' standard genetic code. With \eqn{S_d} synonymous and \eqn{N_d}
#' nonsynonymous iSNVs present in a sample,
#' \deqn{p_S = S_d / S, \quad p_N = N_d / N,}
#' each proportion is corrected for multiple hits with Jukes and Cantor
#' (1969), \eqn{d = -\frac{3}{4}\log(1 - \frac{4}{3}p)}, and
#' \eqn{d_N/d_S} is their ratio. Under neutrality the expected ratio is
#' 1; values below 1 suggest purifying selection and values above 1
#' diversifying selection, but within-host ratios computed from few
#' iSNVs are highly variable and should be read with their counts.
#'
#' Each iSNV present in a sample (non-missing, non-zero frequency, and
#' passing QC when \code{usePassedOnly = TRUE}) is counted once,
#' irrespective of its frequency. Amino acid changes come from
#' \code{\link{annotateCodonChange}}, which should be given the CDS
#' records of the same \code{gff} and the sequence in \code{refFasta}.
#' Variants outside a CDS, or without an amino acid annotation, are not
#' counted. Nonsense changes count as
#' nonsynonymous. Codons shared by overlapping CDS features of the same
#' gene are counted once.
#'
#' @param whe A \code{\link{WithinHostExperiment}} annotated with
#'   \code{\link{annotateCodonChange}}.
#' @param gff Character scalar. Path to the GFF3 file with CDS features
#'   used for the annotation.
#' @param refFasta Character scalar. Path to the reference genome FASTA
#'   used for the annotation.
#' @param gffFeatureCol Character scalar. Column name in
#'   \code{mcols(rowRanges)} containing the gene annotation
#'   (default: \code{"GENEID"}).
#' @param refAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the reference amino acid
#'   (default: \code{"REFAA"}).
#' @param altAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the alternative amino acid
#'   (default: \code{"VARAA"}).
#' @param usePassedOnly Logical. If \code{TRUE} (default), only use
#'   variants where \code{qcPass == TRUE}.
#' @param nameCol Character scalar. GFF3 attribute holding the gene name
#'   (default: \code{"gene"}). It must name the same labels that
#'   \code{\link{annotateCodonChange}} wrote into \code{GENEID}
#'   (its \code{geneCol}); otherwise no iSNV is counted and the
#'   function warns.
#'
#' @return A \code{DataFrame} with one row per sample and gene carrying
#'   at least one counted iSNV (one row with \code{gene = NA} when a
#'   sample has none):
#'   \describe{
#'     \item{sample_id}{Sample identifier.}
#'     \item{nS, nN}{Synonymous and nonsynonymous iSNVs in the sample.}
#'     \item{S_sites, N_sites}{Synonymous and nonsynonymous sites over
#'       all CDS features.}
#'     \item{pS, pN}{Proportions \eqn{S_d/S} and \eqn{N_d/N}.}
#'     \item{dNdS}{Jukes-Cantor corrected \eqn{d_N/d_S} (\code{NA} when
#'       \eqn{d_S} is zero or a proportion reaches 0.75).}
#'     \item{gene, gene_nS, gene_nN, gene_S_sites, gene_N_sites,
#'       gene_dNdS}{The same quantities for each gene.}
#'   }
#'
#' @references
#' Nei M, Gojobori T (1986). Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol Biol Evol} 3:418-426.
#' \doi{10.1093/oxfordjournals.molbev.a040410}
#'
#' Jukes TH, Cantor CR (1969). Evolution of protein molecules. In Munro
#' HN (ed.), \emph{Mammalian Protein Metabolism}, pp. 21-132. Academic
#' Press, New York.
#'
#' @seealso \code{\link{annotateCodonChange}}
#'
#' @export
#' @examples
#' ref <- tempfile(fileext = ".fa")
#' writeLines(c(">seg1", "ATGCTGAAAGGGTAA"), ref)
#' gff <- tempfile(fileext = ".gff3")
#' writeLines(c("##gff-version 3",
#'     "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(1, 6), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "G")
#' S4Vectors::mcols(gr)$alt <- c("G", "A")
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' if (requireNamespace("txdbmaker", quietly = TRUE)) {
#'     cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
#'         strand = "+", ID = "cds-A", gene = "geneA")
#'     whe <- annotateCodonChange(whe, cds,
#'         Biostrings::readDNAStringSet(ref))
#'     dndsWithinHost(whe, gff = gff, refFasta = ref)
#' }
dndsWithinHost <- function(whe, gff, refFasta,
                           gffFeatureCol = "GENEID",
                           refAACol = "REFAA",
                           altAACol = "VARAA",
                           usePassedOnly = TRUE,
                           nameCol = "gene") {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (missing(gff) || missing(refFasta))
        stop("'gff' and 'refFasta' are required: dN/dS is normalised by ",
             "the synonymous and nonsynonymous sites of the coding ",
             "sequences.")
    if (!is.character(gff) || length(gff) != 1L || !file.exists(gff))
        stop("GFF file not found: ", gff)
    if (!is.character(refFasta) || length(refFasta) != 1L ||
        !file.exists(refFasta))
        stop("Reference FASTA not found: ", refFasta)

    rr <- rowRanges(whe)
    mc <- mcols(rr)
    for (column in c(refAACol, altAACol, gffFeatureCol)) {
        if (!column %in% colnames(mc))
            stop("Column '", column, "' not found in mcols(rowRanges). ",
                 "Run annotateCodonChange() first. Available columns: ",
                 paste(colnames(mc), collapse = ", "))
        if (is(mc[[column]], "List"))
            stop("Column '", column, "' holds more than one value per ",
                 "site. dndsWithinHost() needs one transcript per site: ",
                 "annotate with collapse = \"longest\", the default of ",
                 "annotateCodonChange().")
    }

    sites <- .cds_site_counts(.read_gff_cds(gff, nameCol),
                              .read_fasta_simple(refFasta))
    if (nrow(sites) == 0L)
        stop("No complete codons found for the CDS features in 'gff' ",
             "on the sequences in 'refFasta'.")
    S_total <- sum(sites$S_sites)
    N_total <- sum(sites$N_sites)

    ref_aa <- as.character(mc[[refAACol]])
    alt_aa <- as.character(mc[[altAACol]])
    gene_vals <- as.character(mc[[gffFeatureCol]])
    annotated <- !is.na(ref_aa) & !is.na(alt_aa) & !is.na(gene_vals) &
        nzchar(ref_aa) & nzchar(alt_aa)
    coding <- annotated & gene_vals %in% sites$gene
    ## The gene labels of the annotation and of the site counts come from
    ## two sources, so a naming mismatch would silently give every count
    ## as zero.
    if (any(annotated) && !any(coding))
        warning("None of the ", sum(annotated), " annotated site(s) ",
                "carries a gene label present in 'gff', so every count is ",
                "zero.
  labels in the annotation: ",
                paste(utils::head(unique(gene_vals[annotated]), 5L),
                      collapse = ", "),
                "
  labels in 'gff' (nameCol = \"", nameCol, "\"): ",
                paste(utils::head(sites$gene, 5L), collapse = ", "),
                call. = FALSE)
    is_syn <- coding & ref_aa == alt_aa
    is_non <- coding & ref_aa != alt_aa

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids))
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))

    jukes_cantor <- function(p) {
        ifelse(p < 0.75, -0.75 * log(1 - 4 * p / 3), NA_real_)
    }
    ratio <- function(nS, nN, S, N) {
        dS <- jukes_cantor(nS / S)
        dN <- jukes_cantor(nN / N)
        if (is.na(dS) || is.na(dN) || dS == 0) NA_real_ else dN / dS
    }

    results <- lapply(seq_len(ncol(freq_mat)), function(j) {
        freq_j <- freq_mat[, j]
        if (!is.null(qc_mat)) freq_j[!qc_mat[, j]] <- NA
        present <- !is.na(freq_j) & freq_j > 0

        nS <- sum(present & is_syn)
        nN <- sum(present & is_non)
        base_row <- data.frame(
            sample_id = sample_ids[j], nS = nS, nN = nN,
            S_sites = S_total, N_sites = N_total,
            pS = nS / S_total, pN = nN / N_total,
            dNdS = ratio(nS, nN, S_total, N_total),
            stringsAsFactors = FALSE)

        genes <- sort(unique(gene_vals[present & coding]))
        if (length(genes) == 0L) {
            return(cbind(base_row, data.frame(
                gene = NA_character_, gene_nS = NA_integer_,
                gene_nN = NA_integer_, gene_S_sites = NA_real_,
                gene_N_sites = NA_real_, gene_dNdS = NA_real_,
                stringsAsFactors = FALSE)))
        }
        gene_rows <- lapply(genes, function(g) {
            in_gene <- present & coding & gene_vals == g
            gS <- sites$S_sites[match(g, sites$gene)]
            gN <- sites$N_sites[match(g, sites$gene)]
            gnS <- sum(in_gene & is_syn)
            gnN <- sum(in_gene & is_non)
            data.frame(gene = g, gene_nS = gnS, gene_nN = gnN,
                       gene_S_sites = gS, gene_N_sites = gN,
                       gene_dNdS = ratio(gnS, gnN, gS, gN),
                       stringsAsFactors = FALSE)
        })
        cbind(base_row[rep(1L, length(genes)), , drop = FALSE],
              do.call(rbind, gene_rows))
    })

    out <- do.call(rbind, results)
    rownames(out) <- NULL
    DataFrame(out)
}


# ============================================================
# estimateSelectionCoefficient
# ============================================================

#' Estimate selection coefficient from frequency trajectory
#'
#' Given a vector of allele frequencies observed over equally-spaced
#' time points, estimates the selection coefficient \eqn{s} using
#' the change in logit-transformed frequency per unit time:
#' \deqn{s = \frac{\Delta \mathrm{logit}(p)}{\Delta t}}
#' where \eqn{\mathrm{logit}(p) = \log(p / (1 - p))}.
#'
#' This assumes a simple deterministic selection model. Each
#' consecutive pair of time points yields one estimate of \eqn{s}.
#'
#' @param freqs Numeric vector. Allele frequencies at successive
#'   time points (values must be in (0, 1)).
#' @param times Numeric vector (optional). Time values for each
#'   observation. If \code{NULL}, integer time steps 1, 2, ... are
#'   used.
#' @param generation_time Numeric scalar. Time per generation
#'   (default: 1). Used to scale the selection coefficient.
#'
#' @return Numeric vector of selection coefficient estimates, one
#'   per consecutive time interval. Length is
#'   \code{length(freqs) - 1}.
#'
#' @references
#' Feder AF et al. (2016). More effective drugs lead to harder
#' selective sweeps in the evolution of drug resistance in HIV-1.
#' \emph{eLife} 5:e10670. \doi{10.7554/eLife.10670}
#'
#' @export
#' @examples
#' freqs <- c(0.05, 0.10, 0.25, 0.50)
#' estimateSelectionCoefficient(freqs)
#' estimateSelectionCoefficient(freqs, times = c(0, 3, 7, 14))
estimateSelectionCoefficient <- function(freqs, times = NULL,
                                         generation_time = 1) {
    if (!is.numeric(freqs) || length(freqs) < 2L)
        stop("'freqs' must be a numeric vector with at least 2 elements.")
    if (any(is.na(freqs)))
        stop("'freqs' must not contain NA values.")
    if (any(freqs <= 0 | freqs >= 1))
        stop("'freqs' must contain values strictly between 0 and 1.")
    if (!is.numeric(generation_time) || length(generation_time) != 1L ||
        generation_time <= 0)
        stop("'generation_time' must be a single positive number.")

    if (is.null(times)) {
        times <- seq_along(freqs)
    } else {
        if (!is.numeric(times) || length(times) != length(freqs))
            stop("'times' must be a numeric vector with the same ",
                 "length as 'freqs'.")
    }

    logit_p <- log(freqs / (1 - freqs))
    delta_logit <- diff(logit_p)
    delta_t <- diff(times) / generation_time

    if (any(delta_t <= 0))
        stop("'times' must contain strictly increasing values.")

    delta_logit / delta_t
}
