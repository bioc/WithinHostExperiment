# tests/testthat/helper-reference-codon.R
#
# The codon annotator the package carried before annotateCodonChange()
# delegated to VariantAnnotation::predictCoding(), kept here as an
# independent reference implementation for the cross-implementation test.
# It calls nothing in the package: the codon table is written out here
# and the codon arithmetic is the removed function's own, so a test
# comparing the two compares hand-written counting with
# predictCoding()/Biostrings::translate() rather than one implementation
# with itself. It reads a GFF3 and a FASTA from disk, as the removed
# function did; the package no longer does.

.ref_codon_table <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G")

.ref_translate <- function(codon) {
    codon <- toupper(codon)
    if (nchar(codon) != 3L) return(NA_character_)
    aa <- .ref_codon_table[codon]
    if (is.null(aa) || is.na(aa)) NA_character_ else unname(aa)
}

.ref_revcomp <- function(seq_str) {
    comp <- c(A = "T", T = "A", C = "G", G = "C",
              a = "t", t = "a", c = "g", g = "c",
              N = "N", n = "n")
    chars <- strsplit(seq_str, "")[[1L]]
    paste(rev(comp[chars]), collapse = "")
}

.ref_read_fasta <- function(path) {
    lines <- readLines(path)
    header_idx <- grep("^>", lines)
    if (length(header_idx) == 0L)
        stop("No FASTA header lines found in: ", path)
    seqs <- list()
    for (i in seq_along(header_idx)) {
        name <- sub("\\s+.*", "", sub("^>\\s*", "", lines[header_idx[i]]))
        from <- header_idx[i] + 1L
        to <- if (i < length(header_idx)) header_idx[i + 1L] - 1L
              else length(lines)
        body <- if (from <= to) lines[from:to] else character()
        seqs[[name]] <- gsub("[[:space:]]", "", paste(body, collapse = ""))
    }
    seqs
}

.ref_read_gff_cds <- function(gff, nameCol = "gene") {
    lines <- readLines(gff)
    lines <- lines[!startsWith(lines, "#") & nchar(lines) > 0]
    fields <- strsplit(lines, "\t")
    fields <- fields[vapply(fields, length, integer(1)) >= 9L]
    if (length(fields) == 0L) stop("No valid GFF lines found.")
    gff_df <- data.frame(
        seqid = vapply(fields, `[`, character(1), 1L),
        type = vapply(fields, `[`, character(1), 3L),
        start = as.integer(vapply(fields, `[`, character(1), 4L)),
        end = as.integer(vapply(fields, `[`, character(1), 5L)),
        strand = vapply(fields, `[`, character(1), 7L),
        phase = vapply(fields, `[`, character(1), 8L),
        attrs = vapply(fields, `[`, character(1), 9L),
        stringsAsFactors = FALSE)
    cds <- gff_df[gff_df$type == "CDS", , drop = FALSE]
    if (nrow(cds) == 0L) stop("No CDS features found in GFF file.")
    cds$phase <- as.integer(ifelse(cds$phase == ".", "0", cds$phase))
    parse_attr <- function(attr_str, key) {
        for (p in strsplit(attr_str, ";")[[1L]]) {
            kv <- strsplit(trimws(p), "=")[[1L]]
            if (length(kv) == 2L && kv[1L] == key) return(kv[2L])
        }
        NA_character_
    }
    cds$gene <- vapply(cds$attrs, parse_attr, character(1),
                       key = nameCol, USE.NAMES = FALSE)
    if (all(is.na(cds$gene)))
        cds$gene <- vapply(cds$attrs, parse_attr, character(1),
                           key = "ID", USE.NAMES = FALSE)
    rownames(cds) <- NULL
    cds
}

## Returns a data.frame with one row per variant site and the columns the
## removed annotator wrote: GFF_FEATURE, REF_CODON, ALT_CODON, REF_AA,
## ALT_AA and AA_CLASS (Synonymous / Nonsynonymous / Nonsense).
refAnnotateCodonChange <- function(whe, gff, refFasta, nameCol = "gene") {
    cds <- .ref_read_gff_cds(gff, nameCol)
    ref_seqs <- .ref_read_fasta(refFasta)
    cds_gr <- GenomicRanges::GRanges(cds$seqid,
        IRanges::IRanges(cds$start, cds$end))
    rr <- SummarizedExperiment::rowRanges(whe)
    mc <- S4Vectors::mcols(rr)
    n <- length(rr)
    out <- data.frame(
        GFF_FEATURE = rep(NA_character_, n), REF_CODON = NA_character_,
        ALT_CODON = NA_character_, REF_AA = NA_character_,
        ALT_AA = NA_character_, AA_CLASS = NA_character_,
        stringsAsFactors = FALSE)

    hits <- GenomicRanges::findOverlaps(rr, cds_gr)
    for (h in seq_along(S4Vectors::queryHits(hits))) {
        qi <- S4Vectors::queryHits(hits)[h]
        si <- S4Vectors::subjectHits(hits)[h]
        ref_base <- toupper(as.character(mc$ref[qi]))
        alt_base <- toupper(as.character(mc$alt[qi]))
        if (is.na(ref_base) || is.na(alt_base)) next
        if (nchar(ref_base) != 1L || nchar(alt_base) != 1L) next
        if (grepl("[^ACGT]", ref_base) || grepl("[^ACGT]", alt_base)) next
        var_pos <- GenomicRanges::start(rr)[qi]
        seq_name <- as.character(GenomeInfoDb::seqnames(rr)[qi])
        if (!seq_name %in% names(ref_seqs)) next
        genome_seq <- ref_seqs[[seq_name]]
        if (cds$strand[si] == "+") {
            pos_in_cds <- (var_pos - cds$start[si]) - cds$phase[si]
            codon_start <- cds$start[si] + cds$phase[si] +
                (pos_in_cds %/% 3L) * 3L
        } else {
            pos_in_cds <- (cds$end[si] - var_pos) - cds$phase[si]
            codon_start <- cds$end[si] - cds$phase[si] -
                (pos_in_cds %/% 3L) * 3L - 2L
        }
        base_in_codon <- pos_in_cds %% 3L
        if (pos_in_cds < 0L || codon_start < 1L) next
        if (codon_start + 2L > nchar(genome_seq)) next
        ref_codon <- toupper(substr(genome_seq, codon_start,
                                    codon_start + 2L))
        if (nchar(ref_codon) != 3L) next
        if (cds$strand[si] == "-") {
            ref_codon <- .ref_revcomp(ref_codon)
            comp <- c(A = "T", T = "A", C = "G", G = "C")
            ref_base <- comp[ref_base]
            alt_base <- comp[alt_base]
        }
        alt_codon <- ref_codon
        substr(alt_codon, base_in_codon + 1L, base_in_codon + 1L) <-
            alt_base
        r_aa <- .ref_translate(ref_codon)
        a_aa <- .ref_translate(alt_codon)
        if (is.na(r_aa) || is.na(a_aa)) next
        out$GFF_FEATURE[qi] <- cds$gene[si]
        out$REF_CODON[qi] <- ref_codon
        out$ALT_CODON[qi] <- alt_codon
        out$REF_AA[qi] <- r_aa
        out$ALT_AA[qi] <- a_aa
        out$AA_CLASS[qi] <- if (r_aa == a_aa) "Synonymous"
            else if (r_aa == "*" || a_aa == "*") "Nonsense"
            else "Nonsynonymous"
    }
    out
}
