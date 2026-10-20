# R/annotate-codon.R
# CDS-aware codon change annotation from GFF3 + reference FASTA

#' @include annotate-genes.R
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges start findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
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
.translate_codon <- function(codon) {
    codon <- toupper(codon)
    if (nchar(codon) != 3L) return(NA_character_)
    aa <- .CODON_TABLE[codon]
    if (is.null(aa) || is.na(aa)) NA_character_ else unname(aa)
}

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
.revcomp <- function(seq_str) {
    comp <- c(A = "T", T = "A", C = "G", G = "C",
              a = "t", t = "a", c = "g", g = "c",
              N = "N", n = "n")
    chars <- strsplit(seq_str, "")[[1L]]
    paste(rev(comp[chars]), collapse = "")
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

#' Annotate codon changes from GFF3 CDS + reference FASTA
#'
#' For each SNV in a \code{\link{WithinHostExperiment}}, determines
#' the affected codon and amino acid change using CDS coordinates
#' from a GFF3 file and the reference genome sequence from a FASTA
#' file. Adds \code{GFF_FEATURE}, \code{REF_CODON}, \code{ALT_CODON},
#' \code{REF_AA}, \code{ALT_AA}, and \code{AA_CLASS} (Synonymous /
#' Nonsynonymous / Nonsense) columns to \code{mcols(rowRanges)}.
#'
#' This is a pure-R implementation using the standard genetic code
#' (NCBI translation table 1). No external annotation tools are
#' required. For complex genomes with overlapping reading frames,
#' ribosomal frameshift sites, or non-standard codon usage, use
#' SnpEff or ANNOVAR upstream.
#'
#' @param whe A \code{\link{WithinHostExperiment}} object. Must have
#'   \code{ref} and \code{alt} columns in \code{mcols(rowRanges)}.
#' @param gff Character scalar. Path to a GFF3 file containing CDS
#'   features.
#' @param refFasta Character scalar. Path to the reference genome
#'   FASTA file.
#' @param nameCol Character scalar. GFF3 attribute for the gene name
#'   (default: \code{"gene"}). Falls back to \code{"Name"}, then
#'   \code{"ID"}.
#'
#' @return A \code{\link{WithinHostExperiment}} with additional
#'   columns in \code{mcols(rowRanges)}: \code{GFF_FEATURE},
#'   \code{REF_CODON}, \code{ALT_CODON}, \code{REF_AA},
#'   \code{ALT_AA}, \code{AA_CLASS}.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Parse CDS features from GFF3 (requires column 8 = phase).
#'   \item For each variant overlapping a CDS, compute the codon
#'     position within the CDS (accounting for strand and phase).
#'   \item Extract the reference codon from the FASTA, substitute
#'     the variant base, translate both codons.
#'   \item Classify: Synonymous (same AA), Nonsynonymous (different
#'     AA), or Nonsense (stop codon gained/lost).
#' }
#'
#' Only single-nucleotide substitutions are annotated. Indels and
#' multi-nucleotide variants are marked as \code{NA}.
#'
#' @export
#' @examples
#' ## Requires a reference FASTA and GFF3 with CDS entries
#' ref <- system.file("extdata", "test_ref.fa",
#'     package = "WithinHostExperiment")
#' if (file.exists(ref)) {
#'     gff_file <- tempfile(fileext = ".gff3")
#'     writeLines(c("##gff-version 3",
#'         "test_segment\t.\tCDS\t1\t999\t.\t+\t0\tgene=testGene"),
#'         gff_file)
#'     vcf <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(vcf,
#'         colData = S4Vectors::DataFrame(sample_id = "donor"),
#'         caller = "ivar")
#'     whe <- annotateCodonChange(whe, gff_file, ref)
#'     mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))
#'     head(mc[, c("GFF_FEATURE", "REF_AA", "ALT_AA", "AA_CLASS")])
#' }
annotateCodonChange <- function(whe, gff, refFasta,
                                nameCol = "gene") {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!file.exists(gff))
        stop("GFF file not found: ", gff)
    if (!file.exists(refFasta))
        stop("Reference FASTA not found: ", refFasta)

    ## ---- Parse GFF3 CDS features ----
    cds <- .read_gff_cds(gff, nameCol)

    ## ---- Read reference FASTA ----
    ref_seqs <- .read_fasta_simple(refFasta)

    ## ---- Build CDS GRanges ----
    cds_gr <- GRanges(seqnames = cds$seqid,
                      ranges = IRanges(cds$start, cds$end))

    ## ---- Annotate each variant ----
    rr <- rowRanges(whe)
    mc <- mcols(rr)

    if (!"ref" %in% colnames(mc) || !"alt" %in% colnames(mc))
        stop("mcols(rowRanges) must contain 'ref' and 'alt' columns.")

    n <- length(rr)
    gene_vec   <- rep(NA_character_, n)
    ref_codon  <- rep(NA_character_, n)
    alt_codon  <- rep(NA_character_, n)
    ref_aa     <- rep(NA_character_, n)
    alt_aa     <- rep(NA_character_, n)
    aa_class   <- rep(NA_character_, n)

    hits <- findOverlaps(rr, cds_gr)
    q_idx <- queryHits(hits)
    s_idx <- subjectHits(hits)

    for (h in seq_along(q_idx)) {
        qi <- q_idx[h]
        si <- s_idx[h]

        ref_base <- toupper(as.character(mc$ref[qi]))
        alt_base <- toupper(as.character(mc$alt[qi]))
        ## Skip indels and multi-nt
        if (is.na(ref_base) || is.na(alt_base)) next
        if (nchar(ref_base) != 1L || nchar(alt_base) != 1L) next
        if (grepl("[^ACGT]", ref_base) || grepl("[^ACGT]", alt_base)) next

        var_pos <- GenomicRanges::start(rr)[qi]
        cds_start  <- cds$start[si]
        cds_end    <- cds$end[si]
        cds_strand <- cds$strand[si]
        cds_phase  <- cds$phase[si]

        ## Get reference sequence
        seq_name <- as.character(GenomeInfoDb::seqnames(rr)[qi])
        if (!seq_name %in% names(ref_seqs)) next
        genome_seq <- ref_seqs[[seq_name]]

        ## GFF3 phase is the number of bases to skip from the start of the
        ## CDS (the 'end' coordinate on the minus strand) to reach the first
        ## complete codon.
        if (cds_strand == "+") {
            pos_in_cds <- (var_pos - cds_start) - cds_phase
            codon_idx  <- pos_in_cds %/% 3L
            base_in_codon <- pos_in_cds %% 3L
            codon_start <- cds_start + cds_phase + codon_idx * 3L
        } else {
            pos_in_cds <- (cds_end - var_pos) - cds_phase
            codon_idx  <- pos_in_cds %/% 3L
            base_in_codon <- pos_in_cds %% 3L
            codon_start <- cds_end - cds_phase - codon_idx * 3L - 2L
        }

        ## Variants in the phase offset are not in a complete codon
        if (pos_in_cds < 0L) next
        if (codon_start < 1L) next
        codon_end <- codon_start + 2L
        if (codon_end > nchar(genome_seq)) next

        ref_codon_str <- toupper(substr(genome_seq, codon_start,
                                        codon_end))
        if (nchar(ref_codon_str) != 3L) next

        if (cds_strand == "-") {
            ref_codon_str <- .revcomp(ref_codon_str)
            ## Also complement the variant bases for minus strand
            comp <- c(A = "T", T = "A", C = "G", G = "C")
            ref_base <- comp[ref_base]
            alt_base <- comp[alt_base]
        }

        ## Substitute
        alt_codon_str <- ref_codon_str
        substr(alt_codon_str, base_in_codon + 1L,
               base_in_codon + 1L) <- alt_base

        r_aa <- .translate_codon(ref_codon_str)
        a_aa <- .translate_codon(alt_codon_str)

        if (is.na(r_aa) || is.na(a_aa)) next

        gene_vec[qi]  <- cds$gene[si]
        ref_codon[qi] <- ref_codon_str
        alt_codon[qi] <- alt_codon_str
        ref_aa[qi]    <- r_aa
        alt_aa[qi]    <- a_aa
        aa_class[qi]  <- if (r_aa == a_aa) {
            "Synonymous"
        } else if (r_aa == "*" || a_aa == "*") {
            "Nonsense"
        } else {
            "Nonsynonymous"
        }
    }

    mcols(rr)$GFF_FEATURE <- gene_vec
    mcols(rr)$REF_CODON   <- ref_codon
    mcols(rr)$ALT_CODON   <- alt_codon
    mcols(rr)$REF_AA      <- ref_aa
    mcols(rr)$ALT_AA      <- alt_aa
    mcols(rr)$AA_CLASS    <- aa_class

    whe_new <- whe
    SummarizedExperiment::rowRanges(whe_new) <- rr
    whe_new
}
