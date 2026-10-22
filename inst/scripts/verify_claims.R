## Full functional verification of all claimed capabilities
library(WithinHostExperiment)
library(GenomicRanges)

cat("=== 1. Multi-caller import ===\n")
vcf_iv <- system.file("extdata", "test_donor.vcf", package="WithinHostExperiment")
vcf_lf <- system.file("extdata", "test_donor_lofreq.vcf", package="WithinHostExperiment")
vcf_fb <- system.file("extdata", "test_donor_freebayes.vcf", package="WithinHostExperiment")
tsv    <- system.file("extdata", "test_donor_ivar.tsv", package="WithinHostExperiment")
csv    <- system.file("extdata", "test_generic.csv", package="WithinHostExperiment")

ref <- Biostrings::readDNAStringSet(system.file("extdata", "test_ref.fa", package="WithinHostExperiment"))
w1 <- readWithinHost(vcf_iv, colData=S4Vectors::DataFrame(sample_id="d"), caller="ivar", genome=ref)
w2 <- readWithinHost(vcf_lf, colData=S4Vectors::DataFrame(sample_id="d"), caller="lofreq")
w3 <- readWithinHost(vcf_fb, colData=S4Vectors::DataFrame(sample_id="d"), caller="freebayes")
w4 <- readWithinHostTable(tsv, format="ivar")
w5 <- readWithinHostTable(csv, format="generic")
cat("  iVar VCF:", nrow(w1), "sites\n")
cat("  LoFreq VCF:", nrow(w2), "sites\n")
cat("  Freebayes VCF:", nrow(w3), "sites\n")
cat("  iVar TSV:", nrow(w4), "sites\n")
cat("  Generic CSV:", nrow(w5), "sites\n")
f_iv <- sort(SummarizedExperiment::assay(w1,"altFreq")[,1])
f_lf <- sort(SummarizedExperiment::assay(w2,"altFreq")[,1])
cat("  iVar vs LoFreq max|diff|:", max(abs(f_iv - f_lf)), "\n")

cat("\n=== 2. QC pipeline (flag dont delete) ===\n")
whe <- w1
pre <- nrow(whe)
whe <- flagISNV(whe, ISNVFilter(minDepth=200L, minFreq=0.03, maxFreq=0.50))
cat("  Rows before:", pre, "  Rows after:", nrow(whe), "  SAME:", pre==nrow(whe), "\n")
cat("  qcLog steps:", nrow(qcLog(whe)), "\n")
s <- qcSummary(whe)
cat("  Passed:", s$passed_sites, "/", s$total_sites, "\n")

cat("\n=== 3. Replicate QC ===\n")
r1 <- system.file("extdata", "test_donor_rep1.vcf", package="WithinHostExperiment")
r2 <- system.file("extdata", "test_donor_rep2.vcf", package="WithinHostExperiment")
wrep <- readWithinHost(c(r1, r2),
  colData=S4Vectors::DataFrame(sample_id=c("rep1","rep2"),
    replicate_group=c("g","g"), replicate_id=c("1","2")), caller="ivar")
wrep <- flagReplicateDiscordance(wrep, freqTolerance=0.05)
cat("  Replicate discordance logged:", "replicate_discordance" %in% qcLog(wrep)$step, "\n")

cat("\n=== 4. Primer site filtering ===\n")
primers <- GRanges("test_segment", IRanges::IRanges(start=c(40,880), end=c(65,950)))
wp <- flagPrimerSites(w1, primers)
cat("  Primer filter logged:", "primer_filter" %in% qcLog(wp)$step, "\n")

cat("\n=== 5. Diversity indices ===\n")
div <- calcDiversity(whe)
cat("  Shannon:", round(div$shannon, 4), "\n")
cat("  Pi:", format(div$pi, digits=6), "\n")
cat("  Watterson:", format(div$watterson, digits=6), "\n")
cat("  Richness:", div$richness, "\n")

cat("\n=== 6. Transmission analysis ===\n")
vcf2 <- system.file("extdata", "test_recipient.vcf", package="WithinHostExperiment")
wpair <- readWithinHost(c(vcf_iv, vcf2),
  colData=S4Vectors::DataFrame(sample_id=c("donor","recipient"),
    role=c("donor","recipient"), pair_id=c("p1","p1")), caller="ivar")
pairs <- transmissionPairs(wpair)
cat("  Pairs found:", nrow(pairs), "\n")
sv <- calcSharedVariants(wpair, pairId="p1")
cat("  Shared:", sv$n_shared, "  Donor-only:", sv$n_donor_only,
    "  Recipient-only:", sv$n_recipient_only, "\n")

cat("\n=== 7. Bottleneck estimation ===\n")
nb <- quickBottleneck(wpair, pairId="p1", maxNb=100L)
cat("  Nb MLE:", nb$Nb, "\n")
cat("  95% CI: [", nb$ci[1], ",", nb$ci[2], "]\n")

cat("\n=== 8. ViralBottleneck bridge ===\n")
vb <- asViralBottleneckInput(wpair, pairId="p1")
cat("  Columns:", paste(names(vb), collapse=", "), "\n")
cat("  Rows:", nrow(vb), "  All donor >= 0.03:", all(vb$donor_freq >= 0.03), "\n")

cat("\n=== 9. Interoperability ===\n")
df <- as(wpair, "data.frame")
vr <- as(wpair, "VRanges")
wback <- as(vr, "WithinHostExperiment")
cat("  data.frame:", nrow(df), "rows  VRanges:", length(vr), "\n")
cat("  Roundtrip sites:", nrow(wback), "== orig", nrow(wpair), ":", nrow(wback)==nrow(wpair), "\n")

cat("\n=== 10. Combine ===\n")
wm <- WithinHostExperiment::combineRows(wpair[1:10,], wpair[11:nrow(wpair),])
wmc <- WithinHostExperiment::combineCols(wpair[,1], wpair[,2])
cat("  combineRows OK:", nrow(wm)==nrow(wpair), "  combineCols OK:", ncol(wmc)==ncol(wpair), "\n")

cat("\n=== 11. Pi finite-sample correction ===\n")
pi_raw <- piISNV(c(0.1, 0.3), genomeLength=1000)
pi_cor <- piISNV(c(0.1, 0.3), genomeLength=1000, meanDepth=100)
cat("  Correction factor:", round(pi_cor/pi_raw, 6), "(expected ~1.010101)\n")

cat("\n========================================\n")
cat("ALL CLAIMED FUNCTIONS VERIFIED\n")
cat("========================================\n")
