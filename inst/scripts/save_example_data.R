#!/usr/bin/env Rscript
## Generates data/example_whe.rda from the three synthetic VCFs in
## inst/extdata. Run from the package root:
##   Rscript inst/scripts/save_example_data.R

library(S4Vectors)
pkgload::load_all(".", quiet = TRUE)

vcf1 <- "inst/extdata/test_donor.vcf"
vcf2 <- "inst/extdata/test_recipient.vcf"
vcf3 <- "inst/extdata/test_independent.vcf"

example_whe <- readWithinHost(
    c(vcf1, vcf2, vcf3),
    colData = DataFrame(
        sample_id = c("donor", "recipient", "independent"),
        host_id   = c("P1", "P2", "P3"),
        role      = c("donor", "recipient", "independent"),
        pair_id   = c("pair_1", "pair_1", NA_character_),
        ct_value  = c(18.5, 22.1, 15.3)
    ),
    caller = "ivar",
    genome = Biostrings::readDNAStringSet("inst/extdata/test_ref.fa")
)

if (!dir.exists("data")) dir.create("data")
save(example_whe, file = "data/example_whe.rda", compress = "xz")
cat("Saved data/example_whe.rda\n")
cat("  Dimensions:", nrow(example_whe), "x", ncol(example_whe), "\n")
cat("  Pairs:", nrow(transmissionPairs(example_whe)), "\n")
