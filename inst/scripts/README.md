# Provenance of `inst/extdata`

Every file installed in `inst/extdata` is listed below with the script
that produced it and the source it came from. All paths are relative to
the package root, and every script in this directory is written to be
run from the package root.

`inst/extdata` holds nineteen files: four derived from published data,
twelve written by a generator script in this directory, and three that
have no generator in the repository. `make_extdata.R` runs the
generators in order, with the seed each uses, and checks the three
files without one against their recorded checksums; `reproduce_all.R`
goes further and rebuilds every derived file of the repository, from
these fixtures to the README figures, checking each against the
committed copy.

## Real data: the HH6 transmission pair

| File | Produced by | Source |
|------|-------------|--------|
| `MHM1143_rep1.tsv` | `prepare_lauring_data.R` | Bendall *et al.* 2023 |
| `MHM1143_rep2.tsv` | `prepare_lauring_data.R` | Bendall *et al.* 2023 |
| `MHM1184_rep1.tsv` | `prepare_lauring_data.R` | Bendall *et al.* 2023 |
| `MHM1184_rep2.tsv` | `prepare_lauring_data.R` | Bendall *et al.* 2023 |

These are the only files in `inst/extdata` that come from a measurement.

Bendall EE *et al.* (2023) Rapid transmission and tight bottlenecks
constrain the evolution of highly transmissible SARS-CoV-2 variants.
*Nature Communications* 14:272. <https://doi.org/10.1038/s41467-023-36001-5>

The upstream repository is
<https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck>.
`prepare_lauring_data.R` reads two files from it, kept under
`inst/scripts/real_data/`: `all_variants_filtered.tsv` (159 iSNV calls
across 83 samples, both replicates in one row per call) and
`Transmission_pairs.csv` (household transmission pair metadata).

The script selects household HH6, whose pair `HH6_A` is MHM1143
(transmission individual A, the donor) and MHM1184 (individual B, the
recipient). It splits each merged row into the two replicates, writing
an iVar-style table per sample and replicate, and then copies the four
HH6 files into `inst/extdata`. The pair carries four iSNVs in total,
two in each host, none of them at a position shared between the two
hosts. The script also writes the HH17 samples (MHM1435, MHM1436,
MHM1437; MHM1433 has no calls) to `inst/scripts/real_data/`, and those
are not copied into `inst/extdata`.

`inst/scripts/real_data/` is excluded from the built tarball by
`.Rbuildignore`, so the two input files are available only in a clone
of the repository. The script itself needs no network access.

### Terms of use

The upstream repository states no licence. The four files reproduced
here are a small excerpt of the study's filtered call table -- four of
its 159 calls -- kept so that the documentation and the tests can run
against real calls, and attributed to the paper and to the repository
above. Anyone wanting these data for their own work should take them
from the paper and the upstream repository.

## Synthetic fixtures with a generator

`create_test_data.R` writes the following twelve files. It simulates
iSNVs on a 1000 bp reference drawn at random, then emits the same
variants in several caller formats so that the import paths can be
exercised against one another. The seed is fixed at the top of the
script, and re-running it reproduces all twelve files byte for byte
(checked with R 4.5.3).

| File | Description |
|------|-------------|
| `test_ref.fa` | 1000 bp random reference sequence, `test_segment` |
| `test_donor.vcf` | Donor sample, 15 iSNVs, iVar-style VCF |
| `test_recipient.vcf` | Recipient sample, 10 iSNVs, 5 shared with the donor |
| `test_independent.vcf` | Unrelated sample, 8 iSNVs |
| `test_donor_rep1.vcf` | Technical replicate of the donor, 16 sites: two deliberately discordant, one replicate-only site at position 500 |
| `test_donor_rep2.vcf` | Second technical replicate, 15 sites, small perturbations only, position 500 absent |
| `test_donor_ivar.tsv` | The donor variants as an iVar TSV table |
| `test_donor_lofreq.vcf` | The donor variants in LoFreq VCF format, with `DP4` strand counts |
| `test_donor_freebayes.vcf` | The donor variants in Freebayes VCF format |
| `test_generic.csv` | The donor variants as a minimal generic CSV |
| `test_coldata.csv` | Sample metadata for the five synthetic samples |
| `test_primers.bed` | Two primer intervals, for primer-site masking |

Site 13 of the donor (position 812) is given a depth of 150 so that the
depth filter has something to remove, and the two discordant replicate
sites and the replicate-only site exist so that replicate-aware QC has
something to flag.

## Fixtures with no generator

| File | Description |
|------|-------------|
| `test_influenza_donor.tsv` | 12 rows on an influenza HA segment |
| `test_hiv_sample.tsv` | 20 rows on an HIV-1 *pol* region |
| `test_tb_sample.tsv` | 3 rows on an *M. tuberculosis* gene |

No script in this repository produces these three files. They have been
present unchanged since the first commit, and how they were made was
not recorded. They are synthetic: none of them is derived from any
measurement, and none should be read as one. They exist so that the
vignette and the tests can show `checkDiversityPlausibility()` and the
neutrality tests running on something other than SARS-CoV-2, and the
vignette says of them that they are illustrative fixtures rather than
measurements.

What can still be said about how they were made comes from the numbers
themselves. In all three files the read counts are internally
consistent: `REF_DP + ALT_DP` equals `TOTAL_DP` in every row, and
`ALT_DP` equals `TOTAL_DP * ALT_FREQ` rounded to the nearest integer.
In `test_influenza_donor.tsv` and `test_hiv_sample.tsv` the
frequencies carry fourteen to sixteen decimal places and do not equal
`ALT_DP / TOTAL_DP`, which is what happens when frequencies are drawn
from a continuous distribution and the counts are then derived from
them -- so those two files were written by a program, not typed. In
`test_tb_sample.tsv` the three frequencies are 0.04, 0.08 and 0.03
against depths of 800, 1200 and 600, and each frequency is exactly
`ALT_DP / TOTAL_DP`, which is equally consistent with the values having
been chosen by hand.

Regenerating these three files would change the numbers reported in the
vignette and the expectations in `tests/testthat/test-validation.R`, so
they are kept as they are rather than replaced by a new generator.

## Scripts in this directory

| Script | In the tarball | Purpose |
|--------|----------------|---------|
| `create_test_data.R` | yes | Writes the twelve synthetic fixtures listed above into `inst/extdata` |
| `prepare_lauring_data.R` | yes | Builds the per-sample replicate tables from the Bendall *et al.* data and copies the HH6 pair into `inst/extdata` |
| `save_example_data.R` | yes | Builds `data/example_whe.rda` from the three synthetic VCFs in `inst/extdata` |
| `generate_readme_figures.R` | no | Regenerates the six README figures from `inst/scripts/real_data` |
| `fill_readme.R` | no | Fills the README captions from the values the figure script records |
| `README_template.md` | no | The README source that `fill_readme.R` fills |
| `prepare_farjo_metadata.R` | no | Builds the Farjo *et al.* sample metadata; this one does need network access |
| `verify_claims.R` | no | Ad hoc check that the documented capabilities run |
| `real_data/` | no | Inputs for the figure script and for `prepare_lauring_data.R` |

The entries marked "no" are excluded by `.Rbuildignore` and exist only
in a clone of the repository. `README.md` in the package root has the
inventory of `inst/scripts/real_data/`.
