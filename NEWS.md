# WithinHostExperiment 0.99.1

Changes made in response to the Bioconductor review. The package has not
been released, so functions are removed and renamed without a deprecation
cycle.

## Interface changes

* `annotateFromGFF()` is removed, together with the package's own GFF3
  and FASTA readers. `annotateGenes(x, genes)` labels each site with the
  gene it falls in, adding a `GENE` column to `mcols(rowRanges(x))`.
  `genes` is a `GRanges`, such as a GFF3 or GTF file read with
  `rtracklayer::import()`, or a `TxDb`, from which
  `GenomicFeatures::genes()` takes the ranges; the overlap is
  `GenomicRanges::findOverlaps()`. `nameCol` selects the label column, and
  an unknown one is an error that lists those available; `collapse` keeps
  the longest overlapping gene or all of them. When the sites and the
  genes share no sequence name the function stops and prints both sets,
  so a RefSeq/GenBank naming mismatch can no longer return a column of
  `NA`.
* `annotateCodonChange()` is an S4 generic, `annotateCodonChange(x, cds,
  ref)`, taking the coding sequences as a `GRanges` of CDS records, a
  `GRangesList` or a `TxDb`, and the reference as a `DNAStringSet`; a file
  path for either is refused with a message naming the Bioconductor reader
  to use. The consequence is computed by
  `VariantAnnotation::predictCoding()`, and the output columns take its
  names: `TXNAME`, `GENEID`, `CONSEQUENCE`, `CDSLOC`, `PROTEINLOC`,
  `REFCODON`, `VARCODON`, `REFAA` and `VARAA` replace `GFF_FEATURE`,
  `REF_CODON`, `ALT_CODON`, `REF_AA`, `ALT_AA` and `AA_CLASS`. A site in
  several overlapping coding sequences is resolved by `collapse`, the
  longest transcript by default. Three classifications change with the
  delegation: losing a stop codon is nonsynonymous rather than nonsense,
  an alternative start codon (`CTG` or `TTG`) at the first base of a
  coding sequence reads as `M`, and a codon holding an ambiguous base is
  left unannotated, with a warning. The reference length is recorded in
  `seqlengths(rowRanges(x))`.
* `dndsWithinHost()` is replaced by `pNpSWithinHost(x, cds, ref)`, which
  takes the same inputs as `annotateCodonChange()` and calls it to
  classify each iSNV. The quantity is a ratio of nonsynonymous to
  synonymous polymorphism within a host, not of substitutions between
  sequences, so the Jukes-Cantor correction is no longer applied.
  Nei-Gojobori sites are counted over the distinct codons of the genome,
  so codons shared by same-frame overlapping coding sequences, such as
  ORF1a and ORF1ab, are counted once, and a warning names the overlapping
  transcripts. On the README's SARS-CoV-2 call set the pooled ratio moves
  from 0.2857 to 0.2867 with the correction removed.
* `consensusSequences(x, ref, minFreq = 0.5)` is new. It writes each
  sample's substitutions at or above `minFreq` into the reference and
  returns one sequence per sample as a `DNAStringSet`; insertions and
  deletions are not applied, so the sequences stay aligned to the
  reference and to one another. `region` returns one coding sequence,
  joined across its parts and read in its own direction. dN/dS between
  hosts can then be estimated with `MSA2dist::dnastring2kaks()`, which is
  suggested, rather than with a model of this package's own.
* `readWithinHost()` keeps its arguments but reads files with
  `VariantAnnotation::readVcf()` instead of parsing VCF text, and splits
  records with several alternative alleles with `expand()`. A per-allele
  field that the header does not declare `Number=A` on such a record
  stops the import, since `readVcf()` would give every allele the first
  allele's value. `genome` takes a `DNAStringSet`, not a file path; its
  lengths, or those of `##contig` header lines, are recorded in
  `seqlengths()`. A file holding more than one sample is refused;
  `as(vcf, "WithinHostExperiment")` converts a `VCF` object with all of
  its samples, and `as(vr, "WithinHostExperiment")` keeps a `VRanges`'
  sequence lengths.
* `calcDiversity()` and `compareDiversity()` lose the `genomeLength`
  argument. The length that `pi` and `watterson` divide by is read from
  `seqlengths(rowRanges(x))`, summed over the sequences in use, and a
  length that is not recorded is an error saying how to record it.
  `calcDiversity()` gains `n_sites_assessed`, the sites at which a sample
  has any call before QC, and `detection_limit`, the strictest `minFreq`
  of the filters applied, so that a diversity of zero after QC can be
  told apart from a sample with no calls. `example_whe` records its
  1000-base segment length; its values are unchanged.
* `validateDiversity()` is renamed `checkDiversityPlausibility()`. Its
  ranges are the package's own heuristics, not literature values, and the
  man page lists every one. Every sample of a `calcDiversity()` table is
  checked, where only the first row was checked before, and the result
  names the sample.
* `ISNVFilter()` loses the `minBaseQual` and `minMapQual` slots and
  arguments: the object holds no per-site quality, so they could be
  recorded but never applied, and quality filtering belongs in the
  variant caller. `replicateConc` is still recorded without being
  applied, and `flagISNV()` warns when it is set.

## Statistical safeguards

* `WithinHostSFS` gains a `polarization` slot: `"none"` for a folded
  spectrum, `"reference"` for an unfolded spectrum of
  alternative-allele frequencies, which treats the reference allele as
  ancestral, or `"outgroup"`. Validity ties it to `folded`, and `show()`
  prints it and labels the stored depth a median, which it is.
* `neutralityFromSFS()` stops on a folded spectrum, where it printed a
  note and returned a value, because Fay and Wu's H needs derived-allele
  frequencies; its result records the polarisation it used.
* `compareSFS()` stops unless both spectra share their bin breaks, their
  folding and their polarisation and neither is bias-corrected; it used
  to compare any two spectra with the same number of bins.
* `tajimaD()` returns `NA`, with a warning, below three segregating
  sites, where it returned 0 at none and a number at one or two, and
  when the sample size rounds below four, where the variance of D is
  zero and it returned `Inf` or `NaN`.
  `calcNeutralityTests()` reports `n_sites`, the number of sites each
  statistic rests on, and warns once with the number of samples below
  that rule.
* `calcSharedVariants()` gains `minDepth`, by default the strictest
  `minDepth` applied with `flagISNV()` (100 when none was), and splits the
  donor-only count into `n_donor_only_adequate_depth`,
  `n_donor_only_low_depth` and `n_donor_only_unknown_depth`; only the
  first is evidence that a variant was lost in transmission.
  `n_donor_only` is still their sum.
* `exactBottleneck(method = "wright_fisher")` says in a message that its
  result is a heuristic score, not a likelihood; `exact_bb` and
  `presence_absence` are the estimators.

## Dependencies

* Biostrings and GenomicFeatures are imported; BiocStyle, rtracklayer,
  txdbmaker and MSA2dist are suggested. rtracklayer is not imported.
* graphics and tools, which the package already used, are declared in
  Imports; covr, which nothing used, is no longer suggested.
* `LongRead` is removed from biocViews.

## Documentation

* The vignette is rewritten as teaching material with BiocStyle: each
  code chunk is introduced by what it does and followed by how to read
  its output, terms are defined where they are introduced and collected
  in a glossary, the `qcLog` and `qcFilters` slots that
  `WithinHostExperiment` adds to `RangedSummarizedExperiment` are
  described one by one, and installation is from Bioconductor only.
* The README marks installation from GitHub as not recommended and
  states the source and terms of use of the four Bendall et al. (2023)
  excerpts that ship in `inst/extdata`.
* The pkgdown reference index gives every topic a section with a
  one-line description, in place of a catch-all.
* `inst/scripts/make_extdata.R`, which ships, records how each file in
  `inst/extdata` was made and with which seed.
* `readme_figure_output/` is removed from the repository root. The values
  and checksums the README is checked against are kept in
  `inst/scripts/readme_values/`, and `inst/scripts/reproduce_all.R`
  rebuilds every derived file from the raw data and checks it against
  the committed copy.

## Tests

* `tests/testthat/test-corrections.R` holds a regression test for each
  item of 'Corrections made before submission' under 0.99.0 that still
  applies, named after the bug it guards against, and the README's
  statement that QC cannot raise pi as a test. The removed codon
  annotator is kept in `tests/testthat/helper-reference-codon.R` as a
  reference implementation, and `annotateCodonChange()` is tested
  against it.

# WithinHostExperiment 0.99.0

## Corrections made before submission

Before submission every exported function was read against its
documentation and run on fixtures built to break it, and each
discrepancy found was fixed, documented as implemented, or both. The
items below record what changed. In 0.99.1 each item that still applies
was pinned by a regression test in `tests/testthat/test-corrections.R`,
named after the bug it guards against (for example `exact_bb matches
sites by chrom, pos and allele (was: pos only)`), so that none can
return unnoticed.

* `compareSFS()` suppressed every warning from `chisq.test()`, so a p-value
  resting on expected counts below 5 looked as trustworthy as any other. It
  now returns the smallest expected count and `chisq_approx_ok`, which is
  `FALSE` in that case, and silences only that one warning; the statistic and
  p-value are unchanged.
* `flagTemporalInconsistency()` classified each site with the last host
  processed and flagged transient sites in every sample. Sites are now
  classified within each host and flagged only in the samples of hosts where
  they are transient.
* `exactBottleneck(method = "exact_bb")` matched donor variants to rows by
  position alone, which is wrong for segmented genomes, and paired donor
  frequencies with recipient read counts by index when some sites did not
  match. Sites are now matched by chromosome, position and allele. An estimate
  that reaches `maxNb` triggers a warning.
* The documentation of `exactBottleneck()` states how `exact_bb` departs from
  Sobel Leonard et al. (2017) and that `wright_fisher` is a heuristic score.
* `flagISNV()`, `calcDiversity()` and `calcNeutralityTests()` stop on
  arguments they do not use.
* `adaptiveFilter("ont")` uses `minDepth = 60`, the depth above which Bull et
  al. (2020) found over 99% sensitivity and precision for single-nucleotide
  variants; it was 30.
* `correctSFSBias()`, `compareSFS()`, the `WithinHostSFS` bins and depth slot,
  `consensusISNV()`, `flagReplicateDiscordance()` and `validateDiversity()` are
  documented as implemented: truncation rescales the whole spectrum, bins are
  right-closed, the stored depth is a median, the consensus object holds only
  allele frequencies, and the diversity ranges are heuristic plausibility
  bounds rather than literature values.
* `ISNVFilter()` documents `minBaseQual`, `minMapQual` and `replicateConc` as
  recorded but not applied, and `flagISNV()` warns when they are set, instead
  of ignoring them silently. `WithinHostExperiment()` requires an `altFreq`
  assay, as documented. (0.99.1 removes `minBaseQual` and `minMapQual`.)
* `trackFrequency()` builds its table in one vectorised step, which is much
  faster on genome-wide data and gives the same rows.
* `plotFrequencySpectrum()`, `plotPairScatter()` and `plotQCDashboard()` clip the
  view with coordinate limits instead of removing data outside the axis range,
  which dropped jittered points with a warning; `plotFrequencyTrajectory()`
  uses the Okabe-Ito palette, as described, and `plotSFS()` recycles it for
  more than eight samples.

## New features

* `WithinHostExperiment` S4 class extending
  `RangedSummarizedExperiment` for within-host pathogen iSNV data,
  with `setValidity` enforcing structural invariants.
* Multi-caller import: `readWithinHost()` reads VCFs in the LoFreq,
  Freebayes, generic and iVar-style field layouts with auto-detection
  (`ivar variants` writes a TSV, so an iVar-style VCF comes from a
  converter). `readWithinHostTable()` handles the iVar TSV directly,
  and generic CSV.
* Replicate-aware QC: `flagISNV()` applies configurable quality
  filters without removing data. `flagReplicateDiscordance()`
  identifies discordant technical replicates with full provenance
  logging. `flagPrimerSites()` masks amplicon primer regions.
* Diversity indices: `shannonISNV()`, `piISNV()` (with optional
  finite-sample correction per Nei 1987), `wattersonISNV()`,
  `simpsonISNV()`, `chao1ISNV()`, and `calcDiversity()` for
  per-sample within-host diversity.
* Neutrality tests: `tajimaD()`, `fusFs()`, and
  `calcNeutralityTests()` for site-frequency-spectrum-based
  neutrality testing with deep-sequencing-aware sample size
  capping.
* **Site frequency spectrum as a first-class object**:
  `WithinHostSFS` S4 class stores the within-host SFS with
  metadata (sample, genome length, depth, detection threshold).
  `buildSFS()` / `buildSFSList()` construct from a WHE;
  `correctSFSBias()` applies ascertainment bias correction for
  the detection threshold (truncation or binomial method);
  `neutralityFromSFS()` computes Tajima's D, Fu's Fs, and
  **Fay and Wu's H** (new) from a single SFS object;
  `compareSFS()` tests SFS homogeneity between samples;
  `plotSFS()` provides publication-quality SFS visualisation.
* **Temporal QC for longitudinal data**:
  `flagTemporalInconsistency()` identifies transient iSNVs that
  appear at only one timepoint (likely artefacts) versus persistent
  ones detected across multiple consecutive timepoints. Supports
  configurable `minTimepoints` and `minConsecutive` thresholds.
* **Cross-contamination detection**:
  `detectCrossContamination()` identifies index-hopping signatures
  where a high-frequency variant in one sample appears at
  suspiciously low frequency in a co-sequenced sample, consistent
  with adapter exchange during pooled sequencing. Returns a
  structured DataFrame of suspect events with source/sink
  frequencies and run metadata.
* GFF3 annotation: `annotateFromGFF()` overlaps variant sites with
  gene models from a GFF3 file, adding `GFF_FEATURE` to
  `mcols(rowRanges)` for use with `dndsWithinHost()`.
  `annotateCodonChange()` extends this with CDS-aware codon
  translation: given a GFF3 and reference FASTA, it determines
  the affected codon, translates both reference and alternative
  alleles using the standard genetic code, and classifies each
  SNV as Synonymous, Nonsynonymous, or Nonsense. (0.99.1 replaces
  both with functions that take Bioconductor objects and delegates the
  consequence call to `VariantAnnotation::predictCoding()`.)
* Selection analysis: `dndsWithinHost()` estimates within-host dN/dS
  with the Nei-Gojobori counting method, normalising by the synonymous
  and nonsynonymous sites of the CDS features in a GFF3 and reference
  FASTA, with Jukes-Cantor correction; `estimateSelectionCoefficient()`
  for time-series selection coefficient inference.
* Multi-caller consensus: `consensusISNV()` for deriving consensus
  variant calls across callers, `calcCallerConcordance()` for
  inter-caller agreement metrics, `benchmarkCallers()` for caller
  performance comparison.
* Longitudinal analysis: `trackFrequency()` for tracking variant
  allele frequencies over time, `detectEmergingVariants()` for
  identifying newly arising variants, `plotFrequencyTrajectory()`
  for publication-quality trajectory plots.
* Transmission analysis: `transmissionPairs()`,
  `exportPairFrequencies()`, `calcSharedVariants()`, and
  `quickBottleneck()` (approximate beta-binomial MLE, Sobel
  Leonard et al. 2017) for bottleneck estimation.
* Exact bottleneck estimation: `exactBottleneck()` offers a
  beta-binomial likelihood over the transmitted virions that models
  finite recipient depth, in the spirit of Sobel Leonard et al. (2017)
  with two stated departures, a presence-absence likelihood, and a
  Wright-Fisher score, which is a heuristic rather than a likelihood.
* Strand bias detection: `sciStrandBias()` and
  `flagStrandBiasSCI()` for strand bias filtering that accounts
  for unequal strand coverage.
* Adaptive filtering: `adaptiveFilter()` for technology-specific
  QC presets (Illumina, ONT, PacBio HiFi).
* Diversity modelling: `fitDiversityModel()` for regression
  modelling of diversity measures, `multitestCorrection()` for
  multiple testing correction, `permutationTest()` for
  non-parametric hypothesis testing.
* Validation: `validateDiversity()` compares diversity estimates with
  wide heuristic plausibility ranges per pathogen, not with literature
  values.
* ViralBottleneck bridge: `asViralBottleneckInput()` exports
  standardised data for the ViralBottleneck package (Zheng
  et al. 2025).
* Interoperability: `as(whe, "VRanges")`,
  `as(vr, "WithinHostExperiment")`, and
  `as(whe, "data.frame")` for seamless integration.
* Batch operations: methods for the S4Vectors `combineRows()` and
  `combineCols()` generics merge experiments together with their QC
  logs and filters.
* Publication-quality visualisation: `plotFrequencySpectrum()`,
  `plotPairScatter()`, `plotQCDashboard()`, and
  `plotFrequencyTrajectory()` with a colour-blind-safe palette
  (Wong 2011).
* SARS-CoV-2 case study in the vignette, with the same workflow run on
  synthetic examples for influenza A, HIV-1 and M. tuberculosis.
