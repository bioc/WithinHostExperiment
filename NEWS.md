# WithinHostExperiment 0.99.0

## Corrections made before submission

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
  `consensusISNV()`, `flagReplicateDiscordance()` and
  `checkDiversityPlausibility()` are
  documented as implemented: truncation rescales the whole spectrum, bins are
  right-closed, the stored depth is a median, the consensus object holds only
  allele frequencies, and the diversity ranges are heuristic plausibility
  bounds rather than literature values.
* `ISNVFilter()` documents `minBaseQual`, `minMapQual` and `replicateConc` as
  recorded but not applied, and `flagISNV()` warns when they are set, instead
  of ignoring them silently. `WithinHostExperiment()` requires an `altFreq`
  assay, as documented.
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
  converter). The files are read with `VariantAnnotation::readVcf()`,
  and records with several alternative alleles are split with
  `expand()`; a per-allele field that the header does not declare
  `Number=A` on such a record stops the import, since `readVcf()` would
  give every allele the first allele's value. Sequence lengths from
  `##contig` lines, or from a `DNAStringSet` passed as `genome`, are
  recorded in `seqlengths()`. A VCF already in R converts with
  `as(vcf, "WithinHostExperiment")`, taking every sample it holds, and
  `as(vr, "WithinHostExperiment")` now keeps a `VRanges`' sequence
  lengths. `readWithinHostTable()` handles the iVar TSV directly, and
  generic CSV.
* Replicate-aware QC: `flagISNV()` applies configurable quality
  filters without removing data. `flagReplicateDiscordance()`
  identifies discordant technical replicates with full provenance
  logging. `flagPrimerSites()` masks amplicon primer regions.
* Diversity indices: `shannonISNV()`, `piISNV()` (with optional
  finite-sample correction per Nei 1987), `wattersonISNV()`,
  `simpsonISNV()`, `chao1ISNV()`, and `calcDiversity()` for
  per-sample within-host diversity. `calcDiversity()` reads the genome
  length that `pi` and `watterson` divide by from
  `seqlengths(rowRanges(x))`, where `readWithinHost()` and
  `annotateCodonChange()` record it, and stops with instructions when
  it is not recorded; the `genomeLength` argument is gone from it and
  from `compareDiversity()`. Its output gains `n_sites_assessed`, the
  sites at which a sample has any call, and `detection_limit`, the
  strictest `minFreq` of the filters applied, so a diversity of zero
  can be told apart from a sample with no calls. `example_whe` records
  its 1000-base segment length.
* Neutrality tests: `tajimaD()`, `fusFs()`, and
  `calcNeutralityTests()` for site-frequency-spectrum-based
  neutrality testing with deep-sequencing-aware sample size
  capping. `tajimaD()` returns `NA`, with a warning, below three
  segregating sites, where D is a function of one or two frequencies
  rather than of a spectrum; `calcNeutralityTests()` reports
  `n_sites`, the number of sites each statistic rests on.
* **Site frequency spectrum as a first-class object**:
  `WithinHostSFS` S4 class stores the within-host SFS with
  metadata (sample, genome length, depth, detection threshold).
  `buildSFS()` / `buildSFSList()` construct from a WHE;
  `correctSFSBias()` applies ascertainment bias correction for
  the detection threshold (truncation or binomial method);
  `neutralityFromSFS()` computes Tajima's D, Fu's Fs, and
  **Fay and Wu's H** (new) from a single SFS object;
  `compareSFS()` tests SFS homogeneity between samples. The object
  records its `polarization`: `"none"` when folded, `"reference"` for
  an unfolded spectrum of alternative-allele frequencies, which treats
  the reference allele as ancestral, or `"outgroup"`.
  `neutralityFromSFS()` stops on a folded spectrum, since Fay and Wu's
  H needs derived-allele frequencies, and reports the polarisation it
  used. `compareSFS()` stops unless both spectra share their breaks,
  folding and polarisation and neither is bias-corrected, the
  conditions its chi-squared test assumes;
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
* Gene annotation: `annotateGenes()` labels variant sites with the
  genes they overlap, adding `GENE` to `mcols(rowRanges)`. Genes are
  supplied as a `GRanges` (a GFF3 or GTF read with
  `rtracklayer::import()`) or a `TxDb`, and the function stops, naming
  both sets, when the sites and the genes share no sequence name.
  `annotateCodonChange()` reports the codon and amino acid consequence
  of each SNV. It takes the coding sequences as a `GRanges` of CDS
  records, a `GRangesList` or a `TxDb`, and the reference sequence as a
  `DNAStringSet`, and delegates the consequence call to
  `VariantAnnotation::predictCoding()`, so the translation and the
  consequence come from Bioconductor's own annotation machinery rather
  than from a codon table carried in this package. It adds `TXNAME`,
  `GENEID`, `CONSEQUENCE`, `CDSLOC`, `PROTEINLOC`, `REFCODON`,
  `VARCODON`, `REFAA` and `VARAA` to `mcols(rowRanges)`; a site lying
  in several overlapping coding sequences is resolved by the `collapse`
  argument. Two classifications therefore follow `predictCoding()`:
  losing a stop codon counts as nonsynonymous rather than nonsense, and
  an alternative start codon (`CTG` or `TTG`) at the first base of a
  coding sequence reads as `M`. A codon holding an ambiguous base is
  left unannotated, with a warning, rather than translated.
* Selection analysis: `pNpSWithinHost()` computes the within-host
  ratio of nonsynonymous to synonymous polymorphism per site, pN/pS,
  from the same coding sequences (`GRanges`, `GRangesList` or `TxDb`)
  and `DNAStringSet` reference as `annotateCodonChange()`, which it
  calls to classify each iSNV. Sites follow Nei and Gojobori (1986) and
  are counted over the distinct codons of the genome, so codons that
  same-frame overlapping coding sequences share (ORF1a and ORF1ab) are
  counted once, and a warning names the overlapping transcripts. It
  replaces `dndsWithinHost()`, which took GFF3 and FASTA paths, parsed
  them itself and applied a Jukes-Cantor correction; a multiple-hit
  correction belongs to divergence between sequences, not to
  polymorphism within a host, and is no longer applied.
  `consensusSequences()` writes each sample's substitutions at or
  above a frequency into the reference and returns the sequences as a
  `DNAStringSet`, so dN/dS between hosts can be estimated with the
  models of `MSA2dist::dnastring2kaks()` rather than reimplemented
  here. `estimateSelectionCoefficient()`
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
  `calcSharedVariants()` splits the donor-only count by the
  recipient's read depth at the site: absent at adequate depth, absent
  at low depth, or depth unknown, which is the case at every site a
  VCF does not list.
* Exact bottleneck estimation: `exactBottleneck()` offers a
  beta-binomial likelihood over the transmitted virions that models
  finite recipient depth, in the spirit of Sobel Leonard et al. (2017)
  with two stated departures, a presence-absence likelihood, and a
  Wright-Fisher score. The last is a heuristic rather than a
  likelihood, and the function says so with a message when it is run;
  `exact_bb` and `presence_absence` are the estimators.
* Strand bias detection: `sciStrandBias()` and
  `flagStrandBiasSCI()` for strand bias filtering that accounts
  for unequal strand coverage.
* Adaptive filtering: `adaptiveFilter()` for technology-specific
  QC presets (Illumina, ONT, PacBio HiFi).
* Diversity modelling: `fitDiversityModel()` for regression
  modelling of diversity measures, `multitestCorrection()` for
  multiple testing correction, `permutationTest()` for
  non-parametric hypothesis testing.
* Plausibility check: `checkDiversityPlausibility()` (formerly
  `validateDiversity()`) compares diversity estimates with wide,
  heuristic ranges per pathogen, listed on its man page, as a check on
  units and orders of magnitude; the ranges are the package's own, not
  literature values, and a value inside one is no evidence that an
  estimate is correct. Every sample of a `calcDiversity()` table is
  checked, where only the first row used to be.
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
