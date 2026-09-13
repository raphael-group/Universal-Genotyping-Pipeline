# Reference

## Table of Contents
- [Environments](#environments)
- [Configuration](#configuration)
  - [Input Data](#input-data)
    - [Run settings](#run-settings)
    - [File paths](#file-paths)
  - [Parameters](#parameters)
- [Outputs](#outputs)
  - [Genomic unit levels](#genomic-unit-levels)
  - [Final bins](#final-bins)
  - [Unit level](#unit-level)
  - [Intermediates](#intermediates)
  - [TSV columns](#tsv-columns)
  - [QC](#qc-qc_dir)

---

## Environments

Pipeline dependencies live under `workflow/envs/`:

| Environment | Purpose |
|-------------|---------|
| `base.yaml` | Python scientific stack (used by every mode). |
| `bcftools.yaml` | bcftools/tabix: bulk genotyping + het-SNP pileup. |
| `eagle.yaml` | Eagle2 phasing. |
| `shapeit.yaml` | SHAPEIT5 phasing. |
| `longphase.yaml` | LongPhase (long-read) phasing. |
| `cellsnp.yaml` | cellsnp-lite: single-cell genotyping + pileup. |
| `mosdepth.yaml` | mosdepth read-depth counting (bulk). |
| `bedtools.yaml` | bedtools/samtools: per-window read-start counting (bulk). |
| `ucsc.yaml` | UCSC tools (`bigWigToBedGraph`, `liftOver`) for the Repli-seq track. |

---

## Configuration
Defaults in `config/config.yaml`, template in [templates](../resources/templates/). Override with `--config key=value`. Refer to **[spec](sample_sheet.md)** and [templates](../resources/templates/) for sample sheet.

### Input Data

#### Run settings

| Key | Required | Description |
|-----|----------|-------------|
| `workflow_mode` | Yes | `bulk_genotyping` \| `single_cell_genotyping` \| `copytyping_preprocess`. |
| `assay_types` | Yes | Assay types to run, e.g. `["bulkWGS"]`, `["scRNA","scATAC"]`. An assay not valid for `workflow_mode` is dropped; none valid is an error. |
| `sample_id` | Yes | Which `sample_id` of the sample file to process. |
| `dataset_ids` | No | Restrict the run to these `dataset_id`s of that `sample_id`; `[]` (default) runs all of them. Selection precedes the `reference_version` filter: a named `dataset_id` with no record on that build is an error. |
| `chromosomes` | Yes | Chromosomes to run; default `[1..22]`. |
| `species` | Yes | `human` (default) or `mouse`. |
| `reference_version` | Yes | Reference version to select samples. See [Reference version](sample_sheet.md#reference-version). |
| `remote_mode` | Optional | Remote input handling: `storage` (default; download whole file via Snakemake storage) or `stream` (read URLs directly, fetching only `chromosomes`). |
| `genotype_dataset_ids` | Optional | `dataset_id`s piled up to call germline SNPs. Empty -> auto (normal before tumor, short-read before long-read). >1 are pooled in one `mpileup` and must share an `@RG SM` tag. |
| `panel_allele_only` | Optional; `bulk_genotyping`; default `true` | Constrain `bcftools call` to the panel's `REF,ALT` (`--constrain alleles`). `false` uses the panel for positions only, taking `REF` from `reference` and `ALT` from the reads. Forced `true` when a genotyped dataset is a `tumor`. |
| `phase_dataset_ids` | Optional | Datasets used for long-read phasing inputs. |
| `phaser` | Genotyping | `eagle` \| `shapeit` \| `longphase`. |
| `het_snp_vcf_phased` | Optional | Default `true`: the input `het_snp_vcf` is phased or not. |

#### File paths

| Key | Required | Description |
|-----|----------|-------------|
| `sample_file` | Yes | Path to `samples.json`. |
| `reference` | Yes | Genome FASTA. |
| `genome_size` | Yes | Two-column `chrom\tsize` genome size file. |
| `gtf_file` | Yes | Gene annotation GTF (gzipped). |
| `region_bed` | Yes | BED file listing whitelist chromosome arms. |
| `extremity_tsv` | Optional | Headered TSV of upstream SV breakpoints (`#CHR`, `POS0`). Each breakpoint cuts the arm holding it, so no window and no bb spans an SV junction. Setting it ignores any pre-built `window_bed`, since the grid is re-tiled from the cut segments. Template: `resources/templates/extremity.tsv`. |
| `window_bed` | Optional | Pre-built window BED with read depth covariates. Ignored when `extremity_tsv` is set. |
| `mappability_bed` | Optional | BED mappability track (4th column = score). |
| `target_bed` | Optional; `bulkWES` | Hybrid-capture target intervals (BED3+). Marks each window on- or off-target, which `rd_correct` then fits and normalizes separately; see [Whole-exome sequencing (WES)](bulk_genotyping.md#whole-exome-sequencing-wes). Fetch one with `resources/scripts/fetch_capture_targets.sh`. |
| `blacklist_bed` | Optional | ENCODE-style blacklist; pre-built at `resources/data/hg38-blacklist.v2.bed.gz`. |
| `gene_blacklist_file` | Optional | Genes to exclude from AnnData (single-cell), one gene symbol or id per line. |
| `snp_panel` | Genotyping | Population SNP VCF, bgzipped and indexed (`.vcf.gz` + `.tbi`/`.csi`). Bulk passes its positions to `bcftools mpileup -T` and, under `panel_allele_only`, its `REF,ALT` to `bcftools call --constrain alleles`; single-cell passes it to `cellsnp-lite -R`. |
| `phasing_panel` | eagle/shapeit | Per-chromosome BCF reference panel directory. |
| `gmap_path` | eagle/shapeit | Genetic map; `{chrname}` placeholder for per-chromosome maps (SHAPEIT5), literal path for a single map (Eagle2). |
| `het_snp_vcf` | Optional; required for `copytyping_preprocess` | Pre-computed gHET VCF. |
| `bb_file` | copytyping_preprocess | Pre-computed bb annotations TSV. |

> [!IMPORTANT]
> - `genome_size` defines the contig naming convention in `reference` and input alignment files.
> - `snp_panel`, `phasing_panel` and `het_snp_vcf` must follow the same naming
> convention as `genome_size`.
> - Final outputs always use chr-prefix contig naming regardless of input convention.

### Parameters
#### `params_build_windows`
Used by the window-BED build (`build_windows.smk`), which runs in **every** mode.
`build_segment_bed` first cuts the `region_bed` arms at every `extremity_tsv` breakpoint
and subtracts the blacklist into `aux/segment.bed` (region_id + seg_id). One window BED is
then tiled off it, per segment row, so no window and no bin spans a segment bound:
`aux/windows.bed.gz`, shared by every assay of the run.

The windows are the fixed bins of `build_adaptive_bins` in both modes. Bulk additionally
counts them with mosdepth and bias-corrects them in `rd_correct`, which is the only
consumer of the GC/MAP/REPLI columns; a single-cell window BED carries just
`#CHR START END region_id seg_id`, and the Repli-seq fetch is skipped (as it is for any run with `rt_correct: false`).

| Field | Description |
|---|---|
| `window_size` | Window size (bp); fixed tiling of the segment BED, shared by every assay of the run. |

#### `params_bcftools`
Used by `genotype_snps_bulk` and `pileup_snps_bulk_bcftools_chrom` (bulk het-SNP read
counting). Both call one chromosome per job and restrict it with `--regions`.

| Field | Description |
|---|---|
| `min_mapq` | Skip alignments below this mapping quality (genotype + pileup). |
| `min_baseq` | Skip bases below this base quality (genotype + pileup). |
| `min_dp` | Minimum depth to keep a site (genotype only). |
| `max_depth` | Per-file depth cap in `mpileup` (genotype + pileup). |
| `min_qual` | Minimum variant QUAL (genotype only). |
| `extra_params` | Extra `mpileup` flags applied to bulk genotyping + pileup; use a read-type preset, e.g. `-X ont-sup` (ONT), `-X pacbio-ccs`, or `""`/`-X illumina` for short-read. |

#### `params_cellsnp_lite`
Used by `genotype_snps_pseudobulk_mode1b`, `pileup_snps_*` (single-cell).

| Field | Description |
|---|---|
| `UMItag` | UMI tag: `Auto` \| `None` \| a BAM tag. |
| `cellTAG` | Cell-barcode tag: `CB` \| `None` \| a BAM tag. |
| `minMAF_genotype` | Minimum minor-allele frequency when genotyping. |
| `minCOUNT_genotype` | Minimum aggregate count when genotyping. |
| `minMAF_pileup` | Minimum minor-allele frequency when piling up. |
| `minCOUNT_pileup` | Minimum aggregate count when piling up. |

#### `params_genotype_snps`
Genotyping from allele counts, in both modes. One rule: threshold the counts.

| Field | Description |
|---|---|
| `min_dp` | Depth floor; below it a site is not called. |
| `min_het_reads` | Minimum reads on *each* allele for a het call. |
| `min_vaf_thres` | Het VAF must lie in `[min_vaf_thres, 1 - min_vaf_thres]`. |
| `filter_nz_OTH` | Drop SNPs with non-zero OTH (non-ref, non-alt) counts. |
| `filter_hom_ALT` | Drop hom-ALT SNPs. |

The cutoff is applied wherever **any** genotyped dataset is a tumor, one GT rule for the
whole run. When every genotyped dataset is `normal`, bulk skips this section entirely:
`snp_dir/chr{chrname}.vcf.gz` symlinks to `snp_dir/raw/`, so bcftools' germline GT and its
hom-alt calls pass through under the `params_bcftools` gates alone, and every field above,
the two filters included, is inert. Single-cell always applies the cutoff, cellsnp-lite
emitting no GT of its own.

There is no config key naming the mode. `genotype_dataset_ids` and each record's
`sample_type` are its only inputs, and the resolved value is logged at DAG build
(`tumor_genotyping_mode=...`). Left unset, `genotype_dataset_ids` auto-picks the first
`normal`, so a run with a matched normal keeps passing calls through; a sample with no
normal falls through to its tumor, which the cutoff then genotypes.

`panel_allele_only` follows the same input: genotyping a tumor forces it on, so a somatic
allele at a panel position can never become the called ALT. Two consequences of fixing the
alleles: `bcftools` aborts the chromosome when the panel's `REF` disagrees with `reference`,
and a site whose real ALT differs from the panel's is called hom-ref and drops out.

Depth in both sources is `REF + ALT`, with other-allele reads held apart in `OTH`, matching
cellsnp-lite (`src/csp.h`: `DP` is "total counts for ALT and REF"). A bulk site's
`FORMAT/DP` counts every base, so it is not the depth these thresholds see.

#### `params_longphase`
Used by `phase_snps_longphase`.

| Field | Description |
|---|---|
| `min_mapq` | Skip alignments below this mapping quality. |
| `extra_params` | Read chemistry: `--pb` (PacBio) \| `--ont` (ONT). |

#### `params_process_anndata`
Used by `process_rna_anndata` (single-cell RNA / spatial).

| Field | Description |
|---|---|
| `gene_id_colname` | Gene-id column joined against the GTF. |
| `min_frac_barcodes` | Drop a gene expressed in fewer than this fraction of barcodes. |

#### `params_phase_and_concat`
Used by `phase_and_concat_{bulk,nonbulk}`.

| Field | Description |
|---|---|
| `min_depth` | Minimum depth in every sample to keep a SNP (bulk). |
| `gamma` | Credible-interval level of the balanced-het test on the normal; a SNP is kept when its beta posterior interval covers 0.5 (bulk). |
| `exon_only` | Keep exonic SNPs only. |

#### `params_count_reads`
Used by `run_mosdepth` (per-base depth), `count_read_starts_chrom` (read-start counts) and `rd_correct` (HMMcopy-style bias correction), all bulk.

| Field | Description |
|---|---|
| `read_quality` | Skip alignments below this mapping quality. Reaches both `mosdepth -Q` and `samtools view -q`. |
| `exclude_flags` | SAM FLAG bitmask to exclude. Reaches both `mosdepth -F` and `samtools view -F`. |
| `mosdepth_extra_params` | Extra mosdepth flags. |
| `rd_correct_method` | `lowess` \| `median`. |
| `gc_correct` | Model GC content as a covariate during fitting. |
| `rt_correct` | Model replication-timing as a covariate during fitting. |
| `samplesize` | Number of sampled bins during fitting. |
| `routlier` | Upper quantile of read count dropped as outlier. |
| `doutlier` | Top/bottom quantile of the GC/mappability domain dropped as outlier. |
| `min_mappability` | Drop bins with mappability below this cutoff. |

> [!NOTE]
> `read_quality` and `exclude_flags` are deliberately shared so `window.dp.npz` and `{dataset_id}.rdcount.bed.gz` see the same reads.
>
> The default `exclude_flags: 1796` is mosdepth's own (unmapped, secondary, QC-fail, duplicate). It does **not** exclude supplementary alignments, each of which carries its own start position, so a split read is counted once per segment. The inflation scales with read length and concentrates at SV breakpoints. Long-read runs should set `exclude_flags: 3844` (`1796 | 2048`), which moves both files together.

> [!NOTE]
> With `target_bed` set and a `bulkWES` dataset in the run, `rd_correct` fits the bias
> correction on- and off-target apart for that dataset, and applies one library-size factor
> per (dataset, capture group) to the window depth of **every** bulk dataset. Why, and how
> to set it: [Whole-exome sequencing (WES)](bulk_genotyping.md#whole-exome-sequencing-wes).
>
> `{dataset_id}.rdcount.bed.gz` and `bb.rdcount.npz` are raw read counts: neither the bias
> correction nor this rescale touches them.

#### `params_combine_counts`
Used by `combine_counts` (bulk) and `combine_counts_nonbulk` (single-cell).

| Field | Description |
|---|---|
| `min_snp_reads` | SNP-covering reads every tumor column needs to close a bb; a list sweeps `MSR{msr}/`. |
| `min_snp_per_bin` | SNPs needed to close a bb. |
| `min_total_reads` | Read starts every column needs to close a bb (bulk); `0` disables. |
| `detect_loh_tumor_cell_line` | Call clonal-LOH regions from het-SNP density and bin them without the SNP criterion (bulk). Set it for a tumor with no normal cells, e.g. a cell line: a germline het inside LOH is called hom, so the density collapses there. |
| `loh_tile_size` | Unit of the clonal-LOH density chain (bp). A window is too sparse to count on; a tile must hold enough hets for the two rates to separate. |
| `loh_rate_ratio` | LOH-state het rate as a fraction of the fitted neutral rate. |
| `loh_breakpoint_rate` | Poisson breakpoints per bp for that chain; its reciprocal is the mean segment length (`1e-8` -> 100 Mb). The only thing resisting a one-tile flip, so it sets the reported region count. |
| `gene_aware_binning` | Grow bbs by whole genes; never cut inside one. |
| `nu` | Haldane scale turning cM distance into a switch probability. |
| `min_switchprob` | Floor on that switch probability. |
| `switchprob_ps` | Switch probability within one phase set (`PS`); ~0.5 across sets. |
| `nsnp_multi` | SNPs per multi-SNP diagnostic group. |
| `rdr_normalization` | Bulk RDR denominator: `auto` (base else median), `median`, `normal` (base required). |
| `phase_flip_test` | Split a phase cluster failing the haplotype-flip test (bulk). |
| `phase_flip_epsilon` | Effect size of that test (bulk). |
| `phase_flip_alpha` | Significance level of that test (bulk). |

> [!NOTE]
> Adaptive binning merges consecutive windows left to right and closes a bb when BOTH
> hold, the rule HATCHet2 uses (`adaptive_bins_arm`,
> [combine_counts.py](https://github.com/raphael-group/hatchet/blob/master/src/hatchet/utils/combine_counts.py)):
> - every tumor column has `min_snp_reads` SNP reads and the bb holds `min_snp_per_bin` SNPs;
> - `[bulk]` every column has `min_total_reads` read starts.
>
> There is no span cap: a cap can only fire by cutting a bb that has not met these, so
> the two criteria would contradict each other. What bounds a bb instead:
> - a bb never spans two `region_id`, `seg_id`, `loh_id` or `PS` clusters;
> - `[bulk]` nor two `phase_cluster` clusters, under `phase_flip_test`;
> - the next window starts a new gene, under `gene_aware_binning`;
> - a trailing run below threshold merges into the previous bb, and a cluster that never
>   meets the thresholds stays one bb (logged).
>
> Inside a clonal-LOH region the SNP criterion is dropped, since no germline het survives
> there; `min_total_reads` alone sizes those bbs and their allele counts are stamped.
>
> bbs are then post-filtered:
> - `[bulk]` drop a bb whose BAF is NaN: no SNP, or a column with no read over them.
> - `[bulk]` drop a bb whose depth is NaN: every window below `min_mappability`, or
>   (`lowess`) outside the fitted range; or a multi-SNP group overlapping no window.
> - `[bulk]` drop a bb whose RDR is NaN: any depth NaN, a zero or NaN base column, or an
>   all-NaN tumor column.
> - `[nonbulk]` no filter; SNP-free bbs are kept, with all-zero allele rows and their
>   `Xcount`.

#### `threads`
Used by all multi-thread rules.

| Field | Description |
|---|---|
| `genotype` | Threads for genotyping; `bcftools mpileup` is single-threaded, so these size the `call` and `view` output compressors only. |
| `phase` | Threads for phasing. |
| `pileup` | Threads for the pileup step (bulk: the `bgzip` writing the counts; single-cell: cellsnp-lite, which is genuinely parallel). |
| `mosdepth` | Threads for the bulk read-counting rules: `mosdepth -t`, and `samtools view -@` in `count_read_starts_chrom`. |

> [!NOTE]
> `bcftools` applies `--threads` to the output handle only, never to the BAM/CRAM
> readers, so a bulk `mpileup` runs on one core whatever this is set to. Bulk genotyping
> and pileup both fan out per chromosome instead, and small values here leave cores free
> for more concurrent chromosome jobs.

---

## Outputs

Set in `config.yaml`, relative to `snakemake --directory`:

| Key | Path | Contents |
|---|---|---|
| `snp_dir` | `snps` | Genotyped SNP VCFs. |
| `phase_dir` | `phase` | Phased VCFs, parsed genetic map. |
| `pileup_dir` | `pileup` | Per-dataset allele counts and read depth. |
| `allele_dir` | `allele` | The SNP grid and its allele matrices. |
| `bb_dir` | `bb` | The bbs and their matrices. |
| `qc_dir` | `qc` | One multi-page PDF per rule. |
| `log_dir` | `logs` | One log per job. A rule that fans out (per chromosome, dataset or assay) groups its jobs under `{stage}/`; a rule that runs once writes `{rule}.{run_id}.log` at the top level. |
| `aux_dir` | `aux` | Segment BED, window BED, Repli-seq tracks. |
| `bench_dir` | `benchmarks` | Runtime and `max_rss` per job, laid out exactly like `log_dir`. |

### Genomic unit levels

| Level | id | Description |
|---|---|---|
| Region | `region_id` | Chromosome arm, from `region_bed` (`chr1p`). Carried for RDR and QC. |
| Segment | `seg_id` | `{region_id}#{START}-{END}`. |
| Window | `bin_id` | `window_size` tile, shared by every assay. Internal use. |
| bb | `bb_id` | Merged windows; the feature axis of `bb.tsv.gz` and `bb.*.npz`. |

### Final bins

A `min_snp_reads` list writes one `MSR{msr}/` per value.

| Mode | Location |
|---|---|
| `bulk_genotyping` | `bb_dir/MSR{msr}/bulk/` |
| `single_cell_genotyping` | `bb_dir/MSR{msr}/{assay_type}/` |
| `copytyping_preprocess` | `bb_dir/{assay_type}/` |

**`bulk_genotyping`** - one bb set over every bulk assay; columns are samples.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.{depth,rdr}.npz` | Depth, and RDR for the tumor columns only. |
| `bb.rdcount.npz` | Read starts summed from the windows, int32, every column. |
| `sample_ids.tsv` | One row per matrix column. |
| `multi_snp/bulk/` | Multi-SNP diagnostic groups, bb schema, outside `MSR{msr}/`. |

**`single_cell_genotyping`** - one bb set over every assay, copied into each subdir;
columns are cells.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.Xcount.npz` | Native counts: scATAC fragment midpoints, scRNA/VISIUM h5ad UMIs. |
| `sample_ids.tsv` | One row per dataset, not column-aligned. |
| `barcodes.tsv.gz` | The matrix column axis; each barcode is suffixed `_{dataset_id}_{assay_type}`. |
| `multi_snp/{assay_type}/` | Multi-SNP diagnostic groups, bb schema, outside `MSR{msr}/`. |

**`copytyping_preprocess`** - the given bb set, per assay; columns are cells.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis; the given `bb_file` plus `#SNPS`, `feature_id`, (RNA) `#feature`. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.Xcount.npz` | Native counts: scATAC fragment midpoints, scRNA/VISIUM h5ad UMIs. |
| `sample_ids.tsv` | One row per dataset, not column-aligned. |
| `barcodes.tsv.gz` | The matrix column axis; each barcode is suffixed `_{dataset_id}_{assay_type}`. |

### Unit level

`bb_dir/unit/` holds the grids the binning consumes, un-binned and independent of
`min_snp_reads`. Row axes are TSVs and matrices `.npz`, laid out as in `MSR{msr}/`.

**`bulk_genotyping`** - `bb_dir/unit/bulk/`; columns are samples.

| File | Contents |
|---|---|
| `snp.tsv.gz` | The SNPs that landed in a window, the allele matrix row axis. |
| `snp.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `window.tsv.gz` | The windows on the run's chromosomes, the depth matrix row axis. |
| `window.depth.npz` | Bias-corrected depth, windows x datasets. |
| `window.rdcount.npz` | Read starts, windows x datasets, int32. |
| `sample_ids.tsv` | One row per matrix column. |

**`single_cell_genotyping`** and **`copytyping_preprocess`** - `bb_dir/unit/{assay_type}/`; columns are cells.

| File | Contents |
|---|---|
| `snp.tsv.gz` | The SNPs that landed in a window, duplicated into each subdir. |
| `snp.{Tallele,Aallele,Ballele}.npz` | That assay's column slice of the allele counts. |
| `barcodes.tsv.gz` | The matrix column axis. |
| `sample_ids.tsv` | One row per dataset, not column-aligned. |
| `window.{tsv.gz,Xcount.npz}` | scATAC only: fragments counted per window per cell. |
| `gene.{tsv.gz,Xcount.npz}` | RNA assays only: UMIs per gene per cell. A gene is never split, so the gene is the RNA unit and no window matrix is written. |

`copytyping_preprocess` bins onto the given `bb_file`, but tiles a window grid all the
same: in every mode the windows carry no blacklisted span, and scATAC fragments and
SNPs are counted through them, so a blacklisted fragment or SNP inside a bb's span is
dropped rather than swallowed by the hull. A gene is indivisible and is still assigned
whole to the bb it overlaps most.

### Intermediates

| Path | Contents |
|---|---|
| `snp_dir/raw/chr{chrname}.vcf.gz` | Bulk: the caller's own output, before `post_genotype_snps_bulk`. |
| `snp_dir/chr{chrname}.vcf.gz` | Bi-allelic SNPs, what phasing reads. |
| `snp_dir/pseudobulk_{modality}/` | cellsnp-lite pseudobulk output. |
| `phase_dir/chr{chrname}.vcf.gz` | Phased SNPs, concatenated to `phased_het_snps.vcf.gz`. |
| `phase_dir/genetic_map.tsv.gz` | Parsed genetic map (eagle/shapeit). |
| `pileup_dir/{assay_type}_{dataset_id}/` | Bulk `bcftools.counts.tsv.gz`, concatenated from per-chromosome `bcftools.counts.chr{chrname}.tsv.gz` (temporary); single-cell `cellSNP.*`. |
| `pileup_dir/{assay_type}/out_mosdepth/` | Per-dataset mosdepth (bulk). |
| `pileup_dir/{assay_type}/out_rdcount/` | Per-dataset, per-chromosome read-start counts (bulk, temporary). |
| `pileup_dir/{assay_type}/{dataset_id}.rdcount.bed.gz` | Per-window read-start counts (bulk): `#CHR START END COUNT`, headerless, mosdepth's `regions.bed.gz` shape. `combine_counts` joins it onto the window grid. |
| `pileup_dir/bulk/window.raw.dp.npz` | Raw read depth, windows x every bulk dataset. |
| `pileup_dir/bulk/window.dp.npz` | Corrected read depth, windows x every bulk dataset. |
| `allele_dir/` | `snps.tsv.gz`, `snp.{T,A,B}allele.npz`, `sample_ids.tsv`, `barcodes.tsv.gz`. |
| `bb_dir/{assay_type}.h5ad` | Gene x cell AnnData (scRNA / spatial). |
| `aux_dir/segment.bed` | Arms cut at the SV extremities, blacklist subtracted. |
| `aux_dir/windows.bed.gz` | The shared window BED. |
| `aux_dir/repliseq/` | Repli-seq tracks, lifted from hg19 when the run is not hg19. |

### TSV columns

| File | Columns |
|---|---|
| `snps.tsv.gz` | `#CHR POS POS0 START END GT PHASE region_id seg_id feature_id feature_type`; bulk adds `PS`. |
| `bb.tsv.gz` | `#CHR START END #SNPS region_id switchprobs feature_id`, plus `#feature` for an RNA assay and `is_loh` under `detect_loh_tumor_cell_line`: those rows carry a stamped BAF, not a measured one. |
| `snp.tsv.gz` | The `snps.tsv.gz` columns, restricted to the SNPs inside a window. |
| `window.tsv.gz` | `#CHR START END region_id seg_id`, plus `is_loh` under `detect_loh_tumor_cell_line`. |
| `gene.tsv.gz` | `#CHR START END feature_id region_id`. |
| `sample_ids.tsv` | `SAMPLE sample_id dataset_id sample_type assay_type`; bulk adds `rdr_base_dataset_id`. |
| `germline_snp_statistics.tsv` | Per chromosome: het_phased, het_unphased, hom_alt, hom_ref. |
| `depth_statistics.tsv` | Per-dataset depth summary, every bulk dataset in one table. |

> [!NOTE]
> - `PHASE`: 0 = the B-allele is ALT, 1 = the B-allele is REF.
> - `feature_id`: `;`-joined overlapping GTF genes, `intergenic` if none.
> - `sample_ids.tsv` is column-aligned in bulk only; single-cell columns are cells.
> - `SAMPLE` is `{sample_id}_{dataset_id}`, plus `_{assay_type}` for a multiome pair.

### QC (`qc_dir/`)

One multi-page PDF per rule, flat:

| File | Contents |
|---|---|
| `post_genotype_snps.{bulk_or_nonbulk}.pdf` | SNP allele frequency by genotype. |
| `detect_loh.bulk.pdf` | Het-SNP density per tile against the two fitted rates, coloured by the decoded state (`detect_loh_tumor_cell_line` only). |
| `phase_and_concat.{bulk_or_assay}.pdf` | SNP allele frequency and depth. |
| `rd_correction.bulk.pdf` | Depth before/after correction, GC/MAP/RT diagnostics. |
| `combine_counts.{bulk_or_assay}.MSR{msr}.pdf` | Binning QC, one per `min_snp_reads`. |
| `combine_counts_fixed_bins.{assay_type}.pdf` | Per dataset: a SNP BAF page, then a bb page of pseudobulk RDR over BAF. |
