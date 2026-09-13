# TODO

# MSR vector to support joint whole-genome / targeted segmentation

Replace the scalar `min_snp_reads` with a per-dataset vector derived from a target BAF standard
error, so assays of unequal read supply (WGS + WES, WGS + CRISPR-targeted) can share one bb grid.
Design, model, and implementation steps: `.claude/msr_vector.md`.

# slurm support

## Test pipeline

Goal: after a user clones + installs the pipeline, one command runs a real end-to-end
test that exercises a full mode, not a per-rule dry run.

Two layers:
- **Dry-run (DAG) tests** (`tests/test_dryrun.py`, `tests/test_remote.py`): built, all
  modes x JSON/TSV; execute no rule, need no data. This is what CI runs.
- **End-to-end tests**: a real `snakemake` run to final `bb_dir/` outputs. Not in the
  repo. Sample sheets and configs for real cases live outside it, since the inputs are
  whole-genome remote alignments and no hosted runner can execute them.

- [ ] Decide where end-to-end case definitions live (separate repo, or a gitignored
      working dir) and how references are staged.
- [ ] Cases to cover: `single_cell_genotyping`, `copytyping_preprocess`, spatial.

Harness: if end-to-end returns, gate it behind a manual/scheduled trigger and keep the
dry-run tests on every push.

## Tumor-only genotyping (no matched normal)

Today every mode needs a normal to call germline SNPs: `genotype_dataset_ids` defaults to the
first normal, and a tumor-only run genotypes off tumor reads with a WARN.

- [ ] Distinguish germline het from hom-alt SNPs in a high-purity tumor without a matched
      normal. Adapt https://github.com/raphael-group/hetdetect.
    - retrieve population ALT frequency as prior genotype info. high ALT freq indicates likely hom-alt

## Segmentation parameter selection

`min_snp_reads` (and `min_total_reads`) are chosen by hand: a list sweeps the grid and the
user picks a point off the QC PDFs. The `select_segmentation` script that scored a
`(min_snp_reads, max_blocksize)` grid is not in the tree, and `max_blocksize` itself is gone.

- [ ] Automated model selection for the segmentation parameters against sequencing coverage
      and segmentation variance, i.e. recommend one grid point instead of a sweep.
- [ ] Concretely: pick a default MSR at the elbow of lag-1 RDR/BAF dispersion vs bin count,
      and record the pick alongside the outputs (was a TODO comment in `combine_counts.py`,
      pointing at a `docs/combine_counts_pseudocode.md` that is not in the repo).

## Others
- Streaming remote data (DONE for bulk): `remote_mode: stream` reads remote BAM/CRAM directly
  with bcftools/mosdepth, fetching only the config `chromosomes` (index jumps); default stays
  `storage` (whole-file download). Single-cell/copytyping cannot stream (`cellsnp-lite` rejects
  URLs via its `access(F_OK)` guard). Follow-up: validate `##idx##` remote-index support and
  numeric parity on a real URL BAM (see plan verification).

## Skip phasing + within-bb EM phasing

Two linked additions. Full design, integration points, tests and validation:
`.claude/skip_phasing_and_bb_em.md`.

1. `phaser: "none"` - genotype, skip panel and long-read phasing, carry the unphased het SNPs
   into pileup. For runs where no phasing prior exists: no panel for the build, no long reads,
   no external phased VCF.
2. `params_combine_counts.phase_correction` (`none|flip_split|bin_em`) - re-orient SNPs within
   each bb after adaptive binning, before allele counts are summed per bb. Subsumes the
   `phase_flip_test` boolean (`true` -> `flip_split`, `false` -> `none`).

Bulk `combine_counts` collapses CNA/LOH BAF toward 0.5 (seen on `hatchet2_chr22_simulation`
dbSNP151; panel/phaser ruled out). Two causes: (1) `detect_phase_flips` fragments SNPs into
thousands of phase-groups, and the per-group bin-count floor in `_merge_bins_to_bbs` makes
~4-SNP bins so `min_snp_reads` never binds; (2) `apply_phase_to_mat` orients A/B by one
per-SNP `PHASE` bit with no within-bin re-orientation, so unfolded bin BAF is exactly 0.500.
HATCHet2 gets ~220 SNPs/bin and BAF ~0.20 in LOH via a per-bin EM.

Model: naive Bayes EM over the tumor samples, in a new helper
`workflow/scripts/script_utils/phase_em.py`. The per-SNP phase latent is independent across
position; the coupling is across samples. Ports HATCHet2 `multisample_em`
(`hatchet/utils/combine_counts.py`, doi:10.1038/s41587-020-0661-6); same shape as Alleloscope
(doi:10.1038/s41587-021-00911-w).

> [!IMPORTANT]
> The earlier HMM design (`phase_hmm.py`, beta-binomial emission, LD switch/stay transitions,
> tau calibrated on the matched normal, optional `cross_bin_phasing: dp` layer;
> `.claude/phase-hmm-within-bin.md`) is superseded. Its problem statement above still
> holds; its model does not.

- [ ] Phase A: `phaser: "none"` (const, the `phased_snp_vcf` resolution in `parse_workflow`,
      `phase_snps.smk`);
      `phase_correction` config surface; `correct_bin_phases` as a documented no-op wired into
      bulk and single-cell `combine_counts`; the `estimate_switchprobs_PS` `KeyError: 'PS'` fix;
      dry-run tests; docs.
- [ ] Phase B: implement the EM; A/B comparison vs HATCHet2 `bb`; make `bin_em` default.

## RD bias correction (potential over-correction)

Observed on `hatchet2_chr22_simulation` (chr22-only, matched normal + 3 tumors): GC/RT
correction lowers the GC/RT correlation but raises RD dispersion (T2 GC MAD 10.80 -> 11.58),
i.e. removes weak bias while injecting the covariate shape (the RT replication wave).
Contributing factors: weak biases (GC r ~ 0.12, RT r ~ 0.1-0.24), per-sample quadratic fit
over CNA-contaminated bins (`rd_correct_utils.py`; highest-CNA sample hurt most), the
`RT + RT**2` term imprinting the replication-timing wave, and a single-chromosome fit.

- [ ] `gc_correct`, `rt_correct`: options `[true, false, auto]`, resolved per dataset, not
      one global switch.
- [ ] `auto`: choose each covariate from the data - explained-variance / correlation cutoff
      per covariate, matched-normal presence (tumor/normal RDR already cancels shared GC/RT
      bias), min bin count / genome coverage - and drop a covariate below threshold.
- [ ] Fit the correction on CN-neutral bins, or derive it from the normal and apply to the
      tumors, instead of a per-tumor fit over CNA-altered bins.
- [ ] QC gate: flag or revert a correction when it increases a sample's RD dispersion
      (MAD after > MAD before).
- [ ] Fit genome-wide (full GC/RT range, enough bins); warn when a single chromosome or low
      bin count makes the quadratic fit unreliable.
- [ ] Lower-order or regularized model when a covariate's bias is weak.

## Unified bulk grid: WES treated exactly like WGS (IMPLEMENTED)

Done: `build_segment_bed` (region_id arm + seg_id chunk) + one shared window BED
(`aux/windows.bed.gz`, tiled from `segment.bed` at `window_size`). Every bulk assay
(WGS/WGS-lr/WES) bins on that one grid grouped by `seg_id`, with a single read target
`min_snp_reads` gated behind it; one `bb_dir/MSR{msr}/bulk/`. The
earlier WES stream (per-record `wes_targets_bed`, 267 bp exon tiling, WES-onto-WGS depth
projection, `min_snp_reads_wes`) was removed. Verified by DAG tests only. Follow-ups:

- Real-data validation: no execution-level test exists yet; run a real WGS+WES bulk sample
  end-to-end and inspect per-bin WES RDR/BAF.
- WES RD-correction: WES depth carries capture-enrichment structure; confirm the per-sample
  LOWESS fit + `routlier`/`doutlier` handle it, or flag over-correction (see the RD bias
  correction section above).

## streaming panel
`https://ftp.ncbi.nih.gov/snp/organisms/`, `tabix-streams just that region from the remote file — never downloading the ~16 GB whole thing (genotype_snps.py:232): bcftools query -f '%CHROM\t%POS\n' -r chr22 <URL> -> target_chr22.pos.gz.`
what is usual time cost for tabix snp positions from local/online snp panel -> genotyping?